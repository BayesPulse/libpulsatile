#!/usr/bin/env Rscript
# vignettes/precompute.R
#
# Precompute (".Rmd.orig") pattern. Each NN-name.Rmd.orig contains real,
# possibly-slow MCMC. This script knits each one to a static NN-name.Rmd whose
# code blocks are NOT re-executed by R CMD build (they are plain ```r fences,
# not ```{r} chunks), so the package builds in seconds while readers still see
# genuine converged results.
#
# Run from the repo root:
#   Rscript vignettes/precompute.R          # knit all *.Rmd.orig
#   Rscript vignettes/precompute.R getting   # knit only matching files
#
# Commit the resulting *.Rmd and *.png files alongside the *.Rmd.orig source.
#
# Reproducibility note: the committed .Rmd files bake in numeric output from
# one specific R/BLAS/RNG environment. Re-knitting on the same machine is
# deterministic (seeded), but a different R version or BLAS can shift figures
# slightly. If a dependency change alters results, RE-RUN this script -- never
# hand-edit numbers in the generated .Rmd files.

if (!requireNamespace("knitr", quietly = TRUE)) stop("install 'knitr'")
library(bayespulse)

vdir  <- "vignettes"
origs <- list.files(vdir, pattern = "\\.Rmd\\.orig$")
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) origs <- origs[grepl(args[1], origs, fixed = TRUE)]
if (length(origs) == 0) { message("No matching .Rmd.orig files."); quit(save = "no") }

local({
  old <- setwd(vdir)
  on.exit(setwd(old))
  for (o in origs) {
    out <- sub("\\.orig$", "", o)
    message("Knitting ", o, " -> ", out)
    knitr::knit(input = o, output = out)
  }
})
message("Done. Review the generated .Rmd + figure files, then git add them.")
