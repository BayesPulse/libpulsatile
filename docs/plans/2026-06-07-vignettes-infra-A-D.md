# bayespulse Vignettes: Infrastructure + Vignettes A–D — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Wire vignettes into the R package build and ship the first four ladder vignettes (A getting-started, B anatomy-of-the-model, C priors, D convergence) from the approved roadmap (`docs/plans/2026-06-07-vignette-roadmap-design.md`).

**Architecture:** Vignette sources live in repo-root `vignettes/` (source of truth), exposed to the package via a `R-package/vignettes -> ../vignettes` symlink, matching the repo's existing symlink convention (`R-package/inst/include -> ../../include`). Vignettes with real MCMC (A, D) use the **precompute / `.Rmd.orig` pattern**: a `*.Rmd.orig` runs slow MCMC offline and is knitted to a static `*.Rmd` that ships and builds in seconds with no computation. Vignettes without MCMC (B, C) are ordinary live `.Rmd` files. This keeps `R CMD build`/`check` fast despite ~6 ms/iter BD-MCMC, while still showing real converged results.

**Tech Stack:** R package tooling (`devtools`, `R CMD check`), `knitr` + `rmarkdown` (`html_vignette`), `ggplot2`, the `bayespulse` package itself.

**Key API facts (verified this session):**
- `simulate_pulse(...)` returns a `pulse_sim` object: `$data` (cols `observation`, `time`, `concentration`), `$parameters`, `$call`. `plot()` dispatches to `plot.pulse_sim`.
- `pulse_spec(...)` builds the single-subject spec (Strauss location prior etc.).
- `fit_pulse(data, spec, iters=250000, thin=50, burnin=0.1*iters, verbose=FALSE)` — pass `data = simdat$data` (a data frame) or a `pulse_sim`. Returns class `pulse_fit`.
- Single-subject plots: `bp_trace(fit)`, `bp_posteriors(fit, type=c("histogram","density"))`, `bp_location_posterior(fit)`; `predict(fit, cred_interval=0.8)` → data frame with `mean.conc`/`upper`/`lower`; `bp_predicted(fit, predicted)`.
- `patient_chain(fit)` → per-iteration population-equivalent params (`baseline`, `halflife`, `mass_mean`, `width_mean`, ...). `pulse_chain(fit)` → per-pulse draws.
- **Diagnostics that are population-only** (`stop()` on non-`population_fit`): `population_diagnostics`, `plot_trace`, `plot_acf`, `convergence_report`. **Model-agnostic:** `effective_sample_size(x)` (numeric vector), `gelman_rubin(chains)` (list of numeric vectors). `bp_trace(fit)` works on single-subject. → Vignette D uses only the model-agnostic tools + `bp_trace`, and forward-references the population suite to vignette E (decision recorded 2026-06-07).

**Precompute is slow but offline & one-time.** A ≈ 6 min; D ≈ 4 chains × ~4 min ≈ 16 min. `R CMD build` runs **no** MCMC.

---

## Phase 0: Vignette build infrastructure

### Task 0.1: Declare the vignette builder in DESCRIPTION

**Files:**
- Modify: `R-package/DESCRIPTION`

**Step 1: Add knitr/rmarkdown to Suggests and add VignetteBuilder.**

Current `Suggests:` block is:
```
Suggests:
  magrittr,
  RcppArmadillo,
```
Change it to (and add the `VignetteBuilder` line immediately after the `Suggests` block):
```
Suggests:
  magrittr,
  RcppArmadillo,
  knitr,
  rmarkdown
VignetteBuilder: knitr
```
(`ggplot2` is already in `Imports:`, so vignette plotting needs no new dependency.)

**Step 2: Verify the field parses.**

Run: `Rscript -e 'd <- read.dcf("R-package/DESCRIPTION"); cat(d[,"VignetteBuilder"], "\n"); cat(d[,"Suggests"], "\n")'`
Expected: prints `knitr` for VignetteBuilder, and a Suggests string containing `knitr` and `rmarkdown`.

**Step 3: Commit.**
```bash
git add R-package/DESCRIPTION
git commit -m "build: declare knitr vignette builder in DESCRIPTION"
```

---

### Task 0.2: Symlink the package vignettes directory

**Files:**
- Create symlink: `R-package/vignettes -> ../vignettes`

**Step 1: Create the symlink (matches repo convention).**
```bash
ln -s ../vignettes R-package/vignettes
```

**Step 2: Verify it resolves to the root vignettes dir.**

Run: `ls -l R-package/vignettes && ls R-package/vignettes/`
Expected: symlink shown as `R-package/vignettes -> ../vignettes`; listing shows the current root `vignettes/` contents (`my-vignette.Rmd`, etc.).

**Step 3: Commit.**
```bash
git add R-package/vignettes
git commit -m "build: symlink R-package/vignettes to root vignettes (repo convention)"
```

---

### Task 0.3: Configure ignore rules for precompute artifacts

**Files:**
- Modify: `R-package/.Rbuildignore`
- Modify: `vignettes/.gitignore`

**Step 1: Exclude dev-only precompute files from the built package.**

Append to `R-package/.Rbuildignore`:
```
^vignettes/.*\.Rmd\.orig$
^vignettes/precompute\.R$
^vignettes/my-vignette_files$
```
(The first two keep `.Rmd.orig` sources and the precompute driver out of the tarball; the third drops the stale `my-vignette_files/` artifact directory.)

**Step 2: Let `precompute.R` be tracked despite the `*.R` ignore.**

Current `vignettes/.gitignore`:
```
*.html
*.R
```
Change to:
```
*.html
*.R
!precompute.R
```
This keeps generated `.R` (e.g. knitr `purl` output) ignored but tracks our driver script. Generated `*.Rmd` and `*.png` figures are NOT ignored and will be committed.

**Step 3: Commit.**
```bash
git add R-package/.Rbuildignore vignettes/.gitignore
git commit -m "build: ignore precompute sources in build, track precompute.R in git"
```

---

### Task 0.4: Add the precompute driver script

**Files:**
- Create: `vignettes/precompute.R`

**Step 1: Write the driver.**
```r
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

if (!requireNamespace("knitr", quietly = TRUE)) stop("install 'knitr'")
library(bayespulse)

vdir  <- "vignettes"
origs <- list.files(vdir, pattern = "\\.Rmd\\.orig$")
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) origs <- origs[grepl(args[1], origs)]
if (length(origs) == 0) { message("No matching .Rmd.orig files."); quit(save = "no") }

old <- setwd(vdir)
on.exit(setwd(old))
for (o in origs) {
  out <- sub("\\.orig$", "", o)
  message("Knitting ", o, " -> ", out)
  knitr::knit(input = o, output = out)
}
message("Done. Review the generated .Rmd + figure files, then git add them.")
```

**Step 2: Verify it runs as a no-op when there are no `.Rmd.orig` files yet.**

Run: `Rscript vignettes/precompute.R`
Expected: prints `No matching .Rmd.orig files.` and exits 0.

**Step 3: Commit.**
```bash
git add vignettes/precompute.R
git commit -m "build: add vignette precompute driver script"
```

---

### Task 0.5: Smoke-test the vignette pipeline end-to-end

**Files:**
- Create (temporary): `vignettes/_smoke.Rmd`

**Step 1: Write a trivial live vignette (no MCMC) to prove the wiring builds fast.**
````markdown
---
title: "smoke"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{smoke}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r}
1 + 1
```
````

**Step 2: Build vignettes and confirm the engine is wired.**

Run: `cd R-package && Rscript -e 'devtools::build_vignettes()'`
Expected: completes without error; reports building `_smoke.Rmd`; produces output under `R-package/doc/` (or `doc/`). If it errors with "no vignette engine", revisit Task 0.1/0.2.

**Step 3: Delete the smoke vignette and any generated `doc/` output.**
```bash
rm vignettes/_smoke.Rmd
rm -rf R-package/doc R-package/Meta
```

**Step 4: Commit (infra verified; nothing to ship from the smoke test).**
```bash
git add -A vignettes
git commit -m "build: verify vignette pipeline (smoke test, removed)" --allow-empty
```

> Note: with the smoke vignette removed, `R CMD check` may emit a transient NOTE about a declared `VignetteBuilder` with no vignettes. Phase 1 adds the first real vignette and clears it.

---

## Phase 1: Vignette A — Getting started (single-subject), folding in `my-vignette.Rmd`

Uses the precompute pattern (real MCMC). Reading-order slug: `01-getting-started`.

### Task 1.1: Write `01-getting-started.Rmd.orig`

**Files:**
- Create: `vignettes/01-getting-started.Rmd.orig`

**Step 1: Write the source.** Salvage the intro/Section-2 math from `my-vignette.Rmd`; drop all of Rachel's TODO/notes prose. Full source:

````markdown
---
title: "Getting started: single-subject pulse deconvolution"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting started: single-subject pulse deconvolution}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.path = "01-getting-started-",
                      fig.width = 7, fig.height = 5, dpi = 96)
library(bayespulse)
set.seed(2026)
```

<!-- PROSE: 2-3 short paragraphs. What pulsatile hormones are; that bayespulse
     deconvolves an observed series into baseline + pulses + elimination + error;
     that this vignette walks one subject end to end. Keep the deconvolution
     equation from my-vignette.Rmd Section 2 if desired, trimmed. -->

## 1. Simulate a hormone series

```{r simulate}
simdat <- simulate_pulse(num_obs = 60, interval = 10, error_var = 0.005,
                         ipi_mean = 12, ipi_var = 10, ipi_min = 4,
                         mass_mean = 2.5, mass_sd = 0.5,
                         width_mean = 40, width_sd = 5,
                         constant_halflife = 40, constant_baseline = 2.5)
```

```{r plot-sim}
plot(simdat)
```

## 2. Specify the model

<!-- PROSE: 1 paragraph: priors encode physiological knowledge; the Strauss
     location prior discourages implausibly-close pulses. Full prior treatment
     lives in vignette C. -->

```{r spec}
spec <- pulse_spec(
  location_prior_type   = "strauss",
  prior_location_gamma  = 0.1,  prior_location_range = 30,
  prior_halflife_mean   = 38,   prior_halflife_var   = 1000,
  prior_sd_mass         = 0.5,  prior_sd_width       = 5,
  prior_mass_mean       = 0.5,  prior_mass_var       = 2,
  prior_mean_pulse_count = 10,  prior_width_mean     = 35,
  prior_baseline_mean   = 2.25, prior_baseline_var   = 100,
  sv_mass_mean = 2, sv_width_sd = 10, sv_error_var = 0.005, sv_mass_sd = 0.4,
  sv_baseline_mean = 2, sv_halflife_mean = 42,
  pv_mean_pulse_width = 10, pv_indiv_pulse_width = 10, pv_sd_pulse_width = 0.1,
  pv_pulse_location = 5, pv_sd_pulse_mass = 0.45)
```

## 3. Fit the model

<!-- PROSE: 1 paragraph: birth-death MCMC; we run a production-length chain.
     Note this chunk is precomputed offline (see vignettes/precompute.R). -->

```{r fit}
fit <- fit_pulse(data = simdat$data, spec = spec,
                 iters = 60000, thin = 20, burnin = 10000, verbose = FALSE)
```

## 4. Inspect the results

```{r trace}
bp_trace(fit)
```

```{r posteriors}
bp_posteriors(fit, type = "histogram")
```

```{r locations}
bp_location_posterior(fit)
```

```{r predicted}
predicted <- predict(fit, cred_interval = 0.8)
bp_predicted(fit, predicted)
```

<!-- PROSE: 1-2 sentences each: trace = mixing/chains over iterations; posteriors
     = parameter uncertainty; location plot = where pulses occurred; predicted =
     fitted concentration with 80% credible band over the observed data. Close
     with a pointer: vignette D covers convergence; E covers populations. -->
````

**Step 2: Verify the source parses as R (chunks only).**

Run: `Rscript -e 'knitr::purl("vignettes/01-getting-started.Rmd.orig", output=tempfile(), quiet=TRUE); cat("parsed OK\n")'`
Expected: prints `parsed OK` with no error.

**Step 3: Commit the source.**
```bash
git add vignettes/01-getting-started.Rmd.orig
git commit -m "docs(vignette): add getting-started source (.Rmd.orig)"
```

### Task 1.2: Precompute vignette A (run the real MCMC offline)

**Step 1: Knit the source to a static `.Rmd` + figures.**

Run: `Rscript vignettes/precompute.R getting`
Expected: prints `Knitting 01-getting-started.Rmd.orig -> 01-getting-started.Rmd`; takes ~5–6 min; creates `vignettes/01-getting-started.Rmd` and `vignettes/01-getting-started-*.png`. No errors.

**Step 2: Verify the generated `.Rmd` has NO live chunks (so build won't re-run MCMC).**

Run: `grep -c '```{r' vignettes/01-getting-started.Rmd`
Expected: `0` (knitr converts executed chunks to plain ```` ```r ```` fences). If non-zero, the precompute did not run as expected — investigate before shipping.

**Step 3: Confirm the figure files exist.**

Run: `ls vignettes/01-getting-started-*.png`
Expected: several PNGs (sim plot, trace, posteriors, locations, predicted).

**Step 4: Commit generated artifacts.**
```bash
git add vignettes/01-getting-started.Rmd vignettes/01-getting-started-*.png
git commit -m "docs(vignette): precompute getting-started results"
```

### Task 1.3: Build vignette A through the package and verify it's fast

**Step 1: Build vignettes.**

Run: `cd R-package && time Rscript -e 'devtools::build_vignettes()'`
Expected: builds `01-getting-started.Rmd` in **seconds** (no MCMC); produces HTML output; no errors/warnings.

**Step 2: Retire the orphaned WIP vignette (now superseded).**
```bash
git rm vignettes/my-vignette.Rmd
git rm -r --cached vignettes/my-vignette_files 2>/dev/null || true
rm -rf vignettes/my-vignette_files R-package/doc R-package/Meta
```

**Step 3: Commit.**
```bash
git add -A vignettes R-package
git commit -m "docs(vignette): retire WIP my-vignette.Rmd (folded into getting-started)"
```

---

## Phase 2: Vignette B — Anatomy of the model (live, no MCMC)

No MCMC → ordinary live `.Rmd` (builds fast, no precompute). Slug: `02-anatomy-of-the-model`.

### Task 2.1: Write `02-anatomy-of-the-model.Rmd`

**Files:**
- Create: `vignettes/02-anatomy-of-the-model.Rmd`

**Step 1: Write the vignette.**
````markdown
---
title: "Anatomy of the model: the generative process"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Anatomy of the model: the generative process}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.5, dpi = 96)
library(bayespulse)
set.seed(2026)
```

<!-- PROSE: an observed series = baseline + pulsatile secretion (Gaussian
     pulses) convolved with exponential elimination + measurement error. This
     vignette builds intuition for each piece by simulating and varying it. -->

## A single simulated series

```{r base-series}
sim <- simulate_pulse(num_obs = 72, interval = 10,
                      mass_mean = 3.5, width_mean = 35,
                      constant_baseline = 2.6, constant_halflife = 45)
plot(sim)
```

## How pulse mass changes the series

```{r sweep-mass, fig.height = 7}
masses <- c(1.5, 3.5, 6.0)
op <- par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
for (m in masses) {
  s <- simulate_pulse(num_obs = 72, interval = 10, mass_mean = m,
                      width_mean = 35, constant_baseline = 2.6,
                      constant_halflife = 45)
  plot(s$data$time, s$data$concentration, type = "l",
       main = paste("mass_mean =", m), xlab = "Time (min)", ylab = "Conc.")
}
par(op)
```

<!-- Repeat the sweep idiom for width_mean, constant_halflife, and ipi_mean
     (inter-pulse interval), each in its own section with 1 sentence of
     interpretation: wider pulses = broader humps; longer half-life = slower
     decay/higher troughs; larger IPI = fewer pulses. -->
````
(Add the three additional sweep sections for `width_mean`, `constant_halflife`, and `ipi_mean` following the same `for`-loop idiom.)

**Step 2: Verify it builds (fast, live).**

Run: `cd R-package && Rscript -e 'devtools::build_vignettes()'` then `rm -rf R-package/doc R-package/Meta`
Expected: builds `02-anatomy-of-the-model.Rmd` in seconds, no errors.

**Step 3: Commit.**
```bash
git add vignettes/02-anatomy-of-the-model.Rmd
git commit -m "docs(vignette): add anatomy-of-the-model vignette"
```

---

## Phase 3: Vignette C — Priors & specification (live, prior-predictive only)

Prior-predictive simulation only (no MCMC) → live `.Rmd`. Slug: `03-priors-and-specification`.

### Task 3.1: Write `03-priors-and-specification.Rmd`

**Files:**
- Create: `vignettes/03-priors-and-specification.Rmd`

**Step 1: Write the vignette.** Structure:
- **Prose:** what each `pulse_spec` argument means physiologically — `prior_baseline_mean/var`, `prior_mass_mean/var`, `prior_width_mean`, `prior_halflife_mean/var`, `prior_mean_pulse_count`, and the Strauss location prior (`location_prior_type`, `prior_location_gamma`, `prior_location_range`). A short table mapping argument → meaning → "informative vs vague" guidance. Use Unicode escapes for Greek letters per CRAN rule (e.g. `μ`).
- **Code (prior-predictive checks):** draw parameters from a prior and simulate, to *see* what the prior implies. Example chunk:
```{r prior-predictive}
set.seed(2026)
library(bayespulse)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (i in 1:4) {
  bl <- rnorm(1, mean = 2.25, sd = sqrt(1.0))      # baseline prior draw
  hl <- rnorm(1, mean = 38,   sd = sqrt(50))       # half-life prior draw
  s  <- simulate_pulse(num_obs = 60, interval = 10, mass_mean = 3,
                       width_mean = 35, constant_baseline = bl,
                       constant_halflife = hl)
  plot(s$data$time, s$data$concentration, type = "l",
       main = sprintf("baseline=%.2f, halflife=%.0f", bl, hl),
       xlab = "Time (min)", ylab = "Conc.")
}
par(op)
```
- **Prose:** vague vs informative — show a wide vs narrow baseline prior and note the trade-off; tie to the identifiability point (reasonably informative priors stabilize the fit). Note that `population_spec()` mirrors these at the population level.

**Step 2: Verify non-ASCII compliance (CRAN).**

Run: `Rscript -e 'tools::showNonASCIIfile("vignettes/03-priors-and-specification.Rmd")'`
Expected: no output (all ASCII; Greek letters via `\uXXXX`).

**Step 3: Verify it builds (fast, live).**

Run: `cd R-package && Rscript -e 'devtools::build_vignettes()'` then `rm -rf R-package/doc R-package/Meta`
Expected: builds in seconds, no errors.

**Step 4: Commit.**
```bash
git add vignettes/03-priors-and-specification.Rmd
git commit -m "docs(vignette): add priors-and-specification vignette"
```

---

## Phase 4: Vignette D — Assessing convergence (single-subject + generic tools)

Uses the precompute pattern (runs several MCMC chains). Slug: `04-convergence-diagnostics`. Uses only model-agnostic diagnostics + `bp_trace`; forward-references the population suite to E.

### Task 4.1: Write `04-convergence-diagnostics.Rmd.orig`

**Files:**
- Create: `vignettes/04-convergence-diagnostics.Rmd.orig`

**Step 1: Write the source.**
````markdown
---
title: "Assessing convergence and MCMC diagnostics"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Assessing convergence and MCMC diagnostics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.path = "04-convergence-diagnostics-",
                      fig.width = 7, fig.height = 5, dpi = 96)
library(bayespulse)
```

<!-- PROSE: convergence = chains forgot their start and explore the posterior.
     We check it with trace plots, effective sample size (ESS), and the
     Gelman-Rubin statistic across multiple chains. Half-life is the hardest
     parameter to mix here -- we show that honestly. -->

## Run several chains

```{r chains}
simdat <- {set.seed(2026); simulate_pulse(num_obs = 60, interval = 10,
            error_var = 0.005, mass_mean = 2.5, width_mean = 40,
            constant_halflife = 40, constant_baseline = 2.5)}
spec <- pulse_spec(
  location_prior_type = "strauss", prior_location_gamma = 0.1,
  prior_location_range = 30, prior_halflife_mean = 38, prior_halflife_var = 1000,
  prior_sd_mass = 0.5, prior_sd_width = 5, prior_mass_mean = 0.5,
  prior_mass_var = 2, prior_mean_pulse_count = 10, prior_width_mean = 35,
  prior_baseline_mean = 2.25, prior_baseline_var = 100, sv_mass_mean = 2,
  sv_width_sd = 10, sv_error_var = 0.005, sv_mass_sd = 0.4, sv_baseline_mean = 2,
  sv_halflife_mean = 42, pv_mean_pulse_width = 10, pv_indiv_pulse_width = 10,
  pv_sd_pulse_width = 0.1, pv_pulse_location = 5, pv_sd_pulse_mass = 0.45)

fits <- lapply(1:4, function(k) {
  set.seed(100 + k)
  fit_pulse(data = simdat$data, spec = spec,
            iters = 40000, thin = 20, burnin = 10000, verbose = FALSE)
})
```

## Trace plots (visual mixing)

```{r trace}
bp_trace(fits[[1]])
```

## Effective sample size

```{r ess}
params <- c("baseline", "halflife", "mass_mean", "width_mean")
ess_tbl <- sapply(params, function(p)
  effective_sample_size(patient_chain(fits[[1]])[[p]]))
round(ess_tbl, 1)
```

## Gelman-Rubin across chains

```{r gr}
gr <- sapply(params, function(p) {
  chains <- lapply(fits, function(f) patient_chain(f)[[p]])
  gelman_rubin(chains)
})
round(gr, 3)
```

<!-- PROSE: interpret -- ESS > ~100-400 good; R-hat near 1.0 good, >1.1 = not
     converged. Point out half-life as the laggard and list remedies (more
     iterations, more thinning, tighter half-life prior). Then 1 sentence:
     "For population models, bayespulse also provides turnkey convergence_report(),
     plot_trace(), and plot_acf() -- see the population vignette." -->
````

**Step 2: Verify the source parses.**

Run: `Rscript -e 'knitr::purl("vignettes/04-convergence-diagnostics.Rmd.orig", output=tempfile(), quiet=TRUE); cat("parsed OK\n")'`
Expected: `parsed OK`.

**Step 3: Commit source.**
```bash
git add vignettes/04-convergence-diagnostics.Rmd.orig
git commit -m "docs(vignette): add convergence-diagnostics source (.Rmd.orig)"
```

### Task 4.2: Precompute vignette D

**Step 1: Knit (runs 4 chains; ~16 min).**

Run: `Rscript vignettes/precompute.R convergence`
Expected: `Knitting 04-convergence-diagnostics.Rmd.orig -> 04-convergence-diagnostics.Rmd`; creates the `.Rmd` + `04-convergence-diagnostics-*.png`. **Sanity-check the printed values:** ESS positive; most R-hat near 1.0; half-life the worst — that is the expected, honest result. If R-hat is wildly off for every parameter, increase `iters` before shipping.

**Step 2: Verify no live chunks remain.**

Run: `grep -c '```{r' vignettes/04-convergence-diagnostics.Rmd`
Expected: `0`.

**Step 3: Commit artifacts.**
```bash
git add vignettes/04-convergence-diagnostics.Rmd vignettes/04-convergence-diagnostics-*.png
git commit -m "docs(vignette): precompute convergence-diagnostics results"
```

---

## Phase 5: Whole-package verification

### Task 5.1: Full check with all four vignettes

**Step 1: Regenerate docs/exports (per CLAUDE.md pre-commit).**

Run: `cd R-package && Rscript -e 'Rcpp::compileAttributes(); devtools::document()'`
Expected: no errors.

**Step 2: Run R CMD check.**

Run: `cd R-package && Rscript -e 'devtools::check(vignettes = TRUE)'`
Expected: **0 ERRORs, 0 WARNINGs.** Vignettes build fast (precomputed A/D do no MCMC; B/C are trivial). Review NOTEs. Acceptable: the symlinked `vignettes/` may raise a NOTE on some platforms; if `R CMD check` rejects the symlink outright, fall back to making `R-package/vignettes` a real directory whose `.Rmd`/`.png` files are copies (and document the divergence from the symlink convention).

**Step 3: Run the test suite.**

Run: `cd R-package && Rscript -e 'devtools::test()'`
Expected: all tests pass.

**Step 4: Final commit (only if Step 1 changed generated files).**
```bash
git add -A R-package
git commit -m "docs(vignette): regenerate package docs after vignette additions" --allow-empty
```

---

## Notes & deferred items

- **Vignettes E–H** (population, joint, simulation-study, chains) follow in a later plan; F (joint) is the flagship.
- **Future vignettes I–O** are catalogued in the roadmap doc and remain out of scope.
- **Deferred package fix (logged 2026-06-07):** extend `plot_trace()`, `plot_acf()`, `convergence_report()`, and `population_diagnostics()` to accept single-subject `pulse_fit` objects (currently population-only). When done, vignette D can adopt the turnkey suite directly. See the roadmap doc's deferred section.
```
