# poppulsatile
#

setwd("~/Projects/BayesPulse/Software/libpulsatile/R-package/")
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
# library(testthat)

source("./build_package.R", echo=TRUE)

Rcpp::compileAttributes()
devtools::document()
devtools::check()
devtools::install()

devtools::test()

# devtools::use_testthat()
# testthat::use_catch()

# -I"./include/datastructures" -I"./include/mcmc" -I"./include/population" -I"./include/singlesubject" -I"./include/testing" 


library(bayespulse)

set.seed(2018-06-21)
sim <- simulate_pulse()
spec <- pulse_spec()
fit <- fit_pulse(sim, spec = spec, iters = 100, thin = 1, verbose = TRUE)

# Doesn't yet work -- not yet correct class
bp_trace(fit)


