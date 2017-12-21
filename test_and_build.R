# poppulsatile
#

setwd("~/Projects/BayesPulse/poppulsatile/")
library(devtools)
library(Rcpp)
library(RcppArmadillo)
# library(testthat)

devtools::test()

devtools::document()
devtools::check()
devtools::install()

# devtools::use_testthat()
# testthat::use_catch()
