# poppulsatile
#

setwd("~/Projects/BayesPulse/libpulsatile/R-package/")
library(devtools)
library(Rcpp)
library(RcppArmadillo)
# library(testthat)


Rcpp::compileAttributes()
devtools::document()
devtools::check()
devtools::install()

devtools::test()

devtools::use_testthat()
testthat::use_catch()

# -I"./include/datastructures" -I"./include/mcmc" -I"./include/population" -I"./include/singlesubject" -I"./include/testing" 
