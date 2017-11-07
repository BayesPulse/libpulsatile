# poppulsatile
#

setwd("~/Projects/BayesPulse/poppulsatile/")
library(devtools)
# library(testthat)

devtools::document()
devtools::check()
devtools::install()

# devtools::use_testthat()
testthat::use_catch()
