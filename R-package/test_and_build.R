# poppulsatile
#

# setwd("~/Projects/BayesPulse/Software/libpulsatile/R-package/")
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

set.seed(2018-08-28)
sim <- simulate_pulse()
spec <- pulse_spec(location_prior_type = "strauss", 
                   prior_location_gamma = 0,
                   prior_location_range = 30)
fit <- fit_pulse(sim, spec = spec, iters = 250000, thin = 50, burnin = 100000, verbose = TRUE)

# Plot diagnostics
bp_trace(fit)
bp_posteriors(fit)

bp_location_posterior(fit)
plot(sim) +
  geom_vline(data = sim$parameters, aes(xintercept = location))

fit$pulse_chain %>%
  ggplot(aes(x = iteration, y = eta_width)) + geom_path(alpha = 0.2)



#---------------------------------------
# detailed debugging
#---------------------------------------
set.seed(999)
sim    <- bayespulse::simulate_pulse()
myspec <- bayespulse::pulse_spec() # don't like how this stores just indicator of strauss_location_prior instead of string value of prior type
fit    <- bayespulse::fit_pulse(sim, spec = myspec, iters = 51, burnin = 0, thin = 1, verbose = TRUE)


library(tidyverse)
# library(magrittr)
fit$patient_chain %>% .[1:10, ]
fit$patient_chain %>% as.data.frame %>%
  ggplot(aes(x = iteration, y = baseline)) + geom_path()
fit$patient_chain %>% as_data_frame %>%
  ggplot(aes(x = V1, y = V5)) + geom_path()
fit$patient_chain %>% as_data_frame %>%
  ggplot(aes(x  = V4)) + geom_histogram()
fit$patient_chain %>% as_data_frame %>%
  ggplot(aes(x  = V7)) + geom_histogram()

# NOTE: birth death causing malloc error



