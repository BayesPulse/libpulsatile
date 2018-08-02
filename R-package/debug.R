# debug libpulsatile 
#  -- script to call with R -d lldb -e 
#
setwd("~/Projects/BayesPulse/Software/libpulsatile/R-package/")
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)

source("./build_package.R", echo=TRUE)

Rcpp::compileAttributes()
devtools::document()
# devtools::check()
devtools::install()
devtools::test()

library(bayespulse)

set.seed(2018-06-27)
sim <- simulate_pulse()
spec <- pulse_spec()
fit <- fit_pulse(sim, spec = spec, iters = 250000, thin = 50, burnin = 100000, verbose = TRUE)

# Plot diagnostics
plot(sim)
bp_trace(fit)
bp_posteriors(fit)
bp_location_posterior(fit)

# True values
sim$call$mass_mean
sim$call$mass_sd
sim$call$width_mean
sim$call$width_sd
sim$call$constant_baseline
sim$call$constant_halflife
sim$call$mass_sd
sim$call$error_var


# #
# iter <- 1:10000
# output <- function(iteration, burnin = 2000, thin = 50) {
# 
#   outdf <- data.frame(index = integer(), iteration = integer())
# 
#   for (i in iteration)  {
# 
#     if (i > (burnin) & ((i) %% thin == 0))  {
#       index <- ((i - burnin) / thin) - 1
#       outdf[index+1, "index"] <- index
#       outdf[index+1, "iteration"] <- i
#     }
# 
#   }
# 
#   return(outdf)
# 
# }
# 
# output(iter) #%>% length

devtools::unload("bayespulse")
