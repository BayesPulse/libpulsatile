# BayesPulse
#

# setwd("~/Projects/BayesPulse/Software/libpulsatile/R-package/")
setwd("~/Projects/BayesPulse/Software/libpulsatile/R-package/")
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(tidyverse)
# library(testthat)

source("./build_package.R", echo=TRUE)

Rcpp::compileAttributes()
devtools::document()
devtools::check()
devtools::install()

#devtools::test()



library(bayespulse)

# load testing series: test_data
source("./data-raw/test_data.R")

set.seed(2018-09-21)
sim <- simulate_pulse()
plot(sim)
spec <- pulse_spec(location_prior_type = "strauss", #order-statistic") #, 
                   prior_location_gamma = 0,
                   prior_location_range = 40)
fit <- fit_pulse(sim, spec = spec, iters = 250000, thin = 10, burnin = 10000,
                 verbose = TRUE)
bp_trace(fit) + ggtitle("BAYESPULSE")
bp_posteriors(fit, "histogram") + ggtitle("BAYESPULSE")
bp_location_posterior(fit) + ggtitle("BAYESPULSE")

fit_predicted <- predict(fit)
bp_predicted(fit_predicted)







pulse_chain(fit) %>%
  group_by(iteration) %>%
  mutate(timediff = location - lag(location)) #%>%
#   filter(pulse_num == 3) 
#   summarise(any_lt_range = any(timediff < 100)) %>%
#   with(., table(any_lt_range))

# Plot diagnostics
bp_trace(fit) + ggtitle("BAYESPULSE")
bp_posteriors(fit) + ggtitle("BAYESPULSE")

plot(test_data) +
  geom_vline(data = test_data$parameters, aes(xintercept = location)) + ggtitle("BAYESPULSE")

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



conc <- c(0.9864, 1.2863, 1.5225, 1.6848, 1.5266, 1.7643, 1.6280, 1.4952,
          1.4419, 1.3216, 1.3234, 1.3335, 1.2407, 1.2850, 1.0917, 1.1742,
          1.2177, 1.1262, 1.1937, 1.2434, 1.0882, 1.2241, 1.1462, 1.0765,
          1.1123, 1.1509, 1.0083, 1.0543, 1.1745, 1.6471, 1.7496, 1.7042,
          1.7095, 1.5984, 1.5967, 1.9896, 1.9641, 1.9172, 2.0328, 1.8320,
          1.7972, 1.5481, 1.6172, 1.4398, 1.5102, 1.4086, 1.2399, 1.1592,
          1.4520, 1.9262, 1.7840, 1.8651, 1.9138, 1.6058, 1.5749, 1.4671,
          1.4577, 1.4440, 1.2318, 1.2989, 1.0915, 1.2266, 0.9714, 1.3056,
          1.1570, 1.2522, 1.1307, 1.1650, 1.6019, 1.7156, 1.8527, 1.6043,
          1.6282, 1.5820, 1.5795, 1.6706, 1.7033, 1.6826, 1.5780, 1.5759,
          1.5092, 1.4379, 1.5565, 1.8580, 2.0752, 1.9787, 1.8395, 1.7592,
          1.7275, 1.5444, 1.6046, 1.5941, 1.7257, 1.6003, 1.4624, 1.4928,
          1.3952, 1.7805, 1.6982, 1.6382, 1.6487, 1.5158, 1.4884, 1.3650,
          1.3769, 1.2293, 1.2044, 1.1552, 1.2265, 1.0377, 1.1934, 0.8843,
          1.2278, 1.0475, 1.1364, 0.9983, 1.0791, 1.0605, 1.0922, 1.2886,
          1.5316, 1.4970, 1.5719, 1.3826, 1.2423, 1.3035, 1.2633, 1.2709,
          1.2677, 1.3222, 1.0761, 1.4432, 1.5511, 1.3972, 1.4633, 1.3793,
          1.3484, 1.3988, 1.1616, 1.1532, 1.1025, 1.0011, 1.1876, 1.1246)

mean_conc <- c(1.7381, 1.6574, 1.5827, 1.6280, 2.2622, 2.1556, 2.0505, 2.0496,
               2.9221, 2.7956, 2.6677, 2.5433, 2.4227, 2.3064, 2.1946, 2.0876,
               1.9859, 1.8896, 1.7989, 1.7141, 1.6352, 1.5622, 1.4950, 1.4337,
               1.3779, 1.3274, 1.2821, 1.2415, 1.2054, 1.1733, 1.2413, 2.3021,
               2.1915, 2.0847, 1.9831, 1.8869, 1.7965, 1.7118, 1.6330, 1.5602,
               1.4932, 1.4320, 1.3764, 1.3261, 1.2809, 1.2991, 1.9788, 1.8838,
               1.7935, 1.7090, 1.6305, 1.5578, 1.4910, 1.4300, 1.3746, 1.3245,
               1.2794, 1.2391, 1.2032, 1.1714, 1.2867, 2.6507, 2.5271, 2.4070,
               2.2913, 2.1801, 2.0738, 1.9728, 1.8772, 1.7874, 1.7033, 1.6252,
               1.5529, 1.4866, 1.4259, 1.3709, 1.3211, 1.2764, 1.2365, 1.2009,
               1.3060, 1.9359, 1.8436, 1.7558, 1.6739, 1.5980, 1.5279, 1.4637,
               1.4051, 1.3520, 1.3042, 1.2612, 1.2229, 1.1888, 1.1587, 1.1321,
               1.1087, 1.2776, 1.8287, 1.7421, 1.6611, 1.5861, 1.5170, 1.4537,
               1.3961, 1.3439, 1.2968, 1.2547, 1.2171, 1.1837, 1.1541, 1.1281,
               1.1052, 1.0851, 1.0676, 1.1337, 1.7971, 1.7160, 1.6370, 1.5638,
               1.4965, 1.4350, 1.3791, 1.3286, 1.2831, 1.2424, 1.2061, 1.1740,
               1.1456, 1.1206, 1.1877, 1.9580, 1.8642, 1.7751, 1.6919, 1.6146,
               1.5432, 1.4777, 1.4178, 1.3635, 1.3145, 1.2705, 1.2312, 1.1989)

internal_conc <- c(5.237937,  5.156800,  5.773434,  6.364134,  9.872640,
                   8.782145, 8.633523, 7.166886,  6.063834,  6.393283,
                   5.547755,  5.609554, 17.809125, 14.864233, 14.488007,
                   11.012978, 10.074941,  8.579306,  7.925953, 7.505312,
                   7.408293, 6.296038,  6.138372,  6.114483,  4.589039,
                   4.787283, 4.125096,  4.560712, 3.838159,  3.611227,
                   10.149188,  9.114352,  8.475730, 8.166254,  7.115865,
                   6.364551,  6.076659,  4.936654,  4.815760,  5.000790,
                   4.405922,  4.621510, 3.612532,  3.861247,  3.480776,
                   4.382934,  6.908732, 5.677212,  5.697563, 5.049410,
                   5.033481,  5.036818,  4.766425,  4.108461, 4.493277,
                   3.970064, 3.479031,  3.105488,  3.182802,  2.832481,
                   2.904569, 13.393956, 12.160363, 12.604881,  9.282701,
                   8.971699,  7.872894, 7.341750,  6.129950,  6.754877,
                   6.236324,  5.221481,  4.712401,  4.251457, 4.107650,
                   3.642404,  3.409108, 3.673586,  3.573048,  3.006202,
                   3.574850, 2.984541,  2.942794,  6.142817, 6.239127,
                   5.607103,  5.781841,  5.365100, 4.483500,  4.318313,
                   3.729648, 4.186974,  3.661007,  3.587046,  3.420915,
                   3.554909,  2.951518,  3.197826, 2.768718,  4.050486,
                   6.375893,  5.984611, 4.975292,  4.758750,  5.006260,
                   4.137620,  3.676703,  4.133922,  3.445431, 3.813060,
                   3.304948,  3.308399, 3.362024,  3.239059,  5.554801,
                   5.157033, 5.049831,  5.366527,  4.410822, 4.798887,
                   4.386467,  4.013058,  3.749040, 3.580090,  3.434933,
                   3.086843, 3.471470,  3.334008,  2.913723,  2.980340,
                   2.862694,  5.891940,  6.136203, 6.251388,  6.170612,
                   5.191264,  4.810252, 4.989489,  3.951583,  3.793401,
                   4.320815, 10.128100,  9.961510,  8.155648)


data.frame(conc, mean_conc, internal_conc) %>%
  mutate(time = 1:n()) %>%
#   mutate(mean_conc = log(mean_conc)) %>%
  mutate_at(vars(conc, mean_conc), funs(exp)) %>%
  gather(key = type, value = conc, conc, mean_conc, internal_conc) %>%
#   filter(type == "conc") %>%
  ggplot(aes(x = time, y = conc, color = type)) + geom_path() + geom_point()
plot(test_data)
