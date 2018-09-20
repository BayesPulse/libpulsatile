
# library(bayespulse)
set.seed(9999)

# Generate test data for use in c++ unit tests
test_data <-
  bayespulse::simulate_pulse(num_obs = 144,
                             interval = 10,
                             error_var = 0.005,
                             ipi_mean = 12,
                             ipi_var = 40,
                             ipi_min = 4,
                             mass_mean = 3.5,
                             mass_sd = 1.6,
                             width_mean = 42,
                             width_sd = 35,
                             constant_halflife = 45,
                             constant_baseline = 2.6,
                             halflife_mean = NULL,
                             halflife_var = NULL,
                             baseline_mean = NULL,
                             baseline_var = NULL)
plot(test_data)

save(test_data, file = "test_data.RData")

