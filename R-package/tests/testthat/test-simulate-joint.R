context("Coupled driver-response simulation (simulate_pulse_joint)")

# Mean distance from each response pulse to its nearest driver pulse.
# Under stronger coupling, response pulses should sit closer to driver pulses.
nearest_driver_dist <- function(sim) {
  r <- sim$response_parameters$location
  d <- sim$driver_parameters$location
  if (length(r) == 0 || length(d) == 0) return(NA_real_)
  mean(vapply(r, function(x) min(abs(x - d)), numeric(1)))
}

test_that("simulate_pulse_joint() returns a well-formed joint_sim object", {
  sim <- simulate_pulse_joint(rho = 0.8, nu = 100, seed = 1)

  expect_s3_class(sim, "joint_sim")
  expect_true(all(c("time", "concentration") %in% names(sim$driver_data)))
  expect_true(all(c("time", "concentration") %in% names(sim$response_data)))
  expect_true(all(c("rho", "nu") %in% names(sim$association)))
  expect_equal(sim$association$rho, 0.8)
  expect_equal(sim$association$nu, 100)

  # Concentrations are strictly positive (multiplicative log-normal error)
  expect_true(all(sim$driver_data$concentration > 0))
  expect_true(all(sim$response_data$concentration > 0))

  # At least one driver pulse so coupling is defined
  expect_gte(sim$association$n_driver_pulses, 1)
})

test_that("simulate_pulse_joint() is reproducible given a seed", {
  s1 <- simulate_pulse_joint(rho = 1, nu = 100, seed = 123)
  s2 <- simulate_pulse_joint(rho = 1, nu = 100, seed = 123)

  expect_equal(s1$response_data$concentration, s2$response_data$concentration)
  expect_equal(s1$response_parameters$location, s2$response_parameters$location)
  expect_equal(s1$driver_parameters$location, s2$driver_parameters$location)
})

test_that("stronger coupling clusters response pulses nearer driver pulses", {
  # Driver is identical across the two coupling strengths for a given seed
  # (same RNG stream through driver simulation); only rho differs. Average
  # over several seeds to avoid stochastic flakiness.
  strong <- numeric(0)
  weak   <- numeric(0)
  for (s in 1:12) {
    strong <- c(strong,
                nearest_driver_dist(simulate_pulse_joint(rho = 2.5, nu = 50,
                                                          response_pulse_count = 8, seed = s)))
    weak   <- c(weak,
                nearest_driver_dist(simulate_pulse_joint(rho = 0.02, nu = 50,
                                                         response_pulse_count = 8, seed = s)))
  }
  expect_lt(mean(strong, na.rm = TRUE), mean(weak, na.rm = TRUE))
})

test_that("fit_pulse_joint() runs on simulated coupled data with positive, mixing rho/nu", {
  # NOTE: this is a BEHAVIORAL test, not a recovery test. The coupling-parameter
  # likelihood currently under-identifies rho/nu (see ISSUES / known limitation);
  # a recovery test (rho/nu within tolerance of truth) should be added once the
  # coupling log-likelihood is corrected.
  sim <- simulate_pulse_joint(rho = 1.5, nu = 100, response_pulse_count = 10,
                              num_obs = 120, seed = 7)
  fit <- fit_pulse_joint(driver_data = sim$driver_data,
                         response_data = sim$response_data,
                         spec = joint_spec(),
                         iters = 1500, thin = 5, burnin = 500, verbose = FALSE)

  expect_s3_class(fit, "joint_fit")
  expect_true(all(fit$association_chain$rho > 0))
  expect_true(all(fit$association_chain$nu > 0))
  # Chains actually move (not frozen)
  expect_gt(diff(range(fit$association_chain$rho)), 0)
  expect_gt(diff(range(fit$association_chain$nu)), 0)
})
