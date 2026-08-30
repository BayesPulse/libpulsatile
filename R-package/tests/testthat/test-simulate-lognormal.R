context("Log-normal simulate path")

test_that("pulse_distribution = 'lognormal' draws log-scale mass/width", {
  # Fast, deterministic (set.seed) consistency check of the log-normal DGP: with
  # many pulses the sample mean/SD of log(mass) and log(width) match the
  # requested log-scale mean/SD. This is the data-generating side of the
  # papers-default fit (pulse width = exp(N(mu_width, sigma_width^2))).
  set.seed(1)
  s <- simulate_pulse(num_obs = 600, interval = 10, ipi_mean = 12,
                      mass_mean = 1.2, mass_sd = 0.5,
                      width_mean = 3.0, width_sd = 0.5,
                      constant_halflife = 45, constant_baseline = 2.6,
                      pulse_distribution = "lognormal")
  expect_gt(nrow(s$parameters), 20)          # enough pulses to estimate SDs
  expect_true(all(s$parameters$width > 0))
  expect_true(all(s$parameters$mass  > 0))

  lw <- log(s$parameters$width)
  lm <- log(s$parameters$mass)
  # Generous, seed-fixed tolerances -- deterministic, so not flaky.
  expect_equal(mean(lw), 3.0, tolerance = 0.25)
  expect_equal(stats::sd(lw), 0.5, tolerance = 0.20)
  expect_equal(mean(lm), 1.2, tolerance = 0.25)
  expect_equal(stats::sd(lm), 0.5, tolerance = 0.20)

  # Under lognormal there is no per-pulse t-scale mixture: kappa columns are 1.
  expect_true(all(s$parameters$width_kappa == 1))
  expect_true(all(s$parameters$mass_kappa  == 1))
})

test_that("simulate_pulse(pulse_distribution='lognormal') defaults give sane widths", {
  # Regression: with natural-scale defaults the lognormal branch drew
  # width = exp(rnorm(35, 5)) ~ 1e15. Default args must now resolve to LOG-scale
  # (width_mean = 3.0, width_sd = 0.7), so widths are physically reasonable and
  # log(width) ~ N(3.0, 0.7).
  set.seed(3)
  s <- simulate_pulse(num_obs = 600, interval = 10, ipi_mean = 12,
                      pulse_distribution = "lognormal")
  expect_gt(nrow(s$parameters), 20)
  expect_true(all(is.finite(s$parameters$width)))
  expect_true(all(is.finite(s$parameters$mass)))
  # Median width near exp(3.0) ~ 20, and nowhere near the 1e15 blow-up.
  expect_lt(stats::median(s$parameters$width), 200)
  expect_gt(stats::median(s$parameters$width), 2)
  lw <- log(s$parameters$width)
  expect_equal(mean(lw), 3.0, tolerance = 0.3)
  expect_equal(stats::sd(lw), 0.7, tolerance = 0.3)
  # Mass defaults resolve to log-scale too (mean log-mass ~ 1.2).
  expect_equal(mean(log(s$parameters$mass)), 1.2, tolerance = 0.3)
})

test_that("simulate_pulse() default is unchanged (truncnorm)", {
  # Backward compatibility: the default path must reproduce the legacy
  # truncated-Student-t draw byte-for-byte for a fixed seed.
  set.seed(99)
  a <- simulate_pulse(num_obs = 100)
  set.seed(99)
  b <- simulate_pulse(num_obs = 100, pulse_distribution = "truncnorm")
  expect_equal(a$parameters$width, b$parameters$width)
  expect_equal(a$parameters$mass,  b$parameters$mass)
})

test_that("simulate_pulse_population() threads pulse_distribution", {
  ps <- simulate_pulse_population(n_subjects = 2, num_obs = 120,
                                  mass_mean = 1.2, mass_sd = 0.5,
                                  width_mean = 3.0, width_sd = 0.5,
                                  pulse_distribution = "lognormal", seed = 7)
  expect_s3_class(ps, "population_sim")
  expect_length(ps$data, 2)
  expect_true(all(is.finite(ps$data[[1]]$concentration)))
})
