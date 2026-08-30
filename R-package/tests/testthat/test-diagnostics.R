context("Convergence Diagnostics Tests")

test_that("effective_sample_size() works correctly", {
  # Well-mixed chain (low autocorrelation)
  set.seed(123)
  x <- rnorm(1000)
  ess <- effective_sample_size(x)
  
  expect_true(is.numeric(ess))
  expect_true(ess > 0)
  expect_true(ess <= length(x))
  # Well-mixed chain should have high ESS
  expect_true(ess > 500)
  
  # Poorly mixed chain (high autocorrelation)
  y <- cumsum(rnorm(1000))  # Random walk
  ess_poor <- effective_sample_size(y)
  
  # Poor mixing should have much lower ESS
  expect_true(ess_poor < ess)
  expect_true(ess_poor < 100)
})


test_that("effective_sample_size() handles edge cases", {
  # Very short chain
  expect_warning(
    ess <- effective_sample_size(rnorm(5)),
    "too small"
  )
  expect_true(is.na(ess))
  
  # Chain with NA values
  x <- c(rnorm(100), NA, NA, rnorm(100))
  ess <- effective_sample_size(x)
  expect_true(is.numeric(ess))
  expect_true(!is.na(ess))
})


test_that("gelman_rubin() calculates R-hat correctly", {
  set.seed(456)
  
  # Converged chains (same distribution)
  chain1 <- rnorm(1000, mean = 5, sd = 2)
  chain2 <- rnorm(1000, mean = 5, sd = 2)
  chain3 <- rnorm(1000, mean = 5, sd = 2)
  
  rhat <- gelman_rubin(list(chain1, chain2, chain3))
  
  expect_true(is.numeric(rhat))
  expect_true(rhat > 0)
  # Converged chains should have R-hat close to 1
  expect_true(rhat < 1.1)
  
  # Non-converged chains (different means)
  chain1_nc <- rnorm(1000, mean = 0, sd = 1)
  chain2_nc <- rnorm(1000, mean = 5, sd = 1)
  chain3_nc <- rnorm(1000, mean = 10, sd = 1)
  
  rhat_nc <- gelman_rubin(list(chain1_nc, chain2_nc, chain3_nc))
  
  # Non-converged chains should have R-hat >> 1
  expect_true(rhat_nc > 2)
})


test_that("gelman_rubin() validates inputs", {
  # Single chain
  expect_error(
    gelman_rubin(list(rnorm(100))),
    "at least 2 chains"
  )
  
  # Not a list
  expect_error(
    gelman_rubin(rnorm(100)),
    "must be a list"
  )
  
  # Different length chains
  expect_warning(
    rhat <- gelman_rubin(list(rnorm(100), rnorm(200))),
    "different lengths"
  )
  expect_true(is.numeric(rhat))
})


test_that("population_diagnostics() works with population_fit", {
  skip_on_cran()
  skip_on_ci()
  
  set.seed(789)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 20))
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  diag <- population_diagnostics(fit)
  
  expect_true(is.data.frame(diag))
  expect_true(nrow(diag) > 0)
  expect_true(all(c("parameter", "mean", "sd", "ess", "ess_per_sample",
                    "lag1_acf", "converged") %in% names(diag)))
  
  # Check ESS is reasonable
  expect_true(all(diag$ess > 0, na.rm = TRUE))
  expect_true(all(diag$ess <= nrow(fit$population_chain), na.rm = TRUE))
  
  # Check lag-1 ACF is between -1 and 1
  expect_true(all(diag$lag1_acf >= -1 & diag$lag1_acf <= 1, na.rm = TRUE))
})


test_that("population_diagnostics() can select parameters", {
  skip_on_cran()
  skip_on_ci()
  
  set.seed(101)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 20))
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  # Select specific parameters
  diag <- population_diagnostics(fit, parameters = c("mass_mean", "width_mean"))
  
  expect_equal(nrow(diag), 2)
  expect_true(all(diag$parameter %in% c("mass_mean", "width_mean")))
})


test_that("convergence_report() produces output", {
  skip_on_cran()
  skip_on_ci()
  
  set.seed(202)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 20))
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  expect_output(convergence_report(fit), "Convergence Diagnostics")
  expect_output(convergence_report(fit), "Population-Level Parameters")
  expect_output(convergence_report(fit), "ESS")
  
  # Should return diagnostics data frame invisibly
  diag <- convergence_report(fit)
  expect_true(is.data.frame(diag))
})


test_that("plot_trace() creates ggplot", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  
  set.seed(303)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 20))
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  p <- plot_trace(fit)
  
  expect_s3_class(p, "ggplot")
  
  # Test with specific parameters
  p2 <- plot_trace(fit, parameters = c("mass_mean", "baseline_mean"))
  expect_s3_class(p2, "ggplot")
})


test_that("plot_acf() creates ggplot", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("ggplot2")
  
  set.seed(404)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 20))
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  p <- plot_acf(fit)
  
  expect_s3_class(p, "ggplot")
  
  # Test with specific parameters
  p2 <- plot_acf(fit, parameters = c("mass_mean"), max_lag = 20)
  expect_s3_class(p2, "ggplot")
})


test_that("summary.population_fit includes diagnostics", {
  skip_on_cran()
  skip_on_ci()
  
  set.seed(505)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 20))
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  # With diagnostics
  expect_output(summary(fit), "Convergence Diagnostics")
  expect_output(summary(fit), "ESS")
  
  # Without diagnostics: still prints the model summary header, but the
  # convergence diagnostics (ESS) section must be suppressed.
  expect_output(summary(fit, diagnostics = FALSE), "Population Model Summary")
  no_diag_output <- capture.output(summary(fit, diagnostics = FALSE))
  expect_false(any(grepl("Convergence Diagnostics", no_diag_output)))
  expect_false(any(grepl("ESS", no_diag_output)))
})


test_that("identifiability_check() flags posteriors that fill their prior", {
  set.seed(42)

  # An unidentified parameter: nearly uniform over its whole prior [0, 30].
  # A well-identified one: tightly concentrated well inside the prior.
  chain <- data.frame(
    iteration      = seq_len(2000),
    width_sd       = runif(2000, 0, 30),
    mass_mean      = rnorm(2000, mean = 3.5, sd = 0.1)
  )

  prior_bounds <- list(
    width_sd  = c(0, 30),
    mass_mean = c(0, 30)
  )

  res <- identifiability_check(chain, prior_bounds = prior_bounds)

  expect_s3_class(res, "data.frame")
  # iteration column is dropped from the analysis
  expect_false("iteration" %in% res$parameter)
  expect_true(all(c("prior_coverage", "weak_identifiability") %in% names(res)))

  width_row <- res[res$parameter == "width_sd", ]
  mass_row  <- res[res$parameter == "mass_mean", ]

  # The space-filling parameter should cover most of its prior and be flagged.
  expect_true(width_row$prior_coverage > 0.8)
  expect_true(isTRUE(width_row$weak_identifiability))

  # The concentrated parameter covers a tiny fraction and is not flagged.
  expect_true(mass_row$prior_coverage < 0.1)
  expect_false(isTRUE(mass_row$weak_identifiability))
})


test_that("identifiability_check() returns NA verdict without prior bounds", {
  chain <- data.frame(a = rnorm(500), b = rnorm(500))
  res <- identifiability_check(chain)

  expect_true(all(is.na(res$prior_coverage)))
  expect_true(all(is.na(res$weak_identifiability)))
  # Posterior summaries are still populated.
  expect_true(all(is.finite(res$post_sd)))
})


test_that("identifiability_check() accepts a single-subject pulse_fit", {
  set.seed(7)
  sim <- simulate_pulse(num_obs = 30, interval = 10,
                        pulse_distribution = "lognormal")
  spec <- pulse_spec(location_prior_type = "strauss",
                     prior_location_gamma = 0.1, prior_location_range = 30,
                     prior_mean_pulse_count = 5)
  fit <- fit_pulse(data = sim$data, spec = spec,
                   iters = 2000, thin = 10, burnin = 500, verbose = FALSE)

  res <- identifiability_check(fit)

  expect_s3_class(res, "data.frame")
  expect_false("iteration" %in% res$parameter)
  expect_true(all(c("mass_sd", "width_sd") %in% res$parameter))

  # The uniform SD priors give bounded support, so those rows get real
  # verdicts derived automatically from spec$priors (Uniform(0, 10) here).
  sd_rows <- res[res$parameter %in% c("mass_sd", "width_sd"), ]
  expect_true(all(sd_rows$prior_range == 10))
  expect_true(all(is.finite(sd_rows$prior_coverage)))
  expect_true(all(!is.na(sd_rows$weak_identifiability)))

  # Parameters without a bounded prior stay NA (unchecked, not "fine").
  bl_row <- res[res$parameter == "baseline", ]
  expect_true(is.na(bl_row$prior_coverage))

  # Explicit prior_bounds still override/supplement the auto-derived ones.
  res2 <- identifiability_check(fit,
                                prior_bounds = list(baseline = c(0, 20)))
  bl2 <- res2[res2$parameter == "baseline", ]
  expect_equal(bl2$prior_range, 20)
  sd2 <- res2[res2$parameter == "mass_sd", ]
  expect_equal(sd2$prior_range, 10)
})


test_that("identifiability_check() derives no SD bounds under half_cauchy", {
  set.seed(7)
  sim <- simulate_pulse(num_obs = 30, interval = 10,
                        pulse_distribution = "lognormal")
  spec <- pulse_spec(location_prior_type = "strauss",
                     prior_location_gamma = 0.1, prior_location_range = 30,
                     prior_mean_pulse_count = 5,
                     sd_prior = "half_cauchy")
  fit <- fit_pulse(data = sim$data, spec = spec,
                   iters = 2000, thin = 10, burnin = 500, verbose = FALSE)

  res <- identifiability_check(fit)

  # A half-Cauchy prior has unbounded support: sd_max is not a prior bound
  # there, so no automatic verdict may be derived from it.
  sd_rows <- res[res$parameter %in% c("mass_sd", "width_sd"), ]
  expect_true(all(is.na(sd_rows$prior_range)))
  expect_true(all(is.na(sd_rows$weak_identifiability)))
})
