context("pulse_spec")

test_that("pulse_spec() defaults to Gaussian (non-Student-t) pulse random effects", {
  spec <- pulse_spec()
  expect_true("student_t_pulses" %in% names(spec$priors))
  expect_false(isTRUE(spec$priors$student_t_pulses))
})

test_that("pulse_spec() can opt in to Student-t pulse random effects", {
  spec <- pulse_spec(student_t_pulses = TRUE)
  expect_true(isTRUE(spec$priors$student_t_pulses))
})

test_that("pulse_spec() carries a Gaussian (kappa=1) request through the priors", {
  spec <- pulse_spec(student_t_pulses = FALSE)
  expect_false(isTRUE(spec$priors$student_t_pulses))
})

test_that("pulse_spec() validates student_t_pulses", {
  expect_error(pulse_spec(student_t_pulses = "yes"), "single logical")
  expect_error(pulse_spec(student_t_pulses = c(TRUE, FALSE)), "single logical")
  expect_error(pulse_spec(student_t_pulses = NA), "single logical")
})

test_that("pulse_spec() defaults to the lognormal / uniform parameterization", {
  spec <- pulse_spec()
  # Both axis flags are carried through the priors list for C++ to read.
  expect_true(isTRUE(spec$priors$lognormal_pulses))
  expect_true(isTRUE(spec$priors$uniform_sd_prior))
})

test_that("pulse_spec() threads truncnorm / half_cauchy flags", {
  spec <- pulse_spec(pulse_distribution = "truncnorm", sd_prior = "half_cauchy")
  expect_false(isTRUE(spec$priors$lognormal_pulses))
  expect_false(isTRUE(spec$priors$uniform_sd_prior))
})

test_that("pulse_spec() validates the parameterization axes", {
  expect_error(pulse_spec(pulse_distribution = "banana"), "should be one of")
  expect_error(pulse_spec(sd_prior = "flat"), "should be one of")
})

test_that("pulse_spec() lognormal defaults are paper log-scale values", {
  spec <- pulse_spec()  # lognormal, uniform, gaussian
  # Log-scale prior means/vars from the thesis.
  expect_equal(spec$priors$mass_mean, 1.2)
  expect_equal(spec$priors$mass_variance, 25)
  expect_equal(spec$priors$width_mean, 3.5)
  expect_equal(spec$priors$width_variance, 25)
  # Uniform(0, 10) SD bounds threaded as the *_sd_max keys the C++ uniform
  # branch reads.
  expect_equal(spec$priors$mass_sd_max, 10)
  expect_equal(spec$priors$width_sd_max, 10)
  # Log-scale starting values, all inside the Uniform(0, 10) support.
  expect_equal(spec$starting_values$mass_mean, 1.0)
  expect_equal(spec$starting_values$width_mean, 3.0)
  expect_equal(spec$starting_values$mass_sd, 0.5)
  expect_equal(spec$starting_values$width_sd, 0.7)
  # Fixed-effect mean / SD proposals act on LOG-scale parameters, so the
  # natural-scale pv (3700 / 4000) shrink to O(0.01-1).
  expect_true(spec$proposal_variances$width_mean  <= 1)
  expect_true(spec$proposal_variances$width_sd    <= 1)
  # The individual-pulse proposal, however, is a random walk on the NATURAL-scale
  # pulse value (log-normal draw), so it must be scaled to the natural width
  # magnitude exp(width_mean) ~ 20 -- NOT shrunk to log-scale sizes. A sim study
  # (benchmarks/sim_study_lognormal_report.md) confirmed shrinking it stalls the
  # pulse-to-pulse width-SD mixing; the paper default is 9 here.
  expect_equal(spec$proposal_variances$pulse_width, 9)
  expect_equal(spec$proposal_variances$pulse_mass,  0.5)
})

test_that("pulse_spec() rejects a starting SD outside the Uniform prior support", {
  # Under the default uniform SD prior a starting SD must be < its max bound.
  expect_error(pulse_spec(sv_width_sd = 15, prior_sd_width = 10),
               "outside the Uniform")
  # Half-Cauchy has no hard upper bound, so the same values are accepted.
  expect_s3_class(
    pulse_spec(sd_prior = "half_cauchy", sv_width_sd = 15, prior_sd_width = 10),
    "pulse_spec")

  # Parameterization-RESOLVED bound: the case above fixes prior_sd_width
  # explicitly, so it never exercises the lognormal-resolved bound. Under
  # pulse_distribution = "lognormal" (no explicit prior_sd_width) the Uniform
  # bound resolves to 10, so a starting SD of 12 is out of support and must error
  # (under "truncnorm" the bound would resolve to 5).
  expect_error(pulse_spec(pulse_distribution = "lognormal", sv_width_sd = 12),
               "outside the Uniform")
  # In-support log-scale starting SD (below the resolved bound of 10) passes.
  expect_s3_class(pulse_spec(pulse_distribution = "lognormal", sv_width_sd = 5),
                  "pulse_spec")
})

test_that("user-supplied values override the parameterization defaults", {
  spec <- pulse_spec(prior_width_mean = 99, sv_mass_mean = 7)
  expect_equal(spec$priors$width_mean, 99)
  expect_equal(spec$starting_values$mass_mean, 7)
})

test_that("truncnorm + half_cauchy + student-t reproduces the pre-change spec", {
  # REGRESSION CONTRACT: the legacy natural-scale defaults must be recovered
  # exactly when the three axes are set to their pre-change (opt-in) values.
  spec <- pulse_spec(pulse_distribution = "truncnorm",
                     sd_prior = "half_cauchy",
                     student_t_pulses = TRUE)

  # Priors (natural scale)
  expect_equal(spec$priors$mass_mean, 3.5)
  expect_equal(spec$priors$mass_variance, 100)
  expect_equal(spec$priors$width_mean, 42)
  expect_equal(spec$priors$width_variance, 1000)
  expect_equal(spec$priors$mass_sd_param, 5)
  expect_equal(spec$priors$width_sd_param, 5)
  expect_equal(spec$priors$baseline_mean, 2.6)
  expect_equal(spec$priors$baseline_variance, 100)
  expect_equal(spec$priors$halflife_mean, 45)
  expect_equal(spec$priors$halflife_variance, 100)
  expect_equal(spec$priors$error_alpha, 0.0001)
  expect_equal(spec$priors$error_beta, 0.0001)
  expect_equal(spec$priors$pulse_count, 12)
  expect_true(isTRUE(spec$priors$student_t_pulses))
  expect_false(isTRUE(spec$priors$lognormal_pulses))
  expect_false(isTRUE(spec$priors$uniform_sd_prior))

  # Starting values (natural scale)
  expect_equal(spec$starting_values$mass_mean, 3.5)
  expect_equal(spec$starting_values$width_mean, 42)
  expect_equal(spec$starting_values$mass_sd, 1.6)
  expect_equal(spec$starting_values$width_sd, 35)
  expect_equal(spec$starting_values$baseline, 2.6)
  expect_equal(spec$starting_values$halflife, 45)
  expect_equal(spec$starting_values$errorsq, 0.005)

  # Proposal variances (natural scale)
  expect_equal(spec$proposal_variances$mass_mean, 6)
  expect_equal(spec$proposal_variances$width_mean, 3700)
  expect_equal(spec$proposal_variances$mass_sd, 4.5)
  expect_equal(spec$proposal_variances$width_sd, 4000)
  expect_equal(spec$proposal_variances$pulse_mass, 1)
  expect_equal(spec$proposal_variances$pulse_width, 15000)
  expect_equal(spec$proposal_variances$baseline, 0.02)
  expect_equal(spec$proposal_variances$halflife, 1.5)
  expect_equal(spec$proposal_variances$location, 65)
  expect_equal(spec$proposal_variances$sdscale_pulse_mass, 4)
  expect_equal(spec$proposal_variances$sdscale_pulse_width, 4)
})
