context("Single-subject fit")

test_that("student_t_pulses = FALSE holds every pulse t-scale (kappa) at 1", {
  skip_on_cran()

  set.seed(2025)
  n_obs <- 24
  dat <- data.frame(
    time = seq(0, 1440, length.out = n_obs),
    concentration = abs(rnorm(n_obs, mean = 3, sd = 1))
  )

  spec <- pulse_spec(student_t_pulses = FALSE, prior_mean_pulse_count = 4)

  fit <- fit_pulse(data = dat, spec = spec,
                   iters = 1000, thin = 10, burnin = 500, verbose = FALSE)

  # Under Gaussian random effects the per-pulse t-scale is fixed at 1 for every
  # pulse (never drawn; initialized to 1 in the birth-death). The single-subject
  # pulse chain reports these as eta_mass / eta_width.
  expect_true("eta_mass" %in% names(fit$pulse_chain))
  expect_true("eta_width" %in% names(fit$pulse_chain))
  expect_true(all(fit$pulse_chain$eta_mass == 1))
  expect_true(all(fit$pulse_chain$eta_width == 1))

  # width_sd is the PR's primary motivation: its draw was commented out and the
  # parameter sat frozen. It is sampled independently of the t-scale toggle, so
  # assert here that both pulse-to-pulse SDs actually move -- a frozen chain
  # (re-introduced `// draw_sd_widths...`) would be constant, giving sd == 0.
  expect_true(stats::sd(fit$patient_chain$width_sd) > 0)
  expect_true(stats::sd(fit$patient_chain$mass_sd) > 0)
})


test_that("default (lognormal) single-subject fit runs end-to-end on log-normal data", {
  skip_on_cran()

  # End-to-end smoke test of the papers-default single-subject parameterization
  # {lognormal pulse RE, uniform SD prior, gaussian kappa}. Simulate log-normal
  # data with a KNOWN log-scale width_sd, then fit with the DEFAULT pulse_spec()
  # so the data-generating and model scales agree.
  set.seed(4321)
  true_width_sd <- 0.7   # LOG-scale pulse-to-pulse width SD used by the DGP
  sim <- simulate_pulse(pulse_distribution = "lognormal",
                        num_obs   = 48,
                        interval  = 10,
                        ipi_mean  = 12,
                        mass_mean = 1.2,
                        mass_sd   = 0.5,
                        width_mean = 3.0,
                        width_sd  = true_width_sd)

  spec <- pulse_spec()  # default: lognormal / uniform / gaussian
  fit <- fit_pulse(data = sim$data, spec = spec,
                   iters = 2500, thin = 10, burnin = 1250, verbose = FALSE)

  # (a) All chains finite -- no NaN leaking from calc_mean_contribution etc.
  expect_true(all(is.finite(unlist(fit$patient_chain))))
  expect_true(all(is.finite(fit$pulse_chain$mass)))
  expect_true(all(is.finite(fit$pulse_chain$width)))

  # (b) The pulse-to-pulse SD chains mix (positive variance, not frozen).
  expect_true(stats::sd(fit$patient_chain$width_sd) > 0)
  expect_true(stats::sd(fit$patient_chain$mass_sd) > 0)

  # (c) width_sd posterior mean lands in a sane neighborhood of the truth. Width
  # is weakly identified, so the tolerance is deliberately wide (not tight/flaky)
  # -- this only guards against gross mis-scaling, not precise recovery.
  width_sd_hat <- mean(fit$patient_chain$width_sd)
  expect_gt(width_sd_hat, 0.1)
  expect_lt(width_sd_hat, 3.0)
})
