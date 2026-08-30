context("Chain inventory (chains())")

small_pulse_fit <- function() {
  set.seed(42)
  sim <- simulate_pulse(num_obs = 20, interval = 10)
  fit_pulse(data = sim$data, spec = pulse_spec(), iters = 500, thin = 5,
            burnin = 100, verbose = FALSE)
}

# Row lookup helper: the inventory is keyed by the `chain` column.
row_for <- function(inventory, chain) inventory[inventory$chain == chain, ]


# Single-subject ----------------------------------------------------------

test_that("chains() inventories a single-subject fit", {

  skip_on_cran()

  fit <- small_pulse_fit()
  inventory <- chains(fit)

  expect_s3_class(inventory, "data.frame")
  expect_named(inventory, c("chain", "access", "rows", "columns"))
  expect_setequal(inventory$chain, c("patient", "pulse"))
})


test_that("chains() prints an inventory and returns it invisibly", {

  skip_on_cran()

  fit <- small_pulse_fit()

  expect_output(chains(fit), "patient")
  expect_output(chains(fit), "patient_chain\\(fit\\)")

  visibility <- withVisible(chains(fit))
  expect_false(visibility$visible)
})


test_that("chains() reports the documented accessor and the true dimensions", {

  skip_on_cran()

  fit <- small_pulse_fit()
  inventory <- chains(fit)

  patient <- row_for(inventory, "patient")
  expect_equal(patient$access, "patient_chain(fit)")
  expect_equal(patient$rows, nrow(patient_chain(fit)))
  expect_equal(patient$columns, ncol(patient_chain(fit)))

  pulse <- row_for(inventory, "pulse")
  expect_equal(pulse$access, "pulse_chain(fit)")
  expect_equal(pulse$rows, nrow(pulse_chain(fit)))
  expect_equal(pulse$columns, ncol(pulse_chain(fit)))
})


# Population --------------------------------------------------------------

test_that("chains() inventories a population fit, including nested pulse chains", {

  skip_on_cran()
  skip_on_ci()

  pop <- simulate_pulse_population(n_subjects = 2, num_obs = 20, seed = 1)
  fit <- fit_pulse_population(pop$data, spec = population_spec(), iters = 400,
                              thin = 10, burnin = 100, verbose = FALSE)

  inventory <- chains(fit)
  expect_setequal(inventory$chain, c("population", "subject", "pulse"))

  population <- row_for(inventory, "population")
  expect_equal(population$access, "fit$population_chain")
  expect_equal(population$rows, nrow(fit$population_chain))
  expect_equal(population$columns, ncol(fit$population_chain))

  subject <- row_for(inventory, "subject")
  expect_equal(subject$access, "fit$subject_chains[[i]]")
  expect_equal(subject$rows, nrow(fit$subject_chains[[1]]))
  expect_equal(subject$columns, ncol(fit$subject_chains[[1]]))

  # The per-subject pulse chains are a list (subject) of lists (saved draw).
  # Documenting both levels is the point: no user would guess this nesting
  # from the fit object alone.
  pulse <- row_for(inventory, "pulse")
  expect_equal(pulse$access, "fit$pulse_chains[[i]][[j]]")
  expect_equal(pulse$rows, nrow(fit$pulse_chains[[1]][[1]]))
  expect_equal(pulse$columns, ncol(fit$pulse_chains[[1]][[1]]))
})


# Joint -------------------------------------------------------------------

test_that("chains() inventories a joint fit", {

  skip_on_cran()
  skip_on_ci()

  sim <- simulate_pulse_joint(num_obs = 20, interval = 10, seed = 3)
  fit <- fit_pulse_joint(sim$driver_data, sim$response_data,
                         spec = joint_spec(), iters = 400, thin = 10,
                         burnin = 100, verbose = FALSE)

  inventory <- chains(fit)
  expect_setequal(inventory$chain,
                  c("driver", "response", "association", "driver pulse",
                    "response pulse"))

  association <- row_for(inventory, "association")
  expect_equal(association$access, "fit$association_chain")
  expect_equal(association$rows, nrow(fit$association_chain))
  expect_equal(association$columns, ncol(fit$association_chain))

  driver <- row_for(inventory, "driver")
  expect_equal(driver$access, "fit$driver_chain")
  expect_equal(driver$rows, nrow(fit$driver_chain))
})


# Input validation --------------------------------------------------------

test_that("chains() errors informatively on objects that are not model fits", {

  expect_error(chains(list(a = 1)), "pulse_fit")
  expect_error(chains(data.frame(x = 1)), "fit_pulse")
})


#------------------------------------------------------------------------------#
#    End of file                                                               #
#------------------------------------------------------------------------------#
