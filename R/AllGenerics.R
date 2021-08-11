# Definitions found in diagnostics.R -------------------------------------------

#' Creates chain trace plots
#' 
#' @param fit A \code{pulse_fit} or \code{pop_pulse_fit} object
#' @param type String indicating whether to plot as a histogram ("histogram")
#'   or density ("density").
#' @param ... Additional arguments. Use \code{patient} argument to determine
#'   which patient's trace plots are shown when \code{fit} is a
#'   \code{pop_pulse_fit} object.
#' @keywords pulse fit plot diagnostics
#' @examples
#' this_pulse <- simulate_pulse()
#' this_spec  <- pulse_spec()
#' this_fit   <- fit_pulse(data = this_pulse, iters = 1000, thin = 10,
#'                         spec = this_spec)
#' bp_trace(this_fit)
#' @export
bp_trace <- function(fit, type, ...) UseMethod("bp_trace")

#' Creates posterior density plots
#' 
#' @param fit A \code{pulse_fit} or \code{pop_pulse_fit} object
#' @param ... Additional arguments. Use \code{patient} argument to determine
#'   which patient's trace plots are shown when \code{fit} is a
#'   \code{pop_pulse_fit} object.
#' @export
bp_posteriors <- function(fit, ...) UseMethod("bp_posteriors")

#' Creates plot of posterior pulse locations
#' 
#' @param fit A \code{pulse_fit} or \code{pop_pulse_fit} object
#' @param ... Additional arguments. Use \code{patient} argument to determine
#'   which patient's trace plots are shown when \code{fit} is a
#'   \code{pop_pulse_fit} object.
#' @export
bp_location_posterior <- function(fit, ...) UseMethod("bp_location_posterior")

# Definitions found in chains.R ------------------------------------------------
#' Pulse level MCMC chain for a patient
#' 
#' @param fit A \code{pulse_fit} or \code{pop_pulse_fit} object
#' @param ... Additional arguments. Use \code{patient} argument to determine
#'   which patient's pulse chain is returned from a \code{pop_pulse_fit}
#'   object
#' @export
pulse_chain <- function(fit, ...) UseMethod("pulse_chain")

#' Patient level MCMC chain for a patient
#' 
#' @param fit A \code{pulse_fit} or \code{pop_pulse_fit} object
#' @param ... Additional arguments. Use \code{patient} argument to determine
#'   which patient's MCMC chain is returned from a \code{pop_pulse_fit}
#'   object
#' @export
patient_chain <- function(fit, ...) UseMethod("patient_chain")

#' Population level MCMC chain
#' 
#' @param fit A  \code{pop_pulse_fit} object
#' @export
population_chain <- function(fit) UseMethod("population_chain")

# Definitions found in predicted_conc.R ----------------------------------------
#' Plot predicted concentration, credible interval, and observed data
#'
#'
#' @param fit A \code{pulse_fit} or \code{pop_pulse_fit} object
#' @param predicted Result of \code{predict()} run on result of 
#'   predict.pulse_fit(fit)
#' @keywords pulse fit plot diagnostics predicted
#' @examples
#'
#' this_pulse <- simulate_pulse()
#' this_spec  <- pulse_spec(location_prior_type = "strauss",
#'                          prior_location_gamma = 0,
#'                          prior_location_range = 40)
#' this_fit   <- fit_pulse(data = this_pulse, iters = 1000, thin = 10,
#'                         spec = this_spec)
#' fit_predicted <- predict(this_fit, cred_interval = 0.8)
#' bp_predicted(this_fit, fit_predicted)
#'
#' @export
bp_predicted <- function(fit, predicted) UseMethod("bp_predicted")
