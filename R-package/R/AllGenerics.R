# Definitions found in diagnostics.R -------------------------------------------

#' Creates chain trace plots
#' @export
bp_trace <- function(fit, ...) UseMethod("bp_trace")

#' Creates chain trace plots
#' @export
bp_posteriors <- function(fit, ...) UseMethod("bp_posteriors")

#' Creates plot of posterior pulse locations
#' @export
bp_location_posterior <- function(fit, ...) UseMethod("bp_location_posterior")

# Definitions found in chains.R ------------------------------------------------

#' @rdname chains
#' @export
pulse_chain <- function(fit, ...) UseMethod("pulse_chain")

#' @rdname chains
#' @export
patient_chain <- function(fit, ...) UseMethod("patient_chain")

#' @rdname chains
#' @export
population_chain <- function(fit) UseMethod("population_chain")

# Definitions found in predicted_conc.R ----------------------------------------
#' @export
bp_predicted <- function(fit, ...) UseMethod("bp_predicted")
