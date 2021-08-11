#' chains
#' 
#' Functions for extracting and working with the various levels of posterior
#' mcmc chains.
#' 
#' @param fit A model fit from \code{fit_pulse}.
#' @import tibble
#' @keywords pulse fit
#' @examples
#' 
#' pulse <- simulate_pulse()
#' spec  <- pulse_spec()
#' fit   <- fit_pulse(data = pulse, iters = 1000, thin = 10,
#'                    burnin = 100, spec = spec)
#' #chains(fit)
#' #pulse_chain(fit)
#' #patient_chain(fit)
#'
#' @export
chains <- function(fit) {
  cat("Will be a print function for showing which chain objects exist in a fit,
      how to access them, and/or some summary stats.\n")
}

#' @export
population_chain.pop_pulse_fit <- function(fit) {
  fit$population_chain
}

#' @export
pulse_chain.pulse_fit <- function(fit, ...) {
  fit$pulse_chain
}

#' @export
pulse_chain.pop_pulse_fit <- function(fit, ...) {
  args <- list(...)
  if(is.null(args$patient)) { args$patient <- 1}
  fit$pulse_chain[[args$patient]]
}

#' @export
patient_chain.pulse_fit <- function(fit, ...) {
  fit$patient_chain
}

#' @export
patient_chain.pop_pulse_fit <- function(fit, ...) {
  args <- list(...)
  if(is.null(args$patient)) { args$patient <- 1}
  fit$patient_chain[[args$patient]]
}


