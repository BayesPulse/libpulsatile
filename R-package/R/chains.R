#   Extract heirarchical chains
#     - hormone1_chain() ??? driver_chain()
#     - hormone2_chain() ??? response_chain()
#     - pulse_chain()
#     - subject_chain()
#     - pop_chain()
#     - chains() --> list/str chains available and functions to extract them
#

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
#' fit   <- fit_pulse(.data = pulse, iters = 1000, thin = 10,
#'                    burnin = 100, spec = spec)
#' chains(fit)
#' pulse_chain(fit)
#' common_chain(fit)
#'
#' @export
chains <- function(fit) {
  cat("Will be a print function for showing which chain objects exist in a fit,
      how to access them, and/or some summary stats.\n")
}


#' @rdname chains
#' @export
pulse_chain <- function(fit) UseMethod("pulse_chain")

#' @export
pulse_chain.pulse_fit <- function(fit) {
  fit$pulse_chain
}

#' @rdname chains
#' @export
common_chain <- function(fit) UseMethod("common_chain")

#' @export
common_chain.pulse_fit <- function(fit) {
  fit$common_chain
}


