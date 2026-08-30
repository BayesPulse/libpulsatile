# Suppress R CMD check notes for ggplot2 computed aesthetics
utils::globalVariables("density")

#' Diagnostic plots for \code{fit_pulse()} models
#'
#' Plotting functions for mcmc chains from \code{fit_pulse()} models.
#' Includes trace plots and posterior densities of the 'patient' parameters
#' and pulse location density (a set of pulse-specific parameter, from
#' 'pulse' chain).
#'
#' @param fit A model fit from \code{fit_pulse()}.
#' @param type Either histogram or density.  Only applies to
#' \code{bp_posteriors} function
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 after_stat
#' @importFrom rlang .data
#' @keywords pulse fit plot diagnostics
#' @examples
#' this_pulse <- simulate_pulse()
#' this_spec  <- pulse_spec()
#' this_fit   <- fit_pulse(data = this_pulse, iters = 1000, thin = 10,
#'                         spec = this_spec)
#' #bp_trace(this_fit)
#' @export
bp_trace <- function(fit) {

  stopifnot(inherits(fit, "pulse_fit"))

  dat <- patient_chain(fit)
  dat <- tidyr::pivot_longer(dat, cols = -"iteration",
                             names_to = "key", values_to = "value")
  ggplot2::ggplot(dat) +
  ggplot2::aes(x = .data$iteration, y = .data$value) +
  ggplot2::geom_path(linewidth = 0.15) +
  ggplot2::facet_wrap(~ key, ncol = 2, nrow = 4, scales = "free")

}


#' @rdname bp_trace
#' @export
bp_posteriors <- function(fit, type = c("histogram", "density")) {

  type <- match.arg(type)
  stopifnot(inherits(fit, "pulse_fit"))

  dat <- patient_chain(fit)
  dat <- tidyr::pivot_longer(dat, cols = -"iteration",
                             names_to = "key", values_to = "value")

  if (type == "histogram") {
    plt <-
      ggplot2::ggplot(dat) +
      ggplot2::aes(x = .data$value, y = ggplot2::after_stat(density)) +
      ggplot2::geom_histogram() +
      ggplot2::facet_wrap(~ key, ncol = 2, nrow = 4, scales = "free")
  } else if (type == "density") {
    plt <-
      ggplot2::ggplot(dat) +
      ggplot2::aes(x = .data$value) +
      ggplot2::geom_density(alpha = .2) +
      ggplot2::facet_wrap(~ key, ncol = 2, nrow = 4, scales = "free")
  }

  return(plt)

}

#' @rdname bp_trace
#' @export
bp_location_posterior <- function(fit) {

  stopifnot(inherits(fit, "pulse_fit"))
  ggplot2::ggplot(pulse_chain(fit)) +
    ggplot2::aes(x = .data$location) +
    ggplot2::geom_histogram(binwidth = 5)

}

