#' Diagnostic plots for \code{fit_pulse()} models
#' 
#' Plotting functions for mcmc chains from \code{fit_pulse()} models.  Includes
#' trace plots and posterior densities of the 'common' parameters and pulse
#' location density (a set of pulse-specific parameter, from 'pulse' chain).
#' 
#' @useDynLib pulsatile decon_r_interface
#' @param fit A model fit from \code{fit_pulse()}.
#' @param type Either histogram or density.  Only applies to
#' \code{bp_posteriors} function
#' @import tidyr dplyr ggplot2
#' @keywords pulse fit plot diagnostics
#' @examples
#' this_pulse <- simulate_pulse()
#' this_spec  <- pulse_spec()
#' this_fit   <- fit_pulse(.data = this_pulse, iters = 1000, thin = 10,
#'                         spec = this_spec)
#' bp_trace(this_fit)
#' @export
bp_trace <- function(fit) {

  stopifnot(class(fit) == "pulse_fit")

  dat <- common_chain(fit) 
  dat <- tidyr::gather_(dat, key = "key", value = "value",
                        dplyr::select_vars_(names(dat), names(dat),
                                            exclude = "iteration"))
  ggplot2::ggplot(dat) +
  ggplot2::aes_string(x = "iteration", y = "value") +
  ggplot2::geom_path(size = 0.15) +
  ggplot2::facet_wrap( ~ key, ncol = 2, nrow = 4, scales = "free")

}


#' @rdname bp_trace
#' @export
bp_posteriors <- function(fit, type = c("histogram", "density")) {

  stopifnot(class(fit) == "pulse_fit")

  dat <- common_chain(fit) 
  dat <- tidyr::gather_(dat, key = "key", value = "value", 
                        dplyr::select_vars_(names(dat), names(dat),
                                            exclude = "iteration"))

  if (type == "histogram") {
    plt <- 
      ggplot2::ggplot(dat) +
      ggplot2::aes_string(x = "value", y = "..density..") +
      ggplot2::geom_histogram(size = 0.15) + #aes(y = ..density..), size = 0.15) +
      ggplot2::facet_wrap( ~ key, ncol = 2, nrow = 4, scales = "free")
  } else if (type == "density") {
    plt <-
      ggplot2::ggplot(dat) +
      ggplot2::aes_string(x = "value") +
      ggplot2::geom_density(alpha = .2) +
      ggplot2::facet_wrap( ~ key, ncol = 2, nrow = 4, scales = "free")
  }

  return(plt)

}

#' @rdname bp_trace
#' @export
bp_location_posterior <- function(fit) {

  stopifnot(class(fit) == "pulse_fit")
  pulse_chain(fit) %>%
    ggplot(aes(x = location)) +
    geom_histogram(binwidth = 5) 

}

