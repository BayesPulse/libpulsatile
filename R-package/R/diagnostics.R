#' Diagnostic plots for \code{fit_pulse()} models
#' 
#' Plotting functions for mcmc chains from \code{fit_pulse()} models.  Includes
#' trace plots and posterior densities of the 'patient' parameters and pulse
#' location density (a set of pulse-specific parameter, from 'pulse' chain).
#' 
#' @param fit A model fit from \code{fit_pulse()}.
#' @param type Either histogram or density.  Only applies to
#' \code{bp_posteriors} function
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr gather_
#' @importFrom dplyr select_vars_
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 facet_wrap
#' @keywords pulse fit plot diagnostics
#' @examples
#' this_pulse <- simulate_pulse()
#' this_spec  <- pulse_spec()
#' this_fit   <- fit_pulse(data = this_pulse, iters = 1000, thin = 10,
#'                         spec = this_spec)
#' #bp_trace(this_fit)
#' @export
bp_trace.pulse_fit <- function(fit) {

  stopifnot(class(fit) == "pulse_fit")

  # dat <- patient_chain(fit)
  # dat <- tidyr::gather_(dat, key = "key", value = "value",
  #                       dplyr::select_vars_(names(dat), names(dat),
  #                                           exclude = "iteration"))
  # ggplot2::ggplot(dat) +
  # ggplot2::aes_string(x = "iteration", y = "value") +
  # ggplot2::geom_path(size = 0.15) +
  # ggplot2::facet_wrap( ~ key, ncol = 2, nrow = 4, scales = "free")
  
  dat <- patient_chain(fit)
  dat <- tidyr::pivot_longer(dat, !iteration, names_to = "parameter",
                             values_to = "value")
  ggplot(dat) +
    aes_string(x = "iteration", y = "value") +
    geom_path(size = 0.15) +
    facet_wrap(~ parameter, ncol = 2, nrow = 4, scales = "free")

}

#' @importFrom tidyr bind_rows
#' @export
bp_trace.pop_pulse_fit <- function(fit, patient = NULL) {
  
  dat <- population_chain(fit)
  colnames(dat)[2:11] <- paste0("pop_", colnames(dat)[2:11])
  dat <- tidyr::pivot_longer(dat, !iteration, names_to = "parameter",
                                values_to = "value")
  
  if (!is.null(patient)) {
    patDat <- patient_chain(fit, patient)
    colnames(patDat)[2:6] <- paste0("pat_", colnames(patDat)[2:6])
    patDat <- tidyr::pivot_longer(patDat, !iteration, names_to = "parameter",
                               values_to = "value")
    
    dat <- dplyr::bind_rows(dat, patDat)
    
    nCol <- 5
    nRow <- 3
  } else {
    nCol <- 4
    nRow <- 3
  }
  
  ggplot(dat) +
    aes_string(x = "iteration", y = "value") +
    geom_path(size = 0.15) +
    facet_wrap(~ parameter, ncol = nCol, nrow = nRow, scales = "free")
}

#' @export
bp_trace <- function(fit, ...) UseMethod("bp_trace")


#' @rdname bp_trace
#' @export
bp_posteriors <- function(fit, type = c("histogram", "density")) {

  stopifnot(class(fit) == "pulse_fit")

  dat <- patient_chain(fit) 
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
  ggplot2::ggplot(pulse_chain(fit)) +
    ggplot2::aes_string(x = "location") +
    ggplot2::geom_histogram(binwidth = 5)

}

