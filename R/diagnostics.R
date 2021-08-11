#' @importFrom tidyr pivot_longer
#' @importFrom tidyr gather_
#' @importFrom dplyr select_vars_
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 facet_wrap
#' @export
bp_trace.pulse_fit <- function(fit, ...) {

  stopifnot(class(fit) == "pulse_fit")
  
  dat <- patient_chain(fit)
  dat <- tidyr::pivot_longer(dat, !.data$iteration, names_to = "parameter",
                             values_to = "value")
  ggplot(dat) +
    aes_string(x = "iteration", y = "value") +
    geom_path(size = 0.15) +
    facet_wrap(~ parameter, ncol = 2, nrow = 4, scales = "free")

}

#' @importFrom dplyr bind_rows
#' @export
bp_trace.pop_pulse_fit <- function(fit, ...) {
  args <- list(...)
  
  dat <- population_chain(fit)
  colnames(dat)[2:11] <- paste0("pop_", colnames(dat)[2:11])
  dat <- tidyr::pivot_longer(dat, !.data$iteration, names_to = "parameter",
                                values_to = "value")
  
  if (!is.null(args$patient)) {
    patDat <- patient_chain(fit, args$patient)
    colnames(patDat)[2:6] <- paste0("pat_", colnames(patDat)[2:6])
    patDat <- tidyr::pivot_longer(patDat, !.data$iteration, names_to = "parameter",
                               values_to = "value")
    
    dat <- dplyr::bind_rows(dat, patDat)
    
    nCol <- 5
    nRow <- 4
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
bp_posteriors.pulse_fit <- function(fit, type = c("histogram", "density"), ...) {

  type <- match.arg(type)
  dat <- patient_chain(fit)
  dat <- tidyr::pivot_longer(dat, !.data$iteration, names_to = "parameter",
                                values_to = "value")

  if (type == "histogram") {
    plt <- 
      ggplot2::ggplot(dat) +
      ggplot2::aes_string(x = "value", y = "..density..") +
      ggplot2::geom_histogram(size = 0.15) + #aes(y = ..density..), size = 0.15) +
      ggplot2::facet_wrap( ~ parameter, ncol = 2, nrow = 4, scales = "free")
  } else if (type == "density") {
    plt <-
      ggplot2::ggplot(dat) +
      ggplot2::aes_string(x = "value") +
      ggplot2::geom_density(alpha = .2) +
      ggplot2::facet_wrap( ~ parameter, ncol = 2, nrow = 4, scales = "free")
  }

  return(plt)

}

#' @export
bp_posteriors.pop_pulse_fit <- function(fit, type = c("histogram", "density"), 
                                        ...) {
  
  args <- list(...)
  type <- match.arg(type)
  dat <- population_chain(fit)
  colnames(dat)[2:11] <- paste0("pop_", colnames(dat)[2:11])
  dat <- tidyr::pivot_longer(dat, !.data$iteration, names_to = "parameter",
                             values_to = "value")
  
  if (!is.null(args$patient)) {
    patDat <- patient_chain(fit, args$patient)
    colnames(patDat)[2:6] <- paste0("pat_", colnames(patDat)[2:6])
    patDat <- tidyr::pivot_longer(patDat, !.data$iteration, names_to = "parameter",
                                  values_to = "value")
    
    dat <- dplyr::bind_rows(dat, patDat)
    
    nCol <- 5
    nRow <- 4
  } else {
    nCol <- 4
    nRow <- 3
  }
  
  if (type == "histogram") {
    plt <- 
      ggplot2::ggplot(dat) +
      ggplot2::aes_string(x = "value", y = "..density..") +
      ggplot2::geom_histogram(size = 0.15) + #aes(y = ..density..), size = 0.15) +
      ggplot2::facet_wrap( ~ parameter, ncol = nCol, nrow = nRow, scales = "free")
  } else if (type == "density") {
    plt <-
      ggplot2::ggplot(dat) +
      ggplot2::aes_string(x = "value") +
      ggplot2::geom_density(alpha = .2) +
      ggplot2::facet_wrap( ~ parameter, ncol = nCol, nrow = nRow, scales = "free")
  }
  
  return(plt)
}

#' @export
bp_location_posterior.pulse_fit <- function(fit, ...) {
  ggplot2::ggplot(pulse_chain(fit)) +
    ggplot2::aes_string(x = "location") +
    ggplot2::geom_histogram(binwidth = 5)

}

#' @export
bp_location_posterior.pop_pulse_fit <- function(fit, ...) {
  args <- list(...)
  if(!is.null(args$patient)) {
    rtn <- ggplot2::ggplot(pulse_chain(fit, args$patient)) +
      ggplot2::aes_string(x = "location") +
      ggplot2::geom_histogram(binwidth = 5)
  } else {
    location_data <- fit$pulse_chain %>%
      bind_rows(.id = "patient")
    
    rtn <- ggplot(location_data) +
      aes_string(x = "location") +
      geom_histogram(binwidth = 5) +
      facet_wrap(~ patient)
  }
  
  return(rtn)
  
}


