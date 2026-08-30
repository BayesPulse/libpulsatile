#' Predicted concentration and related functions
#'
#' Calculates the predicted concentration and specified credible interval. The
#' plot function plots this predicted concentration with 80% credible intervals
#' (adjustable by argument) and the true concentration.
#'
#' @param object A model fit from \code{fit_pulse()} (class "pulse_object").
#' @param cred_interval Size of the credible interval to calculate.
#' @param ... further arguments passed to or from other methods.
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unnest
#' @importFrom dplyr select
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom rlang sym
#' @importFrom rlang .data
#' @importFrom dplyr full_join
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr arrange
#' @importFrom stats quantile
#' @importFrom stats median
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
#'
#' @export
predict.pulse_fit <- function(object, cred_interval = 0.8, ...) {

  stopifnot(inherits(object, "pulse_fit"))

  if (inherits(object$data, "pulse_sim")) {
    data <- object$data$data
  } else {
    data <- object$data
  }

  diff_from_bound <- (1 - cred_interval) / 2
  lwr             <- diff_from_bound
  upr             <- 1 - diff_from_bound

  onechain <- full_join(pulse_chain(object), patient_chain(object), 
                        by = c("iteration", "total_num_pulses" = "num_pulses"))

  # Only plot within range of data -- saw some werid behavior using full
  # fitstart/fitend
  #   fitstart <- object$time_range[1]
  #   fitend   <- object$time_range[2]
  fitstart <- min(data$time)
  fitend   <- max(data$time)
  time     <- seq(fitstart, fitend, by = 10) # TODO: change 10 to whatever the sampling interval is
  fitstart <- sym(as.character(fitstart))
  fitend <- sym(as.character(fitend))

  onechain <- onechain %>%
    mutate(time = list(time)) %>%
    unnest(cols = c(time)) %>%
    mutate(mean_contrib = calc_mean_contrib(time, .data$location,
                                            .data$halflife, .data$mass,
                                            .data$width))
  onechain <-
    onechain %>%
    select(.data$iteration, .data$pulse_num, .data$baseline, .data$model_error,
           .data$time, .data$mean_contrib)

  wide <- onechain %>%
    pivot_wider(names_from = .data$time, values_from = .data$mean_contrib)

  # Get numeric column names for the time columns
  time_cols <- setdiff(names(wide), c("iteration", "pulse_num", "baseline", "model_error"))

  long <- wide %>%
    group_by(.data$iteration, .data$baseline, .data$model_error) %>%
    summarise(across(all_of(time_cols), sum), .groups = "drop") %>%
    mutate(across(all_of(time_cols),
                  ~ add_baseline_error(.x, baseline = .data$baseline,
                                       model_error = .data$model_error)))

  long <-
    long %>%
    pivot_longer(cols = all_of(time_cols), names_to = "time",
                 values_to = "concentration") %>%
    mutate(time = as.numeric(.data$time))

  rtn <- long %>% 
    ungroup %>% 
    arrange(.data$iteration, .data$time) %>% 
    select(-.data$baseline, -.data$model_error) %>%
    group_by(time) %>%
    summarise(mean.conc = mean(.data$concentration),
              med.conc  = median(.data$concentration),
              upper = quantile(.data$concentration, upr),
              lower = quantile(.data$concentration, lwr))
  return(rtn)

}




#' Creates the integral in the deconvolution part of the model
#'
#' @param time Some time variable
#' @importFrom stats pnorm
erf <- function(time){
  y <- 2 * pnorm(time * sqrt(2), 0, 1) - 1
  y
}


#' Creates the expected hormone concentration for a particular pulse in a
#' particular iteration
#'
#' @param time Vector of times from fitstart to fitend, by minute? sampling
#' interval?
#' @param location Scalar. Location of this pulse
#' @param halflife Scalar. Halflife for this iteration.
#' @param mass Scalar. Mass of this pulse.
#' @param width Scalar. Width of this pulse
#' @rdname predict.pulse_fit
calc_mean_contrib <- function(time, location, halflife, mass, width) {

  decayrate <- log(2)/halflife

  mean_conc <- (mass / 2) *
    pmin(exp((location - time) * decayrate + 0.5 * decayrate^2 * width), 10^100) *
    (1 + erf((time - (location + decayrate * width)) / sqrt(2 * width)))

  return(mean_conc)

}

#' Adds baseline and model error to the mean concentration vector for a given
#' iteration (which is just the elementwise sum of mean_contribution vectors).
#'
#' @param x Mean concentration (elementwise sum of mean_contribution vectors).
#' @param baseline Baseline concentration value for this iteration.
#' @param model_error Model error estimate for this iteration.
#' @importFrom stats rnorm
#' @rdname predict.pulse_fit
add_baseline_error <- function(x, baseline, model_error) {
  x + baseline + rnorm(1, 0, sqrt(model_error))
}





### TODO: so far haven't added eta/tvarscale_ to the prediction functions --
### Don't think you should, this only comes into play in the random draw of the
### SD?

### TODO: add option of passing just the fit object to bp_predicted()


#' Plot predicted concentration, credible interval, and observed data
#'
#'
#' @param fit a fit_pulse object
#' @param predicted result of predict.pulse_fit(fit)
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
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
bp_predicted <- function(fit, predicted) {

  if (inherits(fit$data, "pulse_sim")) {
    data <- fit$data$data
  } else {
    data <- fit$data
  }

  ggplot(data) +
    aes(x = .data$time, y = .data$concentration) +
    geom_path() +
    geom_point() +
    geom_path(data = predicted, aes(y = .data$mean.conc), color = "red") +
    geom_path(data = predicted, aes(y = .data$upper), color = "red", linetype = "dashed") +
    geom_path(data = predicted, aes(y = .data$lower), color = "red", linetype = "dashed") +
    xlab("Time (minutes)") +
    ylab("Concentration (ng/mL)")

}



