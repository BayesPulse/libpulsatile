#' Simulate pulsatile hormone data for a population of patients
#' 
#' @description 
#'   \code{\link{simulate_pop}} simulates multiple time series data sets
#'   representing blood concentration measurements of a pulsatile hormone. A
#'   set of time series as well as \code{data.frame}'s of population, patient,
#'   and individual pulse characteristics are returned.
#'   
#' @param num_pats Number of patients to simulate.
#' @param num_obs Number of observations to simulate per patient. Duration of 
#'    observation window equals \code{num_obs} times \code{interval}.
#' @param interval Time in minutes between observations, typically 6-10.
#' @param error_var Variance of the error added at each observation, error ~ N(0, sqrt(error_var)).
#' @param ipi_mean Mean number of sampling units between pulses (mean inter-pulse interval).
#' @param ipi_var Variance of gamma for drawing interpulse interval
#' @param ipi_min Minimum number of units between pulses
#' @param pop_mass_mean Population pulse mass
#' @param mass_s2s_sd Standard deviation of patient pulse mass means 
#' @param mass_p2p_sd Standard deviation of pulse masses
#' @param pop_width_mean Population mean pulse width (in minutes)
#' @param width_s2s_sd Standard deviation of patient pulse width means (in minutes)
#' @param width_p2p_sd Standard deviation of pulse width (in minutes)
#' @param pop_halflife_mean Population mean of half-life (in minutes)
#' @param halflife_sd Variance of patient half-lives (in minutes)
#' @param pop_baseline_mean Population mean of baseline
#' @param baseline_sd Variance of patient baselines
#' @return A object of class \code{pop_sim} containing time-series data set
#'   and data set of characteristics of each pulse
#' @seealso print.pulse_sim, plot.pulse_sim
#' @keywords pulse simulation
#' @examples
#' this_pop <- simulate_pop()
#' str(this_pop)
#' plot(this_pop)
#' @importFrom truncnorm rtruncnorm
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @export
simulate_pop <- function(num_pats = 10,
                         num_obs = 144,
                         interval = 10,
                         error_var = 0.005,
                         ipi_mean = 12,
                         ipi_var = 40,
                         ipi_min = 4,
                         pop_mass_mean = 3.5,
                         mass_s2s_sd = 1,
                         mass_p2p_sd = 1.6,
                         pop_width_mean = 35,
                         width_s2s_sd = 1,
                         width_p2p_sd = 5,
                         pop_halflife_mean = 45,
                         halflife_sd = 5,
                         pop_baseline_mean = 2.5,
                         baseline_sd = 0.75) {
  
  # Add default args to function call
  args      <- formals(sys.function(sys.parent(1)))
  this_call <- match.call()
  indx      <- match(names(args), names(this_call)[-1], nomatch = 0)
  this_call <- c(as.list(this_call), args[!indx])
  this_call <- as.call(this_call)
  
  # Draw subject level mean mass, mean width, baseline, and halflife from
  # truncated normal distributions
  pulseMass <- rtruncnorm(num_pats, a = 0, b = Inf, 
                          mean = pop_mass_mean, sd = mass_s2s_sd)
  pulseWidth <-rtruncnorm(num_pats, a = 0, b = Inf, 
                          mean = pop_width_mean, sd = width_s2s_sd)
  baseline <- rtruncnorm(num_pats, a = 0, b = Inf, mean = 
                           pop_baseline_mean, sd = baseline_sd)
  halflife <- rtruncnorm(num_pats, a = 0, b = Inf, mean = 
                           pop_halflife_mean, sd = halflife_sd)
  
  # Create objects to store relevant data
  pulseTruth <- data.frame(pat_id = numeric(), pulse_no = numeric(), 
                           mass = numeric(), width = numeric(), 
                           location = numeric(), mass_kappa = numeric(),
                           width_kappa = numeric())
  simData <- list()
  
  # Loop to simulate data for patients
  for (i in 1:num_pats) {
    patData <- simulate_pulse(num_obs = num_obs, interval = interval, 
                              error_var = error_var,
                              ipi_mean = ipi_mean, ipi_var = ipi_var, 
                              ipi_min = ipi_min,
                              mass_mean = pulseMass[i], 
                              mass_sd = mass_p2p_sd,
                              width_mean = pulseWidth[i], 
                              width_sd = width_p2p_sd,
                              constant_halflife = halflife[i],
                              constant_baseline = baseline[i])
    simData[[i]] <- patData$data
    patParams <- cbind(pat_id = rep(i, nrow(patData$parameters)), 
                       patData$parameters)
    pulseTruth <- rbind(pulseTruth, patParams)
  }
  
  # Create DF of subject level values
  patTruth <- pulseTruth %>%
    group_by(.data$pat_id) %>%
    summarise(pulseCount = n(),
              massMean = mean(.data$mass),
              widthMean = mean(.data$width)) %>%
    mutate(halflife = halflife,
           baseline = baseline)
  
  # Create DF of population level values
  populationTruth <- patTruth %>%
    summarise(subjectCount = n(),
              pulseCountMean = mean(.data$pulseCount),
              massMean = mean(.data$massMean),
              widthMean = mean(.data$widthMean),
              halflifeMean = mean(.data$halflife),
              baselineMean = mean(.data$baseline))
  
  # Create structure to be returned
  names(simData) <- paste0("Patient_", 1:num_pats)
  rtn <- structure(list("call"  = this_call,
                        "data"  = simData, 
                        "pop_params" = populationTruth,
                        "patient_params"  = patTruth,
                        "pulse_params" = pulseTruth),
                   class = "pop_sim")
  
  return(rtn)
  
}

#' Plot a simulated population of pulsatile hormone time series.
#' 
#' @description \code{\link{plot.pop_sim}} plots the time-series components of
#' a \code{pop_sim} object
#'
#' @param x An object of class \code{pulse_sim} resulting from the
#'   \code{simulate_pulse()} function.
#' @param ... Other arguments not used by this method.
#' @return A ggplot2 plot of the \code{simulate_pulse} time-series dataset.
#' @seealso simulate_pop, print.pop_sim, summary.pop_sim
#' @keywords pulse simulation
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_line
#' @examples
#' this_pop <- simulate_pop()
#' plot(this_pop)
#' @export 
plot.pop_sim <- function(x, ...) {
  
  # Create data to plot
  if(x$call$num_pats > 12) {
    warning("Only plotting first 12 patients")
    plotData <- bind_rows(x$data, .id = "pat_id")[1:(x$call$num_obs * 12), ]
  } else {
    plotData <- bind_rows(x$data, .id = "pat_id")
  }
  
  # Plot data
  ggplot(plotData) +
    aes_string(x = "time", y = "concentration") +
    geom_line() +
    ggplot2::facet_wrap( ~ pat_id, ncol = 3, nrow = 4, scales = "fixed") +
    labs(title = "Simulated Population Hormone Concentrations",
         x = "Time", y = "Concentration") +
    theme(plot.title = element_text(hjust = 0.5))
  
}