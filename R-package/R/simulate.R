
#' Simulate pulsatile hormone data
#' 
#' @description 
#'   \code{\link{simulate_pulse}} simulates a time series dataset
#'   representing blood concentration measurements of a pulsatile hormone. Both
#'   the time series and a dataset of the individual pulse characteristics are
#'   returned. 
#'    
#' @param num_obs Number of observations to simulate.  Duration of observation
#'    window equals \code{num_obs} times \code{interval}.
#' @param interval Time in minutes between observations, typically 6-10.
#' @param error_var Variance of the error added at each observation, error ~ N(0, sqrt(error_var)).
#' @param ipi_mean Mean number of sampling units between pulses (mean inter-pulse interval).
#' @param ipi_var Variance of gamma for drawing interpulse interval
#' @param ipi_min Minimum number of units between pulses
#' @param mass_mean Mean pulse mass
#' @param mass_sd Standard deviation of pulse mass
#' @param width_mean Mean pulse width (in minutes)
#' @param width_sd Standard deviation of pulse width (in minutes)
#' @param halflife_mean Mean of half-life (in minutes)
#' @param halflife_var Variance of half-life (in minutes)
#' @param baseline_mean Mean of baseline
#' @param baseline_var Variance of baseline
#' @param constant_halflife To use a constant (specified) half-life, set this
#'   to a constant value [0,inf) in minutes. Mean and variance of half-life are
#'   not used if this is non-null.
#' @param constant_baseline To use a constant (specified) baseline, set this to
#'   a constant [0,inf). Mean and variance of baseline are not used if this is
#'   non-null.
#' @return A object of class \code{pulse_sim} containing time-series dataset
#'   and dataset of characteristics of each pulse
#' @seealso print.pulse_sim, plot.pulse_sim
#' @keywords pulse simulation
#' @examples
#' this_pulse <- simulate_pulse()
#' str(this_pulse)
#' plot(this_pulse)
#' @export
simulate_pulse <- function(num_obs           = 144,
                           interval          = 10,
                           error_var         = 0.005,
                           ipi_mean          = 12,
                           ipi_var           = 40,
                           ipi_min           = 4,
                           mass_mean         = 3.5,
                           mass_sd           = 1.6,
                           width_mean        = 35,
                           width_sd          = 5,
                           halflife_mean     = NULL,
                           halflife_var      = NULL,
                           baseline_mean     = NULL,
                           baseline_var      = NULL,
                           constant_halflife = 45,
                           constant_baseline = 2.6) {

  # Add default args to function call
  args      <- formals(sys.function(sys.parent(1)))
  this_call <- match.call()
  indx      <- match(names(args), names(this_call)[-1], nomatch = 0)
  this_call <- c(as.list(this_call), args[!indx])
  this_call <- as.call(this_call)


  #---------------------------------------
  # Helper functions for drawing hormone concentration (w/o error)
  # Normal CDF 
  erfFn <- function(x) 2 * stats::pnorm(x * sqrt(2), 0, 1) - 1

  # Model concentration over time given pulse parameters  
  meanI <- function(interval, b, a, tau1, lam, s2p){
    b + (a / 2) * 
      exp((tau1 - interval) * lam + 0.5 * lam^2 * s2p) * 
      (1 + erfFn((interval - (tau1 + lam * s2p)) / sqrt(2 * s2p)))
  }

  #---------------------------------------
  # Get baseline concentration
  if (is.null(constant_baseline)) {
    while(B <= 0) B <- stats::rnorm(1, baseline_mean, sqrt(baseline_var))
  } else {
    B <- constant_baseline
  }

  #---------------------------------------
  # Get half-life of hormone
  #   H is half-life, H=ln(2)/lambda_x, where lambda_x is the decay constant
  if (is.null(constant_halflife)) {
    while(H <= 8) H <- stats::rnorm(1, halflife_mean, sqrt(halflife_var))
  } else {
    H <- constant_halflife
  }

  #---------------------------------------
  # Generate pulse locations
  #   Using a renewal process, define by interpulsatile interval and variance
  #   then convert gamma parameters
  # - mean = alpha / beta
  # - var = alpha / beta^2
  gammamean <- ipi_mean - ipi_min
  alpha     <- gammamean * gammamean / ipi_var
  beta      <- gammamean / ipi_var

  tau <- rep(0, 25)
  tau[1] <- interval * (stats::rgamma(1, shape = alpha, rate = beta) - 10)

  i <- 1
  while (tau[i] < (num_obs * interval)){
    i      <- i + 1
    tmp    <- ipi_min + stats::rgamma(1, shape = alpha, rate = beta)
    tau[i] <- tau[i-1] + (interval * tmp)
  }

  # Reduce pulse location vector to values within time range (<0, 1440)
  tau <- subset(tau, tau < (num_obs * interval))
  tau <- subset(tau, tau != 0)
  np  <- length(tau) # No. of pulses 

  #---------------------------------------
  # Generate pulse-specific parameters 
  A           <- rep(0, np)             # pulse mass
  s2p         <- rep(0, np)             # pulse width
  mass_kappa  <- rep(0, np)
  width_kappa <- rep(0, np)
  taxis       <- seq(10, (num_obs * interval), interval)
  ytmp        <- rep(0, length(taxis))  # hormone concentration

  for (i in 1:np) {
    # Log-normal 
    #A[i]   <- exp(stats::rnorm(1, mass_mean, mass_sd))
    #s2p[i] <- exp(stats::rnorm(1, width_mean, width_sd))

    # Truncated T (via gamma normal mixture)
    mass_kappa[i]  <- stats::rgamma(1, shape = 2, rate = 2)
    width_kappa[i] <- stats::rgamma(1, shape = 2, rate = 2)
    tvar  <- mass_sd^2  / mass_kappa[i]
    t2var <- width_sd^2 / width_kappa[i]
    while (A[i] < 0.25)   A[i]   <- stats::rnorm(1, mass_mean, sqrt(tvar))
    while (s2p[i] < 0.5)  s2p[i] <- stats::rnorm(1, width_mean, sqrt(t2var))
  }

  #---------------------------------------
  # Draw mean concentration
  ytmp <- 0
  for (i in 1:np) {
    ytmp   <- ytmp + meanI(taxis, 0, A[i], tau[i], log(2) / H, s2p[i])
  }
  # Add baseline and error
  ysim <- ytmp + B
  errors    <- stats::rnorm((length(taxis)), 0, sqrt(error_var))
  ysimerror <- ysim * exp(errors)

  #---------------------------------------
  # Combine into final simulated datasets
  allpulseparms <- cbind("pulse_no" = seq(1, np),
                         "mass"     = A,
                         "width"    = s2p,
                         "location" = tau,
                         "mass_kappa"  = mass_kappa,
                         "width_kappa" = width_kappa)
  allpulseparms <- tibble::as_data_frame(allpulseparms)

  ysim_df   <- cbind("observation"   = 1:length(taxis),
                     "time"          = taxis,
                     "concentration" = ysimerror)
  ysim_df   <- tibble::as_data_frame(ysim_df)

  #---------------------------------------
  # Create return object
  rtn <- structure(list("call"  = this_call,
                        "data"  = ysim_df, 
                        "parameters"  = allpulseparms),
                   class = "pulse_sim")


  return(rtn)

}





#' Plot a simulated pulsatile hormone time series.
#' 
#' @description \code{\link{plot.pulse_sim}} plots the time-series component of
#' a \code{pulse_sim} object
#'
#' @param x An object of class \code{pulse_sim} resulting from the
#'   \code{simulate_pulse()} function.
#' @param ... Other arguments not used by this method.
#' @return A ggplot2 plot of the \code{simulate_pulse} time-series dataset.
#' @seealso simulate_pulse, print.pulse_sim, summary.pulse_sim
#' @keywords pulse simulation
#' @examples
#' this_pulse <- simulate_pulse()
#' plot(this_pulse)
#' @export 
plot.pulse_sim <- function(x, ...) {

  ggplot2::ggplot(data = x$data) +
    ggplot2::aes_string(x = 'time', y = 'concentration') +
    ggplot2::geom_path()

}




# #' Summarizing simulated pulsatile hormone time series.
# #' 
# #' @description \code{\link{plot.pulse_sim}} summarizes the time-series component of
# #' a \code{pulse_sim} object
# #'
# #' @param object An object of class \code{pulse_sim} resulting from the
# #'   \code{simulate_pulse()} function.
# #' @param ... Other arguments not used by this method.
# #' @return Prints summary of the simulation
# #' @seealso simulate_pulse, print.pulse_sim, print.pulse_sim 
# #' @keywords pulse simulation
# #' @examples
# #' this_pulse <- simulate_pulse()
# #' summary(this_pulse)
# #' @export 
# summary.pulse_sim <- function(object, ...) {
# 
#   true_vals <- as.list(stats::getCall(object))
#   est_vals  <- 
#   
#   return(x)
# 
# }

# #' @export 
# print.pulse_sim <- function(x, ...) {
#   print("\n\npulse_sim object\n")
#   print("")
#   invisible(x)
# }


#-------------------------------------------------------------------------------
# End of file  
#-------------------------------------------------------------------------------
