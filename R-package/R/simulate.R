
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
#' @param mass_mean Mean pulse mass. If \code{NULL} (default) it resolves to the
#'   parameterization-appropriate default: natural-scale \code{3.5} under
#'   \code{pulse_distribution = "truncnorm"} or LOG-scale \code{1.2} under
#'   \code{"lognormal"}.
#' @param mass_sd Standard deviation of pulse mass. If \code{NULL} (default),
#'   resolves to \code{1.6} (truncnorm) or LOG-scale \code{0.5} (lognormal).
#' @param width_mean Mean pulse width. If \code{NULL} (default), resolves to
#'   natural-scale \code{35} (truncnorm) or LOG-scale \code{3.0} (lognormal).
#' @param width_sd Standard deviation of pulse width. If \code{NULL} (default),
#'   resolves to \code{5} (truncnorm) or LOG-scale \code{0.7} (lognormal).
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
#' @param pulse_distribution Character, one of \code{"truncnorm"} (default) or
#'   \code{"lognormal"}, selecting the data-generating distribution for the
#'   per-pulse mass and width random effects. \code{"truncnorm"} keeps the legacy
#'   natural-scale truncated Student-t mixture (a Gaussian/gamma scale mixture
#'   with rejection below the physical floors), in which \code{mass_mean},
#'   \code{mass_sd}, \code{width_mean} and \code{width_sd} are natural-scale
#'   parameters. \code{"lognormal"} draws \code{mass = exp(rnorm(1, mass_mean,
#'   mass_sd))} and \code{width = exp(rnorm(1, width_mean, width_sd))} (Gaussian
#'   on the log scale, no t-scale mixture and no truncation), so \code{mass_mean},
#'   \code{mass_sd}, \code{width_mean} and \code{width_sd} are interpreted as
#'   LOG-scale parameters. Use \code{"lognormal"} to generate data consistent
#'   with the papers-default (\code{pulse_distribution = "lognormal"}) fit
#'   produced by \code{\link{pulse_spec}} / \code{\link{population_spec}}.
#' @note This function defaults to \code{pulse_distribution = "truncnorm"},
#'   whereas \code{\link{pulse_spec}} defaults to \code{"lognormal"}. Simulating
#'   with defaults and fitting with defaults therefore fits a lognormal model to
#'   truncnorm-generated data. Pass \code{pulse_distribution = "lognormal"} to
#'   match the papers-default fit.
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
                           mass_mean         = NULL,
                           mass_sd           = NULL,
                           width_mean        = NULL,
                           width_sd          = NULL,
                           halflife_mean     = NULL,
                           halflife_var      = NULL,
                           baseline_mean     = NULL,
                           baseline_var      = NULL,
                           constant_halflife = 45,
                           constant_baseline = 2.6,
                           pulse_distribution = c("truncnorm", "lognormal")) {

  pulse_distribution <- match.arg(pulse_distribution)

  # Conditional defaults for the pulse mass/width location-scale parameters.
  # Under "truncnorm" these are natural-scale (legacy defaults, byte-for-byte);
  # under "lognormal" they are LOG-scale, so exp() of them yields sane natural
  # widths/masses instead of exp(rnorm(35, 5)) ~ 1e15. A user-supplied value
  # always wins. Mirrors the NULL-sentinel pattern in pulse_spec().
  if (pulse_distribution == "lognormal") {
    if (is.null(mass_mean))  mass_mean  <- 1.2
    if (is.null(mass_sd))    mass_sd    <- 0.5
    if (is.null(width_mean)) width_mean <- 3.0
    if (is.null(width_sd))   width_sd   <- 0.7
  } else {
    if (is.null(mass_mean))  mass_mean  <- 3.5
    if (is.null(mass_sd))    mass_sd    <- 1.6
    if (is.null(width_mean)) width_mean <- 35
    if (is.null(width_sd))   width_sd   <- 5
  }

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
    if (pulse_distribution == "lognormal") {
      # Log-normal random effects (Gaussian on the log scale, no t-scale
      # mixture, no truncation). mass_mean/mass_sd/width_mean/width_sd are
      # LOG-scale parameters, matching the papers-default lognormal fit.
      mass_kappa[i]  <- 1
      width_kappa[i] <- 1
      A[i]   <- exp(stats::rnorm(1, mass_mean, mass_sd))
      s2p[i] <- exp(stats::rnorm(1, width_mean, width_sd))
    } else {
      # Truncated T (via gamma normal mixture)
      mass_kappa[i]  <- stats::rgamma(1, shape = 2, rate = 2)
      width_kappa[i] <- stats::rgamma(1, shape = 2, rate = 2)
      tvar  <- mass_sd^2  / mass_kappa[i]
      t2var <- width_sd^2 / width_kappa[i]
      while (A[i] < 0.25)   A[i]   <- stats::rnorm(1, mass_mean, sqrt(tvar))
      while (s2p[i] < 0.5)  s2p[i] <- stats::rnorm(1, width_mean, sqrt(t2var))
    }
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
  allpulseparms <- tibble::as_tibble(allpulseparms)

  ysim_df   <- cbind("observation"   = 1:length(taxis),
                     "time"          = taxis,
                     "concentration" = ysimerror)
  ysim_df   <- tibble::as_tibble(ysim_df)

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


#' Simulate coupled driver-response pulsatile hormone data
#'
#' @description \code{simulate_pulse_joint} simulates a pair of hormone time
#'   series for the joint driver-response model. The driver hormone is generated
#'   as a standard pulsatile series. Response-hormone pulse locations are then
#'   drawn from a non-homogeneous process whose intensity is modulated by the
#'   driver pulses through the coupling kernel, matching the model's response
#'   birth rate \code{base_rate * (1 + lambda(t))}, where
#'   \code{lambda(t) = sum_j (rho / sqrt(2*pi*nu)) * exp(-(t - tau_j)^2 / (2*nu))}
#'   and \code{tau_j} are the driver pulse locations. Larger \code{rho} produces
#'   more response pulses clustered near driver pulses; smaller \code{nu}
#'   produces tighter clustering.
#'
#' @param rho Coupling strength (cluster size); must be > 0.
#' @param nu Coupling temporal spread (cluster width, a variance in minutes^2);
#'   must be > 0.
#' @param num_obs Number of observations. Window duration is
#'   \code{num_obs * interval}.
#' @param interval Time in minutes between observations.
#' @param response_pulse_count Base (driver-independent) expected number of
#'   response pulses over the window.
#' @param driver_mass_mean,driver_mass_sd Driver pulse mass mean and SD. If
#'   \code{NULL} (default) they resolve to natural-scale \code{3.5}/\code{1.0}
#'   under \code{pulse_distribution = "truncnorm"} or LOG-scale
#'   \code{1.2}/\code{0.5} under \code{"lognormal"}.
#' @param driver_width_mean,driver_width_sd Driver pulse width mean and SD. If
#'   \code{NULL} (default) they resolve to natural-scale \code{35}/\code{5}
#'   (truncnorm) or LOG-scale \code{3.0}/\code{0.7} (lognormal).
#' @param driver_baseline,driver_halflife Driver baseline and half-life.
#' @param driver_ipi_mean,driver_ipi_var Driver inter-pulse interval mean and
#'   variance.
#' @param response_mass_mean,response_mass_sd Response pulse mass mean and SD,
#'   with the same \code{NULL} defaults as the driver equivalents
#'   (\code{3.5}/\code{1.0} truncnorm, LOG-scale \code{1.2}/\code{0.5} lognormal).
#' @param response_width_mean,response_width_sd Response pulse width mean and SD,
#'   with the same \code{NULL} defaults as the driver equivalents
#'   (\code{35}/\code{5} truncnorm, LOG-scale \code{3.0}/\code{0.7} lognormal).
#' @param response_baseline,response_halflife Response baseline and half-life.
#' @param error_var Variance of the multiplicative log-normal observation error.
#' @param pulse_distribution Character, one of \code{"truncnorm"} (default) or
#'   \code{"lognormal"}, selecting the data-generating distribution for the
#'   per-pulse mass and width random effects of BOTH the driver and response
#'   hormones. \code{"truncnorm"} keeps the legacy natural-scale truncated
#'   Student-t mixture, in which the \code{*_mass_mean}, \code{*_mass_sd},
#'   \code{*_width_mean} and \code{*_width_sd} arguments are natural-scale.
#'   \code{"lognormal"} draws \code{mass = exp(rnorm(1, mass_mean, mass_sd))} and
#'   \code{width = exp(rnorm(1, width_mean, width_sd))} (Gaussian on the log
#'   scale, no t-scale mixture and no truncation), so those arguments are
#'   LOG-scale. Use \code{"lognormal"} to generate data consistent with the
#'   papers-default (\code{pulse_distribution = "lognormal"}) fit produced by
#'   \code{\link{joint_spec}}.
#' @param seed Optional RNG seed for reproducibility.
#' @note The data-generating default here is \code{pulse_distribution =
#'   "truncnorm"}, whereas \code{\link{joint_spec}} (and thus
#'   \code{\link{fit_pulse_joint}}) defaults to \code{"lognormal"}. Simulating
#'   with defaults and fitting with defaults therefore fits a lognormal model to
#'   truncnorm-generated data. Pass \code{pulse_distribution = "lognormal"} here
#'   to match the papers-default fit.
#' @return An object of class \code{joint_sim}: a list with \code{driver_data}
#'   and \code{response_data} (each a tibble with \code{time} and
#'   \code{concentration}), the \code{driver_parameters} and
#'   \code{response_parameters} pulse tables, and \code{association} (the true
#'   \code{rho}, \code{nu}, and \code{response_pulse_count}).
#' @seealso \code{\link{simulate_pulse}}, \code{\link{fit_pulse_joint}}
#' @keywords pulse simulation
#' @examples
#' sim <- simulate_pulse_joint(rho = 0.8, nu = 100, seed = 1)
#' str(sim$driver_data)
#' @importFrom stats rgamma rnorm runif rpois pnorm
#' @export
simulate_pulse_joint <- function(rho                  = 0.5,
                                 nu                   = 100,
                                 num_obs              = 144,
                                 interval             = 10,
                                 response_pulse_count = 12,
                                 driver_mass_mean     = NULL,
                                 driver_mass_sd       = NULL,
                                 driver_width_mean    = NULL,
                                 driver_width_sd      = NULL,
                                 driver_baseline      = 2.6,
                                 driver_halflife      = 45,
                                 driver_ipi_mean      = 12,
                                 driver_ipi_var       = 40,
                                 response_mass_mean   = NULL,
                                 response_mass_sd     = NULL,
                                 response_width_mean  = NULL,
                                 response_width_sd    = NULL,
                                 response_baseline    = 2.6,
                                 response_halflife    = 45,
                                 error_var            = 0.005,
                                 pulse_distribution   = c("truncnorm",
                                                          "lognormal"),
                                 seed                 = NULL) {

  pulse_distribution <- match.arg(pulse_distribution)

  # Conditional defaults for the driver and response pulse mass/width
  # location-scale parameters. Under "truncnorm" these are natural-scale (legacy
  # defaults, byte-for-byte); under "lognormal" they are LOG-scale, so exp() of
  # them yields sane natural masses/widths (matching pulse_spec()'s log-scale
  # defaults) instead of exp(rnorm(35, 5)) ~ 1e15. A user-supplied value always
  # wins. Mirrors the NULL-sentinel pattern in simulate_pulse().
  if (pulse_distribution == "lognormal") {
    if (is.null(driver_mass_mean))    driver_mass_mean    <- 1.2
    if (is.null(driver_mass_sd))      driver_mass_sd      <- 0.5
    if (is.null(driver_width_mean))   driver_width_mean   <- 3.0
    if (is.null(driver_width_sd))     driver_width_sd     <- 0.7
    if (is.null(response_mass_mean))  response_mass_mean  <- 1.2
    if (is.null(response_mass_sd))    response_mass_sd    <- 0.5
    if (is.null(response_width_mean)) response_width_mean <- 3.0
    if (is.null(response_width_sd))   response_width_sd   <- 0.7
  } else {
    if (is.null(driver_mass_mean))    driver_mass_mean    <- 3.5
    if (is.null(driver_mass_sd))      driver_mass_sd      <- 1.0
    if (is.null(driver_width_mean))   driver_width_mean   <- 35
    if (is.null(driver_width_sd))     driver_width_sd     <- 5
    if (is.null(response_mass_mean))  response_mass_mean  <- 3.5
    if (is.null(response_mass_sd))    response_mass_sd    <- 1.0
    if (is.null(response_width_mean)) response_width_mean <- 35
    if (is.null(response_width_sd))   response_width_sd   <- 5
  }

  if (!is.null(seed)) set.seed(seed)
  stopifnot(rho > 0, nu > 0, response_pulse_count > 0)

  T_end <- num_obs * interval

  # Concentration helpers (mirror simulate_pulse)
  erfFn <- function(x) 2 * stats::pnorm(x * sqrt(2), 0, 1) - 1
  meanI <- function(tt, b, a, tau1, lam, s2p) {
    b + (a / 2) *
      exp((tau1 - tt) * lam + 0.5 * lam^2 * s2p) *
      (1 + erfFn((tt - (tau1 + lam * s2p)) / sqrt(2 * s2p)))
  }
  taxis <- seq(interval, T_end, interval)

  # Build a concentration series from a set of pulse locations
  build_series <- function(tau, mass_mean, mass_sd, width_mean, width_sd,
                           baseline, halflife) {
    np <- length(tau)
    A <- s2p <- mass_kappa <- width_kappa <- rep(0, max(np, 1))
    if (np > 0) {
      for (i in seq_len(np)) {
        if (pulse_distribution == "lognormal") {
          # Log-normal random effects (Gaussian on the log scale, no t-scale
          # mixture, no truncation). mass_mean/mass_sd/width_mean/width_sd are
          # LOG-scale parameters, matching the papers-default lognormal fit and
          # simulate_pulse()'s lognormal branch.
          mass_kappa[i]  <- 1
          width_kappa[i] <- 1
          A[i]   <- exp(stats::rnorm(1, mass_mean,  mass_sd))
          s2p[i] <- exp(stats::rnorm(1, width_mean, width_sd))
        } else {
          mass_kappa[i]  <- stats::rgamma(1, shape = 2, rate = 2)
          width_kappa[i] <- stats::rgamma(1, shape = 2, rate = 2)
          tvar  <- mass_sd^2  / mass_kappa[i]
          t2var <- width_sd^2 / width_kappa[i]
          while (A[i]   < 0.25) A[i]   <- stats::rnorm(1, mass_mean,  sqrt(tvar))
          while (s2p[i] < 0.5)  s2p[i] <- stats::rnorm(1, width_mean, sqrt(t2var))
        }
      }
    }
    ytmp <- 0
    if (np > 0) {
      for (i in seq_len(np)) {
        ytmp <- ytmp + meanI(taxis, 0, A[i], tau[i], log(2) / halflife, s2p[i])
      }
    }
    ysim   <- ytmp + baseline
    errors <- stats::rnorm(length(taxis), 0, sqrt(error_var))
    conc   <- ysim * exp(errors)
    list(
      data = tibble::tibble(observation = seq_along(taxis),
                            time = taxis, concentration = conc),
      parameters = tibble::tibble(
        pulse_no    = if (np > 0) seq_len(np) else integer(0),
        mass        = if (np > 0) A[seq_len(np)] else numeric(0),
        width       = if (np > 0) s2p[seq_len(np)] else numeric(0),
        location    = if (np > 0) tau else numeric(0),
        mass_kappa  = if (np > 0) mass_kappa[seq_len(np)] else numeric(0),
        width_kappa = if (np > 0) width_kappa[seq_len(np)] else numeric(0))
    )
  }

  # --- Driver hormone: standard pulsatile series ---
  driver <- simulate_pulse(num_obs = num_obs, interval = interval,
                           ipi_mean = driver_ipi_mean, ipi_var = driver_ipi_var,
                           mass_mean = driver_mass_mean, mass_sd = driver_mass_sd,
                           width_mean = driver_width_mean, width_sd = driver_width_sd,
                           constant_baseline = driver_baseline,
                           constant_halflife = driver_halflife,
                           pulse_distribution = pulse_distribution)
  tau_driver <- driver$parameters$location
  n_driver   <- length(tau_driver)

  # --- Coupling intensity: base_rate * (1 + lambda(t)) ---
  kernel_sum <- function(t) {
    if (n_driver == 0) return(rep(0, length(t)))
    vapply(t, function(ti)
      sum((rho / sqrt(2 * pi * nu)) * exp(-(ti - tau_driver)^2 / (2 * nu))),
      numeric(1))
  }
  base_rate <- response_pulse_count / T_end

  # --- Response pulse locations via thinning (Lewis-Shedler) ---
  # Upper bound: lambda(t) <= n_driver * peak_kernel for all t.
  peak_kernel <- rho / sqrt(2 * pi * nu)
  mu_max      <- base_rate * (1 + n_driver * peak_kernel)
  n_cand      <- stats::rpois(1, mu_max * T_end)
  if (n_cand > 0) {
    cand         <- sort(stats::runif(n_cand, 0, T_end))
    accept       <- stats::runif(n_cand) < (base_rate * (1 + kernel_sum(cand)) / mu_max)
    tau_response <- cand[accept]
  } else {
    tau_response <- numeric(0)
  }

  response <- build_series(tau_response, response_mass_mean, response_mass_sd,
                           response_width_mean, response_width_sd,
                           response_baseline, response_halflife)

  structure(
    list(
      driver_data         = driver$data[, c("time", "concentration")],
      response_data       = response$data[, c("time", "concentration")],
      driver_parameters   = driver$parameters,
      response_parameters = response$parameters,
      association         = list(rho = rho, nu = nu,
                                 response_pulse_count = response_pulse_count,
                                 n_driver_pulses = n_driver,
                                 n_response_pulses = length(tau_response))
    ),
    class = "joint_sim")

}


#' Simulate a multi-subject pulsatile hormone dataset
#'
#' @description \code{simulate_pulse_population} generates data for the
#'   hierarchical population model: \code{n_subjects} independent pulsatile
#'   series, each produced by \code{\link{simulate_pulse}}. Subject-level means
#'   may vary around the population means via the \code{*_sd} arguments (set them
#'   to 0, the default, for identical subjects). The result plugs directly into
#'   \code{\link{fit_pulse_population}} as its \code{data} argument.
#'
#' @param n_subjects Number of subjects to simulate.
#' @param num_obs Number of observations per subject. Window duration is
#'   \code{num_obs * interval}.
#' @param interval Time in minutes between observations.
#' @param mass_mean,width_mean Population mean pulse mass and width. If
#'   \code{NULL} (default) they resolve to natural-scale \code{3.5}/\code{35}
#'   under \code{pulse_distribution = "truncnorm"} or LOG-scale
#'   \code{1.2}/\code{3.0} under \code{"lognormal"}.
#' @param baseline,halflife Population mean baseline and half-life. When
#'   \code{halflife_sd > 0}, \code{halflife} must exceed 8 (the model's minimum
#'   half-life), otherwise the subject-level rejection sampler cannot converge.
#' @param mass_sd,width_sd Pulse-to-pulse SD of mass and width (within subject).
#'   If \code{NULL} (default) they resolve to \code{1.0}/\code{5} (truncnorm) or
#'   LOG-scale \code{0.5}/\code{0.7} (lognormal).
#' @param ipi_mean,ipi_var Inter-pulse interval mean and variance.
#' @param mass_mean_sd,width_mean_sd Subject-to-subject SD of the mean mass and
#'   mean width (0 = identical subjects).
#' @param baseline_sd,halflife_sd Subject-to-subject SD of baseline and
#'   half-life (0 = identical subjects).
#' @param pulse_distribution Character, one of \code{"truncnorm"} (default) or
#'   \code{"lognormal"}, passed to \code{\link{simulate_pulse}} for each subject
#'   to choose the per-pulse mass/width data-generating distribution. Under
#'   \code{"lognormal"} the \code{mass_sd}/\code{width_sd} (and the subject-mean
#'   \code{mass_mean}/\code{width_mean}) are LOG-scale parameters, matching the
#'   papers-default population fit.
#' @param seed Optional RNG seed for reproducibility.
#' @note This function defaults to \code{pulse_distribution = "truncnorm"},
#'   whereas \code{\link{population_spec}} defaults to \code{"lognormal"}.
#'   Simulating with defaults and fitting with defaults therefore fits a
#'   lognormal model to truncnorm-generated data. Pass \code{pulse_distribution =
#'   "lognormal"} to match the papers-default fit.
#' @return An object of class \code{population_sim}: a list with \code{data} (a
#'   list of per-subject data frames with \code{time} and \code{concentration},
#'   ready for \code{fit_pulse_population}), \code{n_subjects}, and \code{truth}
#'   (the population means and all SD arguments used, for comparison against
#'   recovered posteriors).
#' @seealso \code{\link{simulate_pulse}}, \code{\link{fit_pulse_population}}
#' @keywords pulse simulation
#' @examples
#' sim <- simulate_pulse_population(n_subjects = 4, seed = 1)
#' length(sim$data)
#' @importFrom stats rnorm
#' @export
simulate_pulse_population <- function(n_subjects   = 5,
                                      num_obs       = 144,
                                      interval      = 10,
                                      mass_mean     = NULL,
                                      width_mean    = NULL,
                                      baseline      = 2.6,
                                      halflife      = 45,
                                      mass_sd       = NULL,
                                      width_sd      = NULL,
                                      ipi_mean      = 12,
                                      ipi_var       = 40,
                                      mass_mean_sd  = 0,
                                      width_mean_sd = 0,
                                      baseline_sd   = 0,
                                      halflife_sd   = 0,
                                      pulse_distribution = c("truncnorm",
                                                             "lognormal"),
                                      seed          = NULL) {

  pulse_distribution <- match.arg(pulse_distribution)

  # Conditional defaults for the pulse mass/width mean and SD, forwarded to
  # simulate_pulse() for each subject. Natural-scale under "truncnorm" (legacy),
  # LOG-scale under "lognormal" -- otherwise exp() of a natural-scale width
  # produces absurd values. A user-supplied value always wins.
  if (pulse_distribution == "lognormal") {
    if (is.null(mass_mean))  mass_mean  <- 1.2
    if (is.null(mass_sd))    mass_sd    <- 0.5
    if (is.null(width_mean)) width_mean <- 3.0
    if (is.null(width_sd))   width_sd   <- 0.7
  } else {
    if (is.null(mass_mean))  mass_mean  <- 3.5
    if (is.null(mass_sd))    mass_sd    <- 1.0
    if (is.null(width_mean)) width_mean <- 35
    if (is.null(width_sd))   width_sd   <- 5
  }

  if (!is.null(seed)) set.seed(seed)
  stopifnot(n_subjects >= 1)
  # The deconvolution model's minimum half-life is 8 minutes. When half-life
  # varies across subjects, the subject-level rejection sampler cannot draw a
  # value above that floor unless the population mean exceeds it, so require it.
  if (halflife_sd > 0) stopifnot(halflife > 8)

  # Draw a positive value from N(mean, sd) (sd = 0 returns the mean unchanged)
  rpos <- function(mean, sd, floor = 1e-6) {
    if (sd <= 0) return(mean)
    v <- -1
    while (v <= floor) v <- stats::rnorm(1, mean, sd)
    v
  }

  subjects <- lapply(seq_len(n_subjects), function(i) {
    subj_mass     <- rpos(mass_mean,  mass_mean_sd)
    subj_width    <- rpos(width_mean, width_mean_sd)
    subj_baseline <- rpos(baseline,   baseline_sd)
    subj_halflife <- rpos(halflife,   halflife_sd, floor = 8)
    s <- simulate_pulse(num_obs = num_obs, interval = interval,
                        ipi_mean = ipi_mean, ipi_var = ipi_var,
                        mass_mean = subj_mass, mass_sd = mass_sd,
                        width_mean = subj_width, width_sd = width_sd,
                        constant_baseline = subj_baseline,
                        constant_halflife = subj_halflife,
                        pulse_distribution = pulse_distribution)
    data.frame(time = s$data$time, concentration = s$data$concentration)
  })

  structure(
    list(data = subjects,
         n_subjects = n_subjects,
         truth = list(mass_mean = mass_mean, width_mean = width_mean,
                      baseline = baseline, halflife = halflife,
                      mass_mean_sd = mass_mean_sd, width_mean_sd = width_mean_sd,
                      baseline_sd = baseline_sd, halflife_sd = halflife_sd,
                      mass_sd = mass_sd, width_sd = width_sd)),
    class = "population_sim")

}


#-------------------------------------------------------------------------------
# End of file
#-------------------------------------------------------------------------------
