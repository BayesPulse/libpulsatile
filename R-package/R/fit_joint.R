#-------------------------------------------------------------------------------
# fit_joint.R - Wrapper for joint hormone model estimation
#-------------------------------------------------------------------------------

#' fit_pulse_joint
#'
#' Fits a joint Bayesian deconvolution model for two coupled pulsatile hormones
#' measured from a single subject. The driver hormone pulses influence the
#' intensity of response hormone pulse occurrence through a coupling kernel
#' parameterized by rho (coupling strength) and nu (temporal spread).
#'
#' @param driver_data A data frame containing driver hormone measurements with
#'   columns for time and concentration. Can also be a \code{pulse_sim} object.
#' @param response_data A data frame containing response hormone measurements
#'   with columns for time and concentration. Can also be a \code{pulse_sim} object.
#' @param driver_time A string. Name of the time variable in \code{driver_data}
#'   (default: "time").
#' @param driver_conc A string. Name of the concentration variable in
#'   \code{driver_data} (default: "concentration").
#' @param response_time A string. Name of the time variable in \code{response_data}
#'   (default: "time").
#' @param response_conc A string. Name of the concentration variable in
#'   \code{response_data} (default: "concentration").
#' @param spec An object of class \code{joint_spec}, created by
#'   \code{joint_spec()}, specifying the priors, starting values, and proposal
#'   variances to use for both hormones and their coupling.
#' @param iters Number of MCMC iterations (default: 250000).
#' @param thin Thinning interval - save every \code{thin}th sample (default: 50).
#' @param burnin Number of initial samples to discard as burn-in (default: 0.1 * iters).
#' @param use_tibble Return chains as tibble data frames for prettier printing
#'   (default: TRUE). Requires tibble package.
#' @param verbose Print diagnostic output every 5000 iterations including
#'   parameter estimates and acceptance rates (default: FALSE).
#'
#' @return An object of class \code{joint_fit}, a list containing:
#'   \item{model}{Model type ("joint")}
#'   \item{call}{The matched function call}
#'   \item{driver_chain}{MCMC chain for driver hormone parameters}
#'   \item{response_chain}{MCMC chain for response hormone parameters}
#'   \item{association_chain}{MCMC chain for coupling parameters (rho, nu)}
#'   \item{driver_pulse_chain}{MCMC chain for individual driver pulses}
#'   \item{response_pulse_chain}{MCMC chain for individual response pulses}
#'   \item{driver_data}{Original driver data}
#'   \item{response_data}{Original response data}
#'   \item{time_range}{Time range of observations}
#'   \item{options}{List of fitting options (time/conc column names, thinning, iterations)}
#'   \item{spec}{Model specification used for fitting}
#'
#' @details
#' The joint model estimates:
#'
#' \strong{Driver Hormone:}
#' \itemize{
#'   \item Mean pulse mass and width
#'   \item Baseline concentration and half-life
#'   \item Pulse-to-pulse SD of mass and width
#'   \item Number of pulses (varies by MCMC iteration)
#'   \item Individual pulse locations, masses, and widths
#' }
#'
#' \strong{Response Hormone:}
#' \itemize{
#'   \item Same parameters as driver
#'   \item Pulse occurrence influenced by driver pulses
#' }
#'
#' \strong{Coupling Parameters:}
#' \itemize{
#'   \item rho: Coupling strength (how much driver affects response)
#'   \item nu: Temporal spread of coupling effect (variance)
#' }
#'
#' The coupling is modeled through an intensity function:
#' λ(t) = Σ_j [ρ / √(2πν)] × exp(-(t - τ_j)² / (2ν))
#' where τ_j are driver pulse locations.
#'
#' @importFrom tibble as_tibble
#' @export
#' @keywords pulse joint
#' @examples
#' \dontrun{
#' # Simulate driver and response data
#' set.seed(123)
#' driver <- simulate_pulse(n_obs = 24, true_pulses = 12)
#' response <- simulate_pulse(n_obs = 24, true_pulses = 15)
#'
#' # Create specification
#' joint_spec <- joint_spec(
#'   prior_log_rho_mean = 0.0,
#'   prior_log_nu_mean = 3.0
#' )
#'
#' # Fit joint model
#' fit <- fit_pulse_joint(
#'   driver_data = driver,
#'   response_data = response,
#'   spec = joint_spec,
#'   iters = 10000,
#'   thin = 10,
#'   burnin = 5000,
#'   verbose = TRUE
#' )
#'
#' # Examine coupling parameters
#' summary(fit$association_chain)
#' plot(fit$association_chain$rho, type = "l")
#' }
fit_pulse_joint <- function(driver_data,
                            response_data,
                            driver_time = "time",
                            driver_conc = "concentration",
                            response_time = "time",
                            response_conc = "concentration",
                            spec,
                            iters = 250000,
                            thin = 50,
                            burnin = as.integer(0.1 * iters),
                            use_tibble = TRUE,
                            verbose = FALSE) {

  # Validate spec
  if (!inherits(spec, "joint_spec")) {
    stop("spec must be a joint_spec object created by joint_spec()")
  }

  stopifnot(is.logical(use_tibble), is.logical(verbose))

  if (burnin >= iters) {
    stop("burnin must be < iters")
  }

  # Process driver data
  if (inherits(driver_data, "pulse_sim")) {
    driver_df <- driver_data$data
  } else {
    driver_df <- data.frame(
      time = driver_data[[driver_time]],
      concentration = driver_data[[driver_conc]]
    )
  }

  # Process response data
  if (inherits(response_data, "pulse_sim")) {
    response_df <- response_data$data
  } else {
    response_df <- data.frame(
      time = response_data[[response_time]],
      concentration = response_data[[response_conc]]
    )
  }

  # Validate input data
  stopifnot(
    is.numeric(driver_df$time),
    is.numeric(driver_df$concentration),
    is.numeric(response_df$time),
    is.numeric(response_df$concentration)
  )

  if (nrow(driver_df) == 0 || nrow(response_df) == 0) {
    stop("Driver and response data must have at least one observation")
  }

  # Proposal variance adaptation parameters
  pv_adjust_iter <- 500
  pv_adjust_max_iter <- 25000
  univariate_pv_target_ratio <- 0.35
  bivariate_pv_target_ratio <- 0.25

  # Prepare inputs for C++ function (same pattern as single subject)
  # Process lists with lapply to ensure proper structure
  # Use as.numeric to ensure scalars, not vectors
  driver_priors <- lapply(spec$driver_priors, function(x) as.numeric(ifelse(is.null(x), NA, x)))
  response_priors <- lapply(spec$response_priors, function(x) as.numeric(ifelse(is.null(x), NA, x)))
  proposal_variances <- lapply(spec$proposal_variances, function(x) as.numeric(ifelse(is.null(x), NA, x)))
  driver_startingvals <- lapply(spec$driver_starting_values, function(x) as.numeric(ifelse(is.null(x), NA, x)))
  response_startingvals <- lapply(spec$response_starting_values, function(x) as.numeric(ifelse(is.null(x), NA, x)))

  # Create association lists from scratch (avoid any processing that might break names)
  association_priors <- list(
    log_rho_mean = spec$association_priors$log_rho_mean,
    log_rho_var = spec$association_priors$log_rho_var,
    log_nu_mean = spec$association_priors$log_nu_mean,
    log_nu_var = spec$association_priors$log_nu_var
  )

  association_startingvals <- list(
    rho = spec$association_starting_values$rho,
    nu = spec$association_starting_values$nu
  )

  # Add class attributes
  driver_priors <- structure(driver_priors, class = "bp_priors")
  response_priors <- structure(response_priors, class = "bp_priors")
  # Don't add class to association lists - may interfere with C++ access
  # association_priors <- structure(association_priors, class = "bp_association_priors")
  proposal_variances <- structure(proposal_variances, class = "bp_proposalvariance")
  driver_startingvals <- structure(driver_startingvals, class = "bp_startingvals")
  response_startingvals <- structure(response_startingvals, class = "bp_startingvals")
  # association_startingvals <- structure(association_startingvals, class = "bp_association_startingvals")

  # Store call for output
  Call <- match.call()

  # Debug output removed - was not causing the issue

  # Call C++ function
  fit <- jointsinglesubject_(
    driver_concentration = driver_df$concentration,
    driver_time = driver_df$time,
    response_concentration = response_df$concentration,
    response_time = response_df$time,
    location_prior = spec$location_prior,
    driver_priors = driver_priors,
    response_priors = response_priors,
    association_priors = association_priors,
    proposalvars = proposal_variances,
    driver_startingvals = driver_startingvals,
    response_startingvals = response_startingvals,
    association_startingvals = association_startingvals,
    mcmc_iterations = as.integer(iters),
    thin = as.integer(thin),
    burnin = as.integer(burnin),
    verbose = verbose,
    pv_adjust_iter = as.integer(pv_adjust_iter),
    pv_adjust_max_iter = as.integer(pv_adjust_max_iter),
    bivariate_pv_target_ratio = bivariate_pv_target_ratio,
    univariate_pv_target_ratio = univariate_pv_target_ratio
  )

  # Process output
  # Add column names and convert to data frames
  driver_chain <- as.data.frame(fit$driver_fixed_effects)
  colnames(driver_chain) <- fit$driver_colnames
  driver_chain$num_pulses <- sapply(fit$driver_pulses, nrow)

  response_chain <- as.data.frame(fit$response_fixed_effects)
  colnames(response_chain) <- fit$response_colnames
  response_chain$num_pulses <- sapply(fit$response_pulses, nrow)

  association_chain <- as.data.frame(fit$association)
  colnames(association_chain) <- fit$association_colnames

  # Process pulse chains
  driver_pulse_chain <- as.data.frame(do.call(rbind, fit$driver_pulses))
  colnames(driver_pulse_chain) <- fit$driver_pulse_colnames

  response_pulse_chain <- as.data.frame(do.call(rbind, fit$response_pulses))
  colnames(response_pulse_chain) <- fit$response_pulse_colnames

  # Convert to tibbles if requested
  if (use_tibble) {
    driver_chain <- tibble::as_tibble(driver_chain)
    response_chain <- tibble::as_tibble(response_chain)
    association_chain <- tibble::as_tibble(association_chain)
    driver_pulse_chain <- tibble::as_tibble(driver_pulse_chain)
    response_pulse_chain <- tibble::as_tibble(response_pulse_chain)
  }

  # Create return object
  rtn_obj <- structure(
    list(
      model = "joint",
      call = Call,
      driver_chain = driver_chain,
      response_chain = response_chain,
      association_chain = association_chain,
      driver_pulse_chain = driver_pulse_chain,
      response_pulse_chain = response_pulse_chain,
      driver_data = driver_data,
      response_data = response_data,
      time_range = c(
        min(driver_df$time, response_df$time),
        max(driver_df$time, response_df$time)
      ),
      options = list(
        driver_time = driver_time,
        driver_conc = driver_conc,
        response_time = response_time,
        response_conc = response_conc,
        thinning = thin,
        iterations = iters,
        burnin = burnin
      ),
      spec = spec
    ),
    class = "joint_fit"
  )

  return(rtn_obj)
}


#' @importFrom utils head
#' @export
print.joint_fit <- function(x, ...) {
  cat("\nJoint Hormone Model Fit\n\n")
  cat("Model type:", x$model, "\n")
  cat("MCMC iterations:", x$options$iterations, "\n")
  cat("Thinning:", x$options$thinning, "\n")
  cat("Burnin:", x$options$burnin, "\n")
  cat("Saved samples:", nrow(x$driver_chain), "\n\n")

  cat("Driver Hormone Parameters:\n")
  print(head(x$driver_chain))
  cat("\n")

  cat("Response Hormone Parameters:\n")
  print(head(x$response_chain))
  cat("\n")

  cat("Association Parameters (Driver -> Response):\n")
  print(head(x$association_chain))
  cat("\n")

  invisible(x)
}


#' Summary method for joint_fit objects
#'
#' @param object A joint_fit object
#' @param ... Additional arguments (not used)
#' @export
summary.joint_fit <- function(object, ...) {
  cat("\n=== Joint Hormone Model Summary ===\n\n")
  cat("MCMC samples:", nrow(object$driver_chain), "\n")
  cat("Iterations:", object$options$iterations, "\n")
  cat("Burnin:", object$options$burnin, "\n")
  cat("Thinning:", object$options$thinning, "\n\n")

  cat("Driver Hormone Parameter Estimates (posterior means):\n")
  driver_means <- colMeans(object$driver_chain[, !names(object$driver_chain) %in% c("iteration"), drop = FALSE])
  print(driver_means)
  cat("\n")

  cat("Response Hormone Parameter Estimates (posterior means):\n")
  response_means <- colMeans(object$response_chain[, !names(object$response_chain) %in% c("iteration"), drop = FALSE])
  print(response_means)
  cat("\n")

  cat("Association Parameter Estimates (posterior means):\n")
  assoc_means <- colMeans(object$association_chain[, !names(object$association_chain) %in% c("iteration"), drop = FALSE])
  print(assoc_means)
  cat("\n")

  cat("Coupling Strength (rho):\n")
  cat("  Mean:", mean(object$association_chain$rho), "\n")
  cat("  SD:", sd(object$association_chain$rho), "\n")
  cat("  95% CI: [",
      quantile(object$association_chain$rho, 0.025), ",",
      quantile(object$association_chain$rho, 0.975), "]\n\n")

  cat("Coupling Temporal Spread (nu):\n")
  cat("  Mean:", mean(object$association_chain$nu), "\n")
  cat("  SD:", sd(object$association_chain$nu), "\n")
  cat("  95% CI: [",
      quantile(object$association_chain$nu, 0.025), ",",
      quantile(object$association_chain$nu, 0.975), "]\n\n")

  cat("Driver Pulse Count (posterior mean):",
      mean(object$driver_chain$num_pulses), "\n")
  cat("Response Pulse Count (posterior mean):",
      mean(object$response_chain$num_pulses), "\n\n")

  invisible(object)
}


#------------------------------------------------------------------------------#
#    End of file                                                               #
#------------------------------------------------------------------------------#
