#-------------------------------------------------------------------------------
# fit_population.R - Population model fitting interface
#-------------------------------------------------------------------------------

#' fit_pulse_population
#'
#' Fits a hierarchical Bayesian deconvolution model for pulsatile hormone data
#' from multiple subjects. This function estimates population-level parameters
#' (means and variances) as well as subject-specific and pulse-level parameters.
#'
#' @param data A list of data frames, where each element corresponds to one
#'   subject and contains columns for \code{time} and \code{concentration}.
#'   Alternatively, a single data frame with a \code{subject_id} column.
#' @param time A string. Name of the time variable (default: "time").
#' @param conc A string. Name of the concentration variable (default: "concentration").
#' @param subject_id A string. Name of the subject ID variable if \code{data}
#'   is a single data frame (default: NULL). If provided, data will be split
#'   by this variable.
#' @param spec An object of class \code{population_spec}, created by
#'   \code{population_spec()}, specifying priors, starting values, and proposal
#'   variances for the population model.
#' @param subject_starts Optional list of subject-specific starting values. If
#'   NULL, starting values will be initialized from population starting values
#'   with small random perturbations. Each element should be a list with:
#'   \code{baseline}, \code{halflife}, \code{mass_mean}, \code{width_mean}.
#' @param iters Number of MCMC iterations (default: 250000).
#' @param thin Thinning interval - save every \code{thin}th sample (default: 50).
#' @param burnin Number of initial samples to discard as burn-in (default: 0.1 * iters).
#' @param use_tibble Return chains as tibble data frames for prettier printing
#'   (default: TRUE). Requires tibble package.
#' @param verbose Print diagnostic output every 5000 iterations including
#'   parameter estimates and acceptance rates (default: FALSE).
#'
#' @return An object of class \code{population_fit}, a list containing the
#'   population-level MCMC chain (\code{population_chain}), subject-level chains
#'   (\code{subject_chains}), pulse-level chains (\code{pulse_chains}), the
#'   original data, number of subjects, fitting options, and model specification.
#'
#' @details
#' The population model estimates a three-level hierarchy:
#'
#' \strong{Population Level:}
#' \itemize{
#'   \item Mean pulse mass and width across subjects
#'   \item Mean baseline and halflife across subjects
#'   \item Subject-to-subject SD of mean mass and width
#'   \item Subject-to-subject SD of baseline and halflife
#'   \item Pulse-to-pulse SD of mass and width (shared across subjects)
#'   \item Model error variance
#' }
#'
#' \strong{Subject Level (for each subject):}
#' \itemize{
#'   \item Subject-specific mean pulse mass and width
#'   \item Subject-specific baseline and halflife
#'   \item Number of pulses (varies by MCMC iteration)
#' }
#'
#' \strong{Pulse Level (for each pulse within each subject):}
#' \itemize{
#'   \item Individual pulse mass and width
#'   \item Pulse location (time)
#'   \item T-distribution scale parameters for pulse variability
#' }
#'
#' @importFrom tibble as_tibble
#' @export
#' @keywords pulse population
#' @examples
#' \dontrun{
#' # Simulate data for 5 subjects
#' set.seed(123)
#' sim_data <- lapply(1:5, function(i) {
#'   simulate_pulse(n_obs = 24, true_pulses = 12)$data
#' })
#'
#' # Create specification
#' pop_spec <- population_spec(
#'   prior_mass_mean_mean = 3.5,
#'   prior_width_mean_mean = 42,
#'   prior_mean_pulse_count = 12
#' )
#'
#' # Fit model
#' fit <- fit_pulse_population(
#'   data = sim_data,
#'   spec = pop_spec,
#'   iters = 10000,
#'   thin = 10,
#'   burnin = 5000,
#'   verbose = TRUE
#' )
#'
#' # Examine results
#' summary(fit$population_chain)
#' plot(fit$population_chain$mass_mean, type = "l")
#' }
fit_pulse_population <- function(data,
                                  time = "time",
                                  conc = "concentration",
                                  subject_id = NULL,
                                  spec,
                                  subject_starts = NULL,
                                  iters = 250000,
                                  thin = 50,
                                  burnin = as.integer(0.1 * iters),
                                  use_tibble = TRUE,
                                  verbose = FALSE) {
  
  # Validate inputs
  if (!inherits(spec, "population_spec")) {
    stop("spec must be a population_spec object created by population_spec()")
  }
  
  stopifnot(is.logical(use_tibble), is.logical(verbose))
  
  if (burnin >= iters) {
    stop("burnin must be < iters")
  }
  if (thin < 1) {
    stop("thin must be a positive integer")
  }

  # Process data into list format
  if (!is.null(subject_id)) {
    # Split single data frame by subject_id
    if (!is.data.frame(data)) {
      stop("If subject_id is provided, data must be a data frame")
    }
    if (!(subject_id %in% names(data))) {
      stop(paste("Column", subject_id, "not found in data"))
    }
    
    subject_ids <- unique(data[[subject_id]])
    subject_data_list <- lapply(subject_ids, function(sid) {
      subj_data <- data[data[[subject_id]] == sid, , drop = FALSE]
      list(time = subj_data[[time]],
           concentration = subj_data[[conc]])
    })
    names(subject_data_list) <- as.character(subject_ids)
    
  } else {
    # Data should be a list of data frames
    if (!is.list(data)) {
      stop("data must be either a list of data frames or a single data frame with subject_id")
    }
    
    # Convert each data frame to required format
    subject_data_list <- lapply(data, function(subj_data) {
      if (inherits(subj_data, "pulse_sim")) {
        subj_data <- subj_data$data
      }
      list(time = subj_data[[time]],
           concentration = subj_data[[conc]])
    })
  }
  
  n_subjects <- length(subject_data_list)
  
  if (n_subjects == 0) {
    stop("No subjects found in data")
  }
  
  # Validate each subject has data
  for (i in seq_along(subject_data_list)) {
    if (length(subject_data_list[[i]]$time) == 0) {
      stop(paste("Subject", i, "has no observations"))
    }
    if (length(subject_data_list[[i]]$time) != 
        length(subject_data_list[[i]]$concentration)) {
      stop(paste("Subject", i, ": time and concentration must have same length"))
    }
  }
  
  # Create subject starting values
  if (is.null(subject_starts)) {
    # Initialize from population starting values with small random variation.
    # NOTE: Do NOT re-seed the RNG here. A previous `set.seed(NULL)` call
    # ("ensure randomness") re-seeded from OS entropy on every fit, which
    # overrode the user's set.seed() and made fits non-reproducible. The jitter
    # below now draws from the user's seeded RNG stream so set.seed() controls
    # the whole fit (mirroring the single-subject model).
    subject_starts <- lapply(1:n_subjects, function(i) {
      list(
        baseline = spec$population_starting_values$baseline_mean + 
          rnorm(1, 0, spec$population_starting_values$baseline_sd / 2),
        halflife = spec$population_starting_values$halflife_mean + 
          rnorm(1, 0, spec$population_starting_values$halflife_sd / 2),
        mass_mean = spec$population_starting_values$mass_mean + 
          rnorm(1, 0, spec$population_starting_values$mass_mean_sd / 2),
        width_mean = spec$population_starting_values$width_mean + 
          rnorm(1, 0, spec$population_starting_values$width_mean_sd / 2)
      )
    })
  } else {
    # Validate provided starting values
    if (length(subject_starts) != n_subjects) {
      stop("subject_starts must have one element per subject")
    }
    for (i in seq_along(subject_starts)) {
      required_fields <- c("baseline", "halflife", "mass_mean", "width_mean")
      if (!all(required_fields %in% names(subject_starts[[i]]))) {
        stop(paste("Subject", i, "starting values missing required fields:",
                   paste(required_fields, collapse = ", ")))
      }
    }
  }
  
  # Proposal variance adaptation parameters
  pv_adjust_iter <- 500
  pv_adjust_max_iter <- 25000
  univariate_pv_target_ratio <- 0.35
  bivariate_pv_target_ratio <- 0.25
  
  # Prepare inputs for C++ function
  population_priors <- structure(
    spec$population_priors,
    class = "bp_population_priors"
  )
  
  proposal_variances <- structure(
    spec$proposal_variances,
    class = "bp_proposalvariance"
  )
  
  population_startingvals <- structure(
    spec$population_starting_values,
    class = "bp_population_startingvals"
  )
  
  # Store call for output
  Call <- match.call()
  
  # Call C++ function
  fit <- population_(
    subject_data_list = subject_data_list,
    location_prior = spec$location_prior,
    population_priors = population_priors,
    proposalvars = proposal_variances,
    population_startingvals = population_startingvals,
    subject_startingvals_list = subject_starts,
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
  population_chain <- as.data.frame(fit$population_chain)
  
  # Convert subject chains
  subject_chains <- lapply(fit$subject_chains, as.data.frame)
  
  # Convert pulse chains
  pulse_chains <- lapply(fit$pulse_chains, function(subj_pulses) {
    lapply(subj_pulses, as.data.frame)
  })
  
  # Convert to tibbles if requested
  if (use_tibble) {
    population_chain <- tibble::as_tibble(population_chain)
    subject_chains <- lapply(subject_chains, tibble::as_tibble)
    pulse_chains <- lapply(pulse_chains, function(subj_pulses) {
      lapply(subj_pulses, tibble::as_tibble)
    })
  }
  
  # Create return object
  rtn_obj <- structure(
    list(
      model = "population",
      call = Call,
      population_chain = population_chain,
      subject_chains = subject_chains,
      pulse_chains = pulse_chains,
      data = data,
      num_subjects = n_subjects,
      options = list(
        time = time,
        conc = conc,
        thinning = thin,
        iterations = iters,
        burnin = burnin
      ),
      spec = spec
    ),
    class = "population_fit"
  )
  
  return(rtn_obj)
}


#' @importFrom utils head
#' @export
print.population_fit <- function(x, ...) {
  cat("\nPopulation Pulse Model Fit\n\n")
  cat("Model type:", x$model, "\n")
  cat("Number of subjects:", x$num_subjects, "\n")
  cat("MCMC iterations:", x$options$iterations, "\n")
  cat("Thinning:", x$options$thinning, "\n")
  cat("Burnin:", x$options$burnin, "\n")
  cat("Saved samples:", nrow(x$population_chain), "\n\n")
  
  cat("Population-level parameters:\n")
  print(head(x$population_chain))
  cat("\n")
  
  invisible(x)
}


#' Summary method for population_fit objects
#' 
#' @param object A population_fit object
#' @param diagnostics Logical, whether to include convergence diagnostics
#'   (default: TRUE)
#' @param ... Additional arguments (not used)
#' @export
summary.population_fit <- function(object, diagnostics = TRUE, ...) {
  cat("\n=== Population Model Summary ===\n\n")
  cat("Number of subjects:", object$num_subjects, "\n")
  cat("MCMC samples:", nrow(object$population_chain), "\n")
  cat("Iterations:", object$options$iterations, "\n")
  cat("Burnin:", object$options$burnin, "\n")
  cat("Thinning:", object$options$thinning, "\n\n")
  
  cat("Population-level parameter estimates (posterior means):\n")
  pop_chain <- object$population_chain
  if ("iteration" %in% names(pop_chain)) {
    pop_chain <- pop_chain[, names(pop_chain) != "iteration", drop = FALSE]
  }
  pop_means <- colMeans(pop_chain)
  print(pop_means)
  cat("\n")
  
  # Convergence diagnostics
  if (diagnostics) {
    cat("Convergence Diagnostics:\n")
    cat(rep("-", 70), "\n", sep = "")
    
    diag <- population_diagnostics(object)
    
    # Summary statistics
    n_converged <- sum(diag$converged, na.rm = TRUE)
    n_total <- nrow(diag)
    min_ess <- min(diag$ess, na.rm = TRUE)
    max_ess <- max(diag$ess, na.rm = TRUE)
    
    cat(sprintf("Parameters with ESS > 100: %d / %d (%.1f%%)\n",
                n_converged, n_total, 100 * n_converged / n_total))
    cat(sprintf("ESS range: %.0f - %.0f\n", min_ess, max_ess))
    
    # Flag parameters with low ESS
    low_ess <- diag[diag$ess < 100, ]
    if (nrow(low_ess) > 0) {
      cat("\n[WARNING] Parameters with low ESS:\n")
      for (i in seq_len(nrow(low_ess))) {
        cat(sprintf("  %s: ESS = %.0f\n", low_ess$parameter[i], low_ess$ess[i]))
      }
      cat("\nConsider running longer or adjusting proposal variances.\n")
    } else {
      cat("\n[OK] All parameters have adequate ESS\n")
    }
    cat("\n")
  }
  
  cat("Subject-level parameter estimates:\n")
  cat("(First 5 subjects)\n")
  for (i in seq_len(min(5, object$num_subjects))) {
    cat(sprintf("\nSubject %d:\n", i))
    subj_means <- colMeans(object$subject_chains[[i]])
    print(subj_means)
  }
  
  if (diagnostics) {
    cat("\nFor detailed diagnostics, use: convergence_report(fit)\n")
    cat("For trace plots, use: plot_trace(fit)\n")
    cat("For ACF plots, use: plot_acf(fit)\n")
  }
  
  invisible(object)
}


#------------------------------------------------------------------------------#
#    End of file                                                               #
#------------------------------------------------------------------------------#
