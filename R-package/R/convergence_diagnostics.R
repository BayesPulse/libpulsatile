#-------------------------------------------------------------------------------
# convergence_diagnostics.R - MCMC convergence diagnostic functions
#-------------------------------------------------------------------------------

# Declare global variables to avoid R CMD check NOTEs for NSE in ggplot2
utils::globalVariables(c("iteration", "value", "lag", "acf"))

#' Calculate Effective Sample Size (ESS)
#'
#' Computes the effective sample size for an MCMC chain, accounting for
#' autocorrelation. Higher ESS indicates better mixing and more independent
#' samples.
#'
#' @param x Numeric vector of MCMC samples
#' @param max_lag Maximum lag for autocorrelation calculation (default: length(x)/2)
#'
#' @return Effective sample size (numeric)
#'
#' @details
#' ESS is calculated as: n / (1 + 2 * sum(autocorrelation))
#' where the sum is over positive autocorrelations up to max_lag.
#'
#' Rules of thumb:
#' \itemize{
#'   \item ESS > 100: Generally acceptable
#'   \item ESS > 400: Good
#'   \item ESS > 1000: Excellent
#' }
#'
#' @importFrom stats acf na.pass sd var
#' @export
effective_sample_size <- function(x, max_lag = NULL) {
  
  n <- length(x)
  
  if (n < 10) {
    warning("Sample size too small for reliable ESS calculation")
    return(NA_real_)
  }
  
  # Remove NA values
  x <- x[!is.na(x)]
  n <- length(x)
  
  if (is.null(max_lag)) {
    max_lag <- min(n - 1, floor(n / 2))
  }
  
  # Calculate autocorrelation
  acf_vals <- stats::acf(x, lag.max = max_lag, plot = FALSE, na.action = na.pass)$acf
  
  # Sum positive autocorrelations (skip lag 0 which is always 1)
  # Stop when autocorrelation becomes negative or very small
  rho_sum <- 0
  for (k in 2:length(acf_vals)) {
    if (is.na(acf_vals[k]) || acf_vals[k] < 0.05) break
    rho_sum <- rho_sum + acf_vals[k]
  }
  
  # ESS formula
  ess <- n / (1 + 2 * rho_sum)
  
  return(ess)
}


#' Calculate Gelman-Rubin Diagnostic (R-hat)
#'
#' Computes the Gelman-Rubin potential scale reduction factor for assessing
#' convergence across multiple MCMC chains. Requires at least 2 chains.
#'
#' @param chains List of numeric vectors, where each element is an MCMC chain
#'   for the same parameter
#'
#' @return R-hat statistic (numeric)
#'
#' @details
#' R-hat compares between-chain and within-chain variance. Values close to 1
#' indicate convergence.
#'
#' Rules of thumb:
#' \itemize{
#'   \item R-hat < 1.01: Excellent convergence
#'   \item R-hat < 1.05: Good convergence
#'   \item R-hat < 1.1: Acceptable convergence
#'   \item R-hat > 1.1: Poor convergence, run longer
#' }
#'
#' @references
#' Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
#' using multiple sequences. Statistical Science, 7(4), 457-472.
#'
#' @export
gelman_rubin <- function(chains) {
  
  if (!is.list(chains)) {
    stop("chains must be a list of numeric vectors")
  }
  
  if (length(chains) < 2) {
    stop("Need at least 2 chains for Gelman-Rubin diagnostic")
  }
  
  # Number of chains and samples per chain
  m <- length(chains)
  n <- sapply(chains, length)
  
  # Check all chains have same length
  if (length(unique(n)) > 1) {
    warning("Chains have different lengths, using minimum length")
    n_min <- min(n)
    chains <- lapply(chains, function(x) x[1:n_min])
    n <- n_min
  } else {
    n <- n[1]
  }
  
  # Calculate chain means
  chain_means <- sapply(chains, mean)
  
  # Overall mean
  grand_mean <- mean(chain_means)
  
  # Between-chain variance
  B <- n * sum((chain_means - grand_mean)^2) / (m - 1)
  
  # Within-chain variance
  chain_vars <- sapply(chains, stats::var)
  W <- mean(chain_vars)
  
  # Variance estimate
  var_hat <- ((n - 1) / n) * W + (1 / n) * B
  
  # R-hat
  R_hat <- sqrt(var_hat / W)
  
  return(R_hat)
}


#' Compute Convergence Diagnostics for Population Fit
#'
#' Calculates ESS and autocorrelation for all population-level parameters
#' in a population_fit object.
#'
#' @param fit A \code{population_fit} object from \code{fit_pulse_population()}
#' @param parameters Optional character vector of parameter names to include.
#'   If NULL, computes for all parameters except 'iteration'.
#'
#' @return A data frame with columns:
#'   \item{parameter}{Parameter name}
#'   \item{mean}{Posterior mean}
#'   \item{sd}{Posterior standard deviation}
#'   \item{ess}{Effective sample size}
#'   \item{ess_per_sample}{ESS divided by total samples}
#'   \item{lag1_acf}{Lag-1 autocorrelation}
#'   \item{converged}{Logical, TRUE if ESS > 100}
#'
#' @export
population_diagnostics <- function(fit, parameters = NULL) {
  
  if (!inherits(fit, "population_fit")) {
    stop("fit must be a population_fit object")
  }
  
  pop_chain <- fit$population_chain
  
  # Remove iteration column if present
  if ("iteration" %in% names(pop_chain)) {
    pop_chain <- pop_chain[, names(pop_chain) != "iteration", drop = FALSE]
  }
  
  # Select parameters
  if (is.null(parameters)) {
    parameters <- names(pop_chain)
  } else {
    parameters <- intersect(parameters, names(pop_chain))
    if (length(parameters) == 0) {
      stop("No valid parameters found")
    }
  }
  
  # Calculate diagnostics for each parameter
  results <- lapply(parameters, function(param) {
    x <- pop_chain[[param]]
    
    ess <- effective_sample_size(x)
    acf_vals <- stats::acf(x, lag.max = 10, plot = FALSE)$acf
    
    data.frame(
      parameter = param,
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      ess = ess,
      ess_per_sample = ess / length(x),
      lag1_acf = acf_vals[2],  # Lag 1 autocorrelation
      converged = ess > 100,
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  
  return(results_df)
}


#' Plot MCMC Trace Plots
#'
#' Creates trace plots for population-level parameters to visually assess
#' convergence and mixing.
#'
#' @param fit A \code{population_fit} object
#' @param parameters Character vector of parameter names to plot. If NULL,
#'   plots all population parameters.
#' @param ncol Number of columns in plot grid (default: 2)
#'
#' @return A ggplot2 object
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs theme_minimal
#' @importFrom tidyr pivot_longer
#' @export
plot_trace <- function(fit, parameters = NULL, ncol = 2) {
  
  if (!inherits(fit, "population_fit")) {
    stop("fit must be a population_fit object")
  }
  
  pop_chain <- fit$population_chain
  
  # Remove iteration column for selection, but keep for plotting
  param_names <- setdiff(names(pop_chain), "iteration")
  
  if (is.null(parameters)) {
    parameters <- param_names
  } else {
    parameters <- intersect(parameters, param_names)
    if (length(parameters) == 0) {
      stop("No valid parameters found")
    }
  }
  
  # Prepare data for plotting
  if (!"iteration" %in% names(pop_chain)) {
    pop_chain$iteration <- seq_len(nrow(pop_chain))
  }
  
  plot_data <- pop_chain[, c("iteration", parameters), drop = FALSE]
  plot_data_long <- tidyr::pivot_longer(
    plot_data,
    cols = -iteration,
    names_to = "parameter",
    values_to = "value"
  )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = iteration, y = value)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = ncol) +
    ggplot2::labs(
      title = "MCMC Trace Plots",
      subtitle = paste("Population Model -", fit$num_subjects, "subjects"),
      x = "Iteration",
      y = "Parameter Value"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}


#' Plot Autocorrelation Functions
#'
#' Creates autocorrelation plots for population-level parameters.
#'
#' @param fit A \code{population_fit} object
#' @param parameters Character vector of parameter names to plot. If NULL,
#'   plots all population parameters.
#' @param max_lag Maximum lag to display (default: 50)
#' @param ncol Number of columns in plot grid (default: 2)
#'
#' @return A ggplot2 object
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_segment facet_wrap labs theme_minimal
#' @export
plot_acf <- function(fit, parameters = NULL, max_lag = 50, ncol = 2) {
  
  if (!inherits(fit, "population_fit")) {
    stop("fit must be a population_fit object")
  }
  
  pop_chain <- fit$population_chain
  
  # Remove iteration column
  if ("iteration" %in% names(pop_chain)) {
    pop_chain <- pop_chain[, names(pop_chain) != "iteration", drop = FALSE]
  }
  
  if (is.null(parameters)) {
    parameters <- names(pop_chain)
  } else {
    parameters <- intersect(parameters, names(pop_chain))
    if (length(parameters) == 0) {
      stop("No valid parameters found")
    }
  }
  
  # Calculate ACF for each parameter
  acf_data <- lapply(parameters, function(param) {
    x <- pop_chain[[param]]
    acf_vals <- stats::acf(x, lag.max = max_lag, plot = FALSE)
    data.frame(
      parameter = param,
      lag = acf_vals$lag,
      acf = acf_vals$acf,
      stringsAsFactors = FALSE
    )
  })
  
  acf_df <- do.call(rbind, acf_data)
  
  # Create plot
  p <- ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0), color = "steelblue") +
    ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = ncol) +
    ggplot2::labs(
      title = "Autocorrelation Functions",
      subtitle = paste("Population Model -", fit$num_subjects, "subjects"),
      x = "Lag",
      y = "Autocorrelation"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}


#' Comprehensive Convergence Report
#'
#' Generates a comprehensive text report of convergence diagnostics.
#'
#' @param fit A \code{population_fit} object
#'
#' @return Invisibly returns the diagnostics data frame, prints report
#'
#' @export
convergence_report <- function(fit) {
  
  cat("\n")
  cat("=" %R% 70, "\n")
  cat("MCMC Convergence Diagnostics Report\n")
  cat("=" %R% 70, "\n\n")
  
  cat("Model:", fit$model, "\n")
  cat("Subjects:", fit$num_subjects, "\n")
  cat("MCMC samples:", nrow(fit$population_chain), "\n")
  cat("(after thinning and burnin)\n\n")
  
  # Calculate diagnostics
  diag <- population_diagnostics(fit)
  
  cat("Population-Level Parameters:\n")
  cat("-" %R% 70, "\n")
  
  # Format output
  for (i in seq_len(nrow(diag))) {
    row <- diag[i, ]
    
    cat(sprintf("%-20s", row$parameter))
    cat(sprintf("  Mean: %8.3f", row$mean))
    cat(sprintf("  SD: %8.3f", row$sd))
    cat(sprintf("  ESS: %6.0f", row$ess))
    cat(sprintf("  ACF[1]: %5.3f", row$lag1_acf))
    
    # Convergence flag
    if (row$converged) {
      cat("  [OK]\n")
    } else {
      cat("  [LOW ESS]\n")
    }
  }
  
  cat("\n")
  cat("Convergence Summary:\n")
  cat("-" %R% 70, "\n")
  
  n_converged <- sum(diag$converged)
  n_total <- nrow(diag)
  pct_converged <- 100 * n_converged / n_total
  
  cat(sprintf("Parameters with ESS > 100: %d / %d (%.1f%%)\n",
              n_converged, n_total, pct_converged))
  
  min_ess <- min(diag$ess, na.rm = TRUE)
  max_ess <- max(diag$ess, na.rm = TRUE)
  cat(sprintf("ESS range: %.0f - %.0f\n", min_ess, max_ess))
  
  avg_acf <- mean(diag$lag1_acf, na.rm = TRUE)
  cat(sprintf("Average lag-1 autocorrelation: %.3f\n", avg_acf))
  
  cat("\n")
  cat("Interpretation:\n")
  cat("-" %R% 70, "\n")

  if (min_ess > 400) {
    cat("[OK] Excellent: All parameters have high effective sample sizes\n")
  } else if (min_ess > 100) {
    cat("[OK] Good: All parameters meet minimum ESS threshold\n")
  } else {
    cat("[WARNING] Some parameters have low ESS\n")
    cat("  Consider running longer or adjusting proposal variances\n")
  }

  if (avg_acf < 0.3) {
    cat("[OK] Good mixing: Low autocorrelation\n")
  } else if (avg_acf < 0.7) {
    cat("~ Moderate mixing\n")
  } else {
    cat("[WARNING] Poor mixing: High autocorrelation\n")
    cat("  Consider thinning more or adjusting proposal variances\n")
  }

  # Identifiability / space-filling check for bounded-prior variance components.
  # A parameter whose posterior fills its prior can pass ESS/ACF checks yet be
  # essentially uninformed by the data.
  ident <- tryCatch(identifiability_check(fit), error = function(e) NULL)
  if (!is.null(ident)) {
    flagged <- ident[ident$weak_identifiability %in% TRUE, , drop = FALSE]
    if (nrow(flagged) > 0) {
      cat("[WARNING] Weakly identified: posterior fills most of the prior range\n")
      for (i in seq_len(nrow(flagged))) {
        cat(sprintf("  %-16s posterior 95%% width covers %.0f%% of its prior\n",
                    flagged$parameter[i], 100 * flagged$prior_coverage[i]))
      }
      cat("  These parameters are likely poorly informed by the data;\n")
      cat("  their posteriors approach the prior. Interpret with caution.\n")
    }
  }

  cat("\n")
  cat("=" %R% 70, "\n\n")

  invisible(diag)
}


#' Flag Weakly Identified (Space-Filling) Parameters
#'
#' Detects parameters whose posterior spans most of their (bounded) prior
#' support -- the signature of a parameter that the data barely informs, so its
#' posterior simply reproduces the prior. This is a common, easily-missed
#' failure mode for variance components with little information in the data
#' (e.g. the pulse-to-pulse SD of pulse width in the population model), which
#' can wander across the whole prior range while every standard convergence
#' statistic (ESS, R-hat) still looks acceptable.
#'
#' For each parameter the check compares the width of the central 95\% posterior
#' interval to the width of the prior support. A large ratio (\code{coverage})
#' means the posterior fills the prior and the parameter is likely weakly
#' identified. Only parameters with a finite, bounded prior range contribute a
#' non-\code{NA} verdict. For a \code{population_fit} the uniform-prior variance
#' components (the \code{*_sd} columns) are mapped to their prior bounds
#' automatically. For a single-subject \code{pulse_fit} the chain is taken from
#' \code{\link{patient_chain}} and, when the spec uses the default
#' \code{sd_prior = "uniform"}, the \code{mass_sd}/\code{width_sd} bounds are
#' derived from the spec's Uniform(0, max) priors; under
#' \code{sd_prior = "half_cauchy"} the SD priors are unbounded and no bounds
#' can be derived.
#'
#' @param fit A \code{population_fit} or \code{pulse_fit} object, or a data
#'   frame / matrix of posterior draws (one column per parameter).
#' @param prior_bounds Optional named list giving the prior support for
#'   parameters as length-2 numeric vectors \code{c(lower, upper)}. Overrides /
#'   supplements the bounds derived automatically from a fitted model.
#' @param coverage_threshold Fraction of the prior range above which a parameter
#'   is flagged as weakly identified (default 0.5).
#'
#' @return A data frame with one row per parameter and columns:
#'   \item{parameter}{Parameter name}
#'   \item{post_mean}{Posterior mean}
#'   \item{post_sd}{Posterior standard deviation}
#'   \item{q2.5, q97.5}{Central 95\% posterior interval endpoints}
#'   \item{prior_range}{Width of the prior support (\code{NA} if unknown/unbounded)}
#'   \item{prior_coverage}{Posterior 95\% width divided by prior range}
#'   \item{weak_identifiability}{Logical; \code{TRUE} if \code{prior_coverage}
#'     exceeds \code{coverage_threshold}, \code{NA} if the prior range is unknown}
#'
#' @details
#' A high \code{prior_coverage} is a symptom, not a proof: strongly informative
#' data can still yield a wide posterior if the prior is tight. Read it alongside
#' the trace plot. A parameter that both fills its prior and drifts/space-fills
#' in the trace is very likely unidentified by the data at hand.
#'
#' @importFrom stats quantile sd
#' @export
identifiability_check <- function(fit, prior_bounds = NULL,
                                  coverage_threshold = 0.5) {

  # Resolve the chain and, for a fitted model, the uniform-prior bounds.
  if (inherits(fit, "population_fit")) {
    chain <- fit$population_chain
    auto_bounds <- .population_prior_bounds(fit)
    # Explicit prior_bounds take precedence over the auto-derived ones.
    user_bounds <- if (is.null(prior_bounds)) list() else prior_bounds
    prior_bounds <- utils::modifyList(auto_bounds, user_bounds)
  } else if (inherits(fit, "pulse_fit")) {
    chain <- patient_chain(fit)
    auto_bounds <- .pulse_prior_bounds(fit)
    user_bounds <- if (is.null(prior_bounds)) list() else prior_bounds
    prior_bounds <- utils::modifyList(auto_bounds, user_bounds)
  } else if (is.data.frame(fit) || is.matrix(fit)) {
    chain <- as.data.frame(fit)
    if (is.null(prior_bounds)) prior_bounds <- list()
  } else {
    stop(paste("fit must be a population_fit, a pulse_fit, or a",
               "data.frame/matrix of posterior draws"))
  }

  if ("iteration" %in% names(chain)) {
    chain <- chain[, names(chain) != "iteration", drop = FALSE]
  }

  params <- names(chain)
  results <- lapply(params, function(param) {

    x <- as.numeric(chain[[param]])
    x <- x[is.finite(x)]

    # Degrade gracefully on an all-NA / empty column: quantile() would error,
    # and because the caller wraps this in a single tryCatch, one bad column
    # would otherwise blank the entire report instead of just its own row.
    if (length(x) == 0) {
      qs <- c(NA_real_, NA_real_)
      post_range <- NA_real_
    } else {
      qs <- stats::quantile(x, c(0.025, 0.975), names = FALSE)
      post_range <- qs[2] - qs[1]
    }

    pb <- prior_bounds[[param]]
    if (!is.null(pb) && length(pb) == 2 && is.finite(post_range) &&
        all(is.finite(pb)) && (pb[2] - pb[1]) > 0) {
      prior_range <- pb[2] - pb[1]
      coverage    <- post_range / prior_range
      # 0.5 is a deliberately conservative trigger: a posterior 95% interval
      # spanning more than half the prior support is a strong "not learning
      # from data" signal. Tighten it if it over-flags in practice.
      weak        <- coverage > coverage_threshold
    } else {
      prior_range <- NA_real_
      coverage    <- NA_real_
      weak        <- NA
    }

    data.frame(
      parameter            = param,
      post_mean            = mean(x),
      post_sd              = stats::sd(x),
      q2.5                 = qs[1],
      q97.5                = qs[2],
      prior_range          = prior_range,
      prior_coverage       = coverage,
      weak_identifiability = weak,
      stringsAsFactors     = FALSE
    )
  })

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


# Derive prior bounds for a single-subject pulse_fit's chain columns. Only the
# pulse-to-pulse SDs have bounded priors, and only under sd_prior = "uniform"
# (Uniform(0, max)); a half-Cauchy prior has unbounded support, so the stored
# *_sd_max values are not prior bounds there and nothing can be derived.
.pulse_prior_bounds <- function(fit) {

  pri <- fit$spec$priors
  if (is.null(pri) || !isTRUE(as.logical(pri$uniform_sd_prior))) return(list())

  # chain column -> priors upper-bound name (Uniform(0, max))
  bound_names <- c(
    mass_sd  = "mass_sd_max",
    width_sd = "width_sd_max"
  )

  bounds <- list()
  for (param in names(bound_names)) {
    mx <- pri[[bound_names[[param]]]]
    if (!is.null(mx) && is.finite(mx) && mx > 0) {
      bounds[[param]] <- c(0, mx)
    }
  }
  bounds
}


# Map a population_fit's uniform-prior variance components to their prior
# support c(0, max). Returns an empty list when the spec is unavailable.
.population_prior_bounds <- function(fit) {

  pri <- fit$spec$population_priors
  if (is.null(pri)) return(list())

  # chain column -> population_priors upper-bound name (all are Uniform(0, max))
  bound_names <- c(
    mass_mean_sd  = "mass_mean_sd_max",
    width_mean_sd = "width_mean_sd_max",
    baseline_sd   = "baseline_sd_max",
    halflife_sd   = "halflife_sd_max",
    mass_sd       = "mass_sd_max",
    width_sd      = "width_sd_max"
  )

  bounds <- list()
  for (param in names(bound_names)) {
    mx <- pri[[bound_names[[param]]]]
    if (!is.null(mx) && is.finite(mx) && mx > 0) {
      bounds[[param]] <- c(0, mx)
    }
  }
  bounds
}


# Helper for string repetition
`%R%` <- function(x, n) {
  paste(rep(x, n), collapse = "")
}


#------------------------------------------------------------------------------#
#    End of file                                                               #
#------------------------------------------------------------------------------#
