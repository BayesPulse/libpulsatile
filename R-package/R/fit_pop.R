#-------------------------------------------------------------------------------------------
# fit_pop.R - Wrappers for estimating population features of a pulsatile time-series dataset
#-------------------------------------------------------------------------------------------

#' fit_pop_pulse
#' 
#' @param data Placeholder
#' @param num_patients Number of patients in the analysis
#' @param time A string. Name of the time variable in \code{data}
#' @param spec An object of class \code{spec}, created by
#'   \code{pulse_spec()}, specifying the priors, starting values, and proposal
#'   variances to use.
#' @param iters Number of iterations for MCMC
#' @param thin Thinning to apply to MCMC chains (i.e. Keep every 'thin'th
#'   sample).
#' @param burnin Burn-in to apply to MCMC chains (i.e. remove first 'burnin'
#'  samples). Applied prior to thinning.
#' @param use_tibble Return chains as tbl_df class data frames, from the tibble
#'   package.  Mostly used for the print.tbl_df method, which limits the rows and
#'   columns printed to those which fit in the console.
#' @param verbose Prints diagnostics and estimates every 5000th iteration.
#'   Default is \code{FALSE}.
#'   
#' @export


fit_pop_pulse <- function(data, #1
                          time = "time",
                          num_patients, #3
                          spec,
                          iters = 250000,
                          thin = 50,
                          burnin = as.integer(0.1 * iters),
                          use_tibble = TRUE, 
                          verbose = FALSE) {
  
  #if(num_patients != length(concs)) stop("Number of patients and concentration columns unequal")
  
  indata <- list("time" = data$time, "concentrations" = as.matrix(data[,-1]))
  
  #stopifnot(is.numeric(indata[[time]]), is.numeric(indata[[conc]]),
   #         is.logical(use_tibble), is.logical(verbose))
  #if (burnin >= iters) stop("burnin >= iters")
  
  pv_adjust_iter <- 500
  pv_adjust_max_iter <- 25000
  univariate_pv_target_ratio <- 0.35
  bivariate_pv_target_ratio  <- 0.25
  
  # priors class type -- Definitely in pulse_spec
  priors            <- lapply(spec$priors, function(x) ifelse(is.null(x), NA, x))
  proposalvariances <- lapply(spec$proposal_variances, function(x) ifelse(is.null(x), NA, x))
  startingvalues    <- lapply(spec$starting_values, function(x) ifelse(is.null(x), NA, x))
  priors            <- structure(priors, class = "bp_priors")
  proposalvariances <- structure(proposalvariances, class = "bp_proposalvariance")
  startingvalues    <- structure(startingvalues, class = "bp_startingvals")
  
  # ideas via survival::coxph
  Call  <- match.call()
  
  # not entirely sure what this line does... it is not used again
 # arg_indx <- match(c("data", "time", "conc", "iters", "thin",
                      #"spec"), 
                    #names(Call), nomatch = 0)
  
  if (class(spec) != "pop_pulse_spec") {
    stop("spec is invalid -- see the fit_pop_pulse() and pop_spec()
         documentation.")
  }
  
  # Call RCPP population function
  fit <- population_(indata$concentration,
                        indata$time,
                        "temp_location_prior",
                        priors,
                        proposalvariances,
                        startingvalues,
                        iters, thin, burnin, verbose,
                        pv_adjust_iter, pv_adjust_max_iter,
                        bivariate_pv_target_ratio, univariate_pv_target_ratio)
  
  # temp return line
  return(fit)
  
}

