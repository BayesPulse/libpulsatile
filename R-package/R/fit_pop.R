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


fit_pop_pulse <- function(data,
                          time,
                          num_patients, #Currently not used
                          location_prior = "strauss",
                          spec,
                          iters = 250000,
                          thin = 50,
                          burnin = as.integer(0.1 * iters),
                          use_tibble = TRUE, 
                          verbose = FALSE,
                          test_birthdeath = TRUE,
                          test_sd_masses = TRUE,
                          test_sd_width = FALSE,
                          test_s2s_sd_width = TRUE,
                          test_s2s_sd_mass = TRUE,
                          test_s2s_sd_baseline = TRUE,
                          test_s2s_sd_halflife = TRUE,
                          test_pop_means_width = TRUE,
                          test_pop_means_mass = TRUE,
                          test_pop_means_baseline = TRUE,
                          test_pop_means_halflife = TRUE,
                          test_fixeff_mass = TRUE,
                          test_fixeff_width = TRUE,
                          test_blhl = TRUE,
                          test_error = TRUE,
                          test_locations = TRUE,
                          test_masses = TRUE,
                          test_widths = TRUE,
                          test_tvarscale_mass = TRUE,
                          test_tvarscale_width = TRUE,
                          testMassVec = NULL,
                          testWidthVec = NULL,
                          testMKappaVec = NULL,
                          testWKappaVec = NULL,
                          testLocVec = NULL
                          ) {
  
  #if(num_patients != length(concs)) stop("Number of patients and concentration columns unequal")
  
  #indata <- list("time" = data$time, "concentrations" = as.matrix(data[,-1]))
  
  #stopifnot(is.numeric(indata[[time]]), is.numeric(indata[[conc]]),
   #         is.logical(use_tibble), is.logical(verbose))
  #if (burnin >= iters) stop("burnin >= iters")
  
  if(thin <= 0) stop("thin must be greater than 0")
  
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
  fit <- population_(data,
                     time,
                     location_prior,
                     priors,
                     proposalvariances,
                     startingvalues,
                     iters, thin, burnin, verbose,
                     pv_adjust_iter, pv_adjust_max_iter,
                     bivariate_pv_target_ratio, univariate_pv_target_ratio,
                     test_birthdeath,
                     test_sd_masses,
                     test_sd_width,
                     test_s2s_sd_width,
                     test_s2s_sd_mass,
                     test_s2s_sd_baseline,
                     test_s2s_sd_halflife,
                     test_pop_means_width,
                     test_pop_means_mass,
                     test_pop_means_baseline,
                     test_pop_means_halflife,
                     test_fixeff_mass,
                     test_fixeff_width,
                     test_blhl,
                     test_error,
                     test_locations,
                     test_masses,
                     test_widths,
                     test_tvarscale_mass,
                     test_tvarscale_width,
                     testMassVec,
                     testWidthVec,
                     testMKappaVec,
                     testWKappaVec,
                     testLocVec)



  population_chain <- as.data.frame(fit$pop_chain)

  patient_chains <- lapply(fit$patient_chains, as.data.frame);
  pulse_chains <- lapply(fit$pulse_chains, 
                        FUN = function(x){as.data.frame(do.call(rbind, x))})
  #patient_chains <- as.data.frame(fit$patient_chains)
  #pulse_chains <- as.data.frame(do.call(rbind, fit$pulse_chains))
  #pulse_chains <- fit$pulse_chains


  # Convert doubles to ints -- not strictly necessary -- consider.
  #pulse_chain$iteration        <- as.integer(pulse_chain$iteration)
  #pulse_chain$total_num_pulses <- as.integer(pulse_chain$total_num_pulses)
  #pulse_chain$pulse_num        <- as.integer(pulse_chain$pulse_num)
  #patient_chain$iteration       <- as.integer(patient_chain$iteration)
  #patient_chain$num_pulses      <- as.integer(patient_chain$num_pulses)

  #if (use_tibble) {
  #  patient_chain <- tibble::as_data_frame(patient_chain)
  #  pulse_chain  <- tibble::as_data_frame(pulse_chain)
  #}
  # temp return line

  rtn_obj <-
    structure(list("model"            = "population",
                   "call"             = Call,
                   "population_chain" = population_chain,
                   "patient_chain"    = patient_chains,
                   "pulse_chain"      = pulse_chains,
                   "data"             = data,
                   "time_range"       = fit$time_range,
                   "options"          = list(#"time"       = time,
                                             #"conc"       = conc,
                                             "thinning"   = thin,
                                             "iterations" = iters),
                   "spec"             = spec),
              class = "pop_pulse_fit")

  return(rtn_obj)
                   
}

