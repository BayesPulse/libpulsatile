#-------------------------------------------------------------------------------------------
# fit_pop.R - Wrappers for estimating population features of a pulsatile time-series dataset
#-------------------------------------------------------------------------------------------

#' fit_pop_pulse
#' 
#' @param data Numeric matrix specifying pulsatile hormone concentration for each
#'   patient in population at each time of observation, with each patient as a 
#'   column and each row as an observation time. No additional columns should be
#'   present.
#' @param time Numeric vector for observation times. Must be the same length as
#'   the number of rows in \code{data}
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
#' @param verbose Prints diagnostics and estimates at iterations determined by
#'   'verbose_iters' argument. Default is \code{FALSE}.
#' @param verbose_patient Expands diagnostic output to include patient level
#'   parameters. Default is \code{FALSE}.
#' @param verbose_iters Integer that sets diagnostic output to printout at the
#'   given interval if 'verbose' is \code{TRUE}. Default is 5000.
#' @param fix_params Character vector enabling parameters to fixed (i.e. not 
#'   estimated). Vector may include options:
#'   "pop_mass_mean", "pop_width_mean", "pop_baseline_mean",
#'   "pop_halflife_mean", "mass_s2s_sd", "width_s2s_sd", "baseline_s2s_sd",
#'   "halflife_s2s_sd", "mass_p2p_sd", "width_p2p_sd", "pat_mass_mean",
#'   "pat_width_mean", "pat_bl_hl", "pat_error", "pulse_count", "pulse_location",
#'   "pulse_mass", "pulse_width", "pulse_mass_sdscale", "pulse_width_sdscale", 
#'   or may be left empty to estimate all parameters. For population level 
#'   parameters, values are fixed at the starting values specified by
#'   \code{pop_spec()}. For patient and pulse parameters, vectors must be
#'   provided using the function arguments below to specify their fixed values.
#' @param pat_mass_means Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for patient mass means. Must
#'   be the same length as the number of patients being analyzed (the number
#'   of columns passed to \code{data}).
#' @param pat_width_means Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for patient width means. Must
#'   be the same length as the number of patients being analyzed (the number
#'   of columns passed to \code{data}).
#' @param pat_baselines Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for patient baseline means. Must
#'   be the same length as the number of patients being analyzed (the number
#'   of columns passed to \code{data}).
#' @param pat_halflives Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for patient halflife means. Must
#'   be the same length as the number of patients being analyzed (the number
#'   of columns passed to \code{data}).
#' @param pulse_counts Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for patient pulse numbers. Must
#'   be the same length as the number of patients being analyzed (the number
#'   of columns passed to \code{data}).
#' @param pulse_masses Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse masses. 
#'   If used, \code{pulse_count} must also be fixed, and the length of the 
#'   vector must equal the sum of \code{pulse_counts}.
#' @param pulse_widths Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse widths. If
#'   used, \code{pulse_count} must also be fixed, and length of the vector 
#'   must equal the sum of \code{pulse_counts}.
#' @param pulse_mass_sdscale Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse mass standard
#'   deviation scales. If used, \code{pulse_count} must also be fixed, and 
#'   length of the vector must equal the sum of \code{pulse_counts}.
#' @param pulse_width_sdscale Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse width standard
#'   deviation scales. If used, \code{pulse_count} must also be fixed, and 
#'   length of the vector must equal the sum of \code{pulse_counts}.
#' @param pulse_locations Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse locations. 
#'   If used, \code{pulse_count} must also be fixed, and 
#'   length of the vector must equal the sum of \code{pulse_counts}.
#'   
#' @export
fit_pop_pulse <- function(data,
                          time,
                          location_prior = "strauss",
                          spec,
                          iters = 250000,
                          thin = 50,
                          burnin = as.integer(0.1 * iters),
                          use_tibble = TRUE, 
                          verbose = FALSE,
                          verbose_patient = FALSE,
                          verbose_iter = 5000,
                          fix_params = NULL,
                          pat_mass_means = numeric(),
                          pat_width_means = numeric(),
                          pat_baselines = numeric(),
                          pat_halflives = numeric(),
                          pulse_counts = numeric(),
                          pulse_masses = numeric(),
                          pulse_widths = numeric(),
                          pulse_mass_sdscales = numeric(),
                          pulse_width_sdscales = numeric(),
                          pulse_locations = numeric()
                          ) {
  
  # Prepare named list that is used to fix parameters
  param_list <- c("pop_mass_mean", "pop_width_mean", "pop_baseline_mean", 
                  "pop_halflife_mean", "mass_s2s_sd", "width_s2s_sd", 
                  "baseline_s2s_sd","halflife_s2s_sd", "mass_p2p_sd", 
                  "width_p2p_sd", "pat_mass_mean", "pat_width_mean", 
                  "pat_bl_hl", "pat_error", "pulse_count", "pulse_location", 
                  "pulse_mass", "pulse_width", "pulse_mass_sdscale", 
                  "pulse_width_sdscale")
  if(!all(fix_params %in% param_list)) {stop("Invalid 'fix_params' element(s)")}
  fix_params <- as.list(param_list %in% fix_params)
  names(fix_params) <- param_list
  
  pop_param_validation(data, time, location_prior, spec, iters, thin, burnin, 
                       use_tibble, verbose, verbose_patient, verbose_iter, 
                       fix_params, pat_mass_means, pat_width_means, 
                       pat_baselines, pat_halflives, pulse_counts, 
                       pulse_masses, pulse_widths, pulse_mass_sdscales, 
                       pulse_width_sdscales, pulse_locations)
  
  pv_adjust_iter <- 500
  pv_adjust_max_iter <- 25000
  univariate_pv_target_ratio <- 0.35
  bivariate_pv_target_ratio  <- 0.25
  
  # priors class type -- Definitely in pulse_spec
  priors            <- lapply(spec$priors, 
                              function(x) ifelse(is.null(x), NA, x))
  proposalvariances <- lapply(spec$proposal_variances, 
                              function(x) ifelse(is.null(x), NA, x))
  startingvalues    <- lapply(spec$starting_values, 
                              function(x) ifelse(is.null(x), NA, x))
  priors            <- structure(priors, class = "bp_priors")
  proposalvariances <- structure(proposalvariances, 
                                 class = "bp_proposalvariance")
  startingvalues    <- structure(startingvalues, 
                                 class = "bp_startingvals")
  
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
  fit <- population_(data, time, location_prior, priors, proposalvariances,
                     startingvalues, iters, thin, burnin, verbose,
                     verbose_patient, verbose_iter, pv_adjust_iter, 
                     pv_adjust_max_iter, bivariate_pv_target_ratio, 
                     univariate_pv_target_ratio, fix_params, pat_mass_means,
                     pat_width_means, pat_baselines, pat_halflives, 
                     pulse_counts, pulse_masses, pulse_widths, 
                     pulse_mass_sdscales, pulse_width_sdscales, pulse_locations)

  population_chain <- as.data.frame(fit$pop_chain)
  patient_chains <- lapply(fit$patient_chains, as.data.frame)
  pulse_chains <- lapply(fit$pulse_chains, 
                        FUN = function(x){as.data.frame(do.call(rbind, x))})
  
  names(patient_chains) <- colnames(data)
  names(pulse_chains) <- colnames(data)

  rtn_obj <-
    structure(list("model"            = "population",
                   "call"             = Call,
                   "population_chain" = population_chain,
                   "patient_chain"    = patient_chains,
                   "pulse_chain"      = pulse_chains,
                   "data"             = as.data.frame(cbind(time, data)),
                   "time_range"       = fit$time_range,
                   "options"          = list(#"time"       = time,
                                             #"conc"       = conc,
                                             "thinning"   = thin,
                                             "iterations" = iters),
                   "spec"             = spec),
              class = "pop_pulse_fit")

  return(rtn_obj)
                   
}

# Useful utility
`%notin%` <- Negate(`%in%`)

pop_param_validation <- function(data, time, location_prior, spec, iters, thin,
                                 burnin, use_tibble, verbose, verbose_patient,
                                 verbose_iter, fix_params, pat_mass_means,
                                 pat_width_means, pat_baselines,
                                 pat_halflives, pulse_counts, pulse_masses,
                                 pulse_widths, pulse_mass_sdscales, pulse_width_sdscales,
                                 pulse_locations) {
  
  param_list <- c("pop_mass_mean", "pop_width_mean", 
                  "pop_baseline_mean", "pop_halflife_mean",
                  "mass_s2s_sd", "width_s2s_sd",
                  "baseline_s2s_sd","halflife_s2s_sd",
                  "mass_p2p_sd", "width_p2p_sd", "pat_mass_mean", 
                  "pat_width_mean", "pat_bl_hl",  "pat_error", "pulse_count", 
                  "pulse_location", "pulse_mass", "pulse_width", 
                  "pulse_mass_sdscale", "pulse_width_sdscale")
  
  # Check general arguments
  if(thin <= 0) stop("thin must be greater than 0")
  if(iters <= 1) stop("iters must be greater than 1")
  #if(num_patients != length(concs)) stop("Number of patients and concentration columns unequal")
  #stopifnot(is.numeric(indata[[time]]), is.numeric(indata[[conc]]),
  #         is.logical(use_tibble), is.logical(verbose))
  if(burnin >= iters) stop("burnin >= iters")
  
  # Check patient parameter fixing arguments
  if(fix_params[["pat_mass_mean"]] & (length(pat_mass_means) == 0)) {
    stop("pat_mass_means vector must be provided if pat_mass_mean is fixed")
  }
  if(fix_params[["pat_width_mean"]] & (length(pat_width_means) == 0)) {
    stop("pat_width_means vector must be provided if pat_width_mean is fixed")
  }
  if(fix_params[["pat_bl_hl"]] & 
     (length(pat_baselines) == 0 | length(pat_halflives) == 0)) {
    stop("pat_baselines and pat_halflives vectors must be provided if pat_bl_hl is fixed")
  }
  
  if(fix_params[["pat_mass_mean"]] & 
     (length(pat_mass_means) != ncol(data))) {
    stop("length(pat_mass_means) != ncol(data)")
  }
  if(fix_params[["pat_width_mean"]] & 
     (length(pat_width_means) != ncol(data))) {
    stop("length(pat_width_means) != ncol(data)")
  }
  if(fix_params[["pat_bl_hl"]] & 
     (length(pat_baselines) != ncol(data))) {
    stop("length(pat_baselines) != ncol(data)")
  }
  if(fix_params[["pat_bl_hl"]] & 
     (length(pat_halflives) != ncol(data))) {
    stop("length(pat_halflives) != ncol(data)")
  }
  
  # Check pulse parameter fixing arguments
  if(any(fix_params[["pulse_location"]], fix_params[["pulse_mass"]],
         fix_params[["pulse_width"]], fix_params[["pulse_mass_sdscale"]],
         fix_params[["pulse_width_sdscale"]]) & 
     (!fix_params[["pulse_count"]])) {
    stop("pulse_count must be fixed if any other pulse parameters are fixed")
  }
  
  if(fix_params[["pulse_location"]] & (length(pulse_locations) == 0)) {
    stop("pulse_locations vector must be provided if pulse_location is fixed")
  }
  if(fix_params[["pulse_mass"]] & (length(pulse_masses) == 0)) {
    stop("pulse_masses vector must be provided if pulse_mass is fixed")
  }
  if(fix_params[["pulse_width"]] & (length(pulse_widths) == 0)) {
    stop("pulse_widths vector must be provided if pulse_width is fixed")
  }
  if(fix_params[["pulse_mass_sdscale"]] & (length(pulse_mass_sdscales) == 0)) {
    stop("pulse_mass_sdscales vector must be provided if pulse_mass_sdscale is fixed")
  }
  if(fix_params[["pulse_width_sdscale"]] & (length(pulse_width_sdscales) == 0)) {
    stop("pulse_width_sdscales vector must be provided if pulse_width_sdscale is fixed")
  }
  
  if(fix_params[["pulse_location"]] & 
     (length(pulse_locations) != sum(pulse_counts))) {
    stop("length(pulse_locations) != sum(pulse_counts)")
  }
  if(fix_params[["pulse_mass"]] & 
     (length(pulse_masses) != sum(pulse_counts))) {
    stop("length(pulse_masses) != sum(pulse_counts)")
  }
  if(fix_params[["pulse_width"]] & 
     (length(pulse_widths) != sum(pulse_counts))) {
    stop("length(pulse_widths) != sum(pulse_counts)")
  }
  if(fix_params[["pulse_mass_sdscale"]] & 
     (length(pulse_mass_sdscales) != sum(pulse_counts))) {
    stop("length(pulse_mass_sdscales) != sum(pulse_counts)")
  }
  if(fix_params[["pulse_width_sdscale"]] & 
     (length(pulse_width_sdscales) != sum(pulse_counts))) {
    stop("length(pulse_width_sdscales) != sum(pulse_counts)")
  }
  
}