#-------------------------------------------------------------------------------
# fit.R - Wrappers for estimating features of a pulsatile time-series dataset
#-------------------------------------------------------------------------------


#' fit_pulse
#'
#' Primary function for fitting deconvolution model for time series of pulsatile
#' hormone data.  Expects a tidy, long-format dataset (columns for time,
#' concentration, and subject ID; row for every concentration measurement).
#'
#' @param data A dataset containing \code{time} and \code{conc}.  Can
#'   also be a \code{pulse_sim} object.
#' @param time A string. Name of the time variable in \code{data}.
#' @param conc A string. Name of the hormone concentration variable in
#'   \code{data}.
#' @param subject_id A string. Name of the subject id variable in \code{data}.
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
#' @param verbose_iter Integer that sets diagnostic output to printout at the
#'   given interval if 'verbose' is \code{TRUE}. Default is 5000.
#' @param fix_params Character vector enabling parameters to fixed (i.e. not 
#'   estimated). Vector may include options:
#'   "pulse_count", "mass_mean", "width_mean", "mass_sd", 
#'   "width_sd", "bl_hl", "pulse_location", "error_var",
#'   "pulse_mass", "pulse_width", "pulse_mass_sdscale", 
#'   "pulse_width_sdscale" 
#'   or may be left empty to estimate all parameters. For patient level 
#'   parameters, values are fixed at the starting values specified by
#'   \code{pulse_spec()}. For pulse level parameters, vectors must be
#'   provided using the function arguments below to specify their fixed values.
#' @param pulse_masses Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse masses. 
#'   If used, \code{pulse_count} must also be fixed, and the length of the 
#'   vector must equal the sum of \code{pulse_counts}.
#' @param pulse_widths Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse widths. If
#'   used, \code{pulse_count} must also be fixed, and length of the vector 
#'   must equal the sum of \code{prior_mean_pulse_count}.
#' @param pulse_mass_sdscales Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse mass standard
#'   deviation scales. If used, \code{pulse_count} must also be fixed, and 
#'   length of the vector must equal the sum of \code{prior_mean_pulse_count}.
#' @param pulse_width_sdscales Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse width standard
#'   deviation scales. If used, \code{pulse_count} must also be fixed, and 
#'   length of the vector must equal the sum of \code{prior_mean_pulse_count}.
#' @param pulse_locations Numeric vector used in conjunction with 
#'   \code{fix_params} to choose fixed values for individual pulse locations. 
#'   If used, \code{pulse_count} must also be fixed, and 
#'   length of the vector must equal the sum of \code{prior_mean_pulse_count}.
#' @importFrom tibble as_data_frame
#' @keywords pulse fit
#' @examples
#' this_pulse <- simulate_pulse()
#' this_spec  <- pulse_spec()
#' #this_fit   <- fit_pulse(data = this_pulse, iters = 1000, thin = 10, 
#' #                        burnin = 100, spec = this_spec)
#' @export
fit_pulse <- function(data,
                      time = "time",
                      conc = "concentration",
                      subject_id = NULL,
                      spec,
                      iters = 250000,
                      thin       = 50,
                      burnin     = as.integer(0.1 * iters),
                      use_tibble = TRUE,
                      verbose    = FALSE,
                      verbose_iter = 5000,
                      fix_params = c("width_sd"),
                      pulse_masses = numeric(),
                      pulse_widths = numeric(),
                      pulse_mass_sdscales = numeric(),
                      pulse_width_sdscales = numeric(),
                      pulse_locations = numeric()) {
  
  if (all(class(data) == "pulse_sim")) {
    indata <- data$data
  } else {
    indata <- data.frame("time" = data[[time]], "concentration" = data[[conc]])
  }
  stopifnot(is.numeric(indata[[time]]), is.numeric(indata[[conc]]),
            is.logical(use_tibble), is.logical(verbose))
  
  # Prepare named list that is used to fix parameters
  param_list <- c("pulse_count", "mass_mean", "width_mean", "mass_prec",
                  "width_prec", "bl_hl", "pulse_location", "error_var",
                  "pulse_mass", "pulse_width", "pulse_mass_sdscale", 
                  "pulse_width_sdscale")
  if(!all(fix_params %in% param_list)) {stop("Invalid 'fix_params' element(s)")}
  fix_params <- as.list(param_list %in% fix_params)
  names(fix_params) <- param_list
  
  param_validation(data, time, conc, subject_id, spec, iters, thin,
                   burnin, use_tibble, verbose, verbose_iter,
                   fix_params,pulse_masses, pulse_widths,
                   pulse_mass_sdscales, pulse_width_sdscales,
                   pulse_locations)

  #---------------------------------------
  # TODO: Temporary work-arounds  (What todo for a non-temp work around??)
  #---------------------------------------
  # model type argument -- put in pulse_spec? 
  #model_type <- match.arg(model_type)
  # proposal variance tuning parameters -- in pulse_spec?
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
  #---------------------------------------

  # ideas via survival::coxph
  Call  <- match.call()
  arg_indx <- match(c("data", "time", "conc", "iters", "thin",
                      "spec"), 
                    names(Call), nomatch = 0)

  if (class(spec) != "pulse_spec") {
    stop("spec is invalid -- see the fit_pulse() and pulse_spec()
         documentation.")
  }

  fit <- singlesubject_(indata$concentration,
                        indata$time,
                        spec$location_prior,
                        priors,
                        proposalvariances,
                        startingvalues,
                        iters, thin, burnin, verbose, verbose_iter,
                        pv_adjust_iter, pv_adjust_max_iter,
                        bivariate_pv_target_ratio, univariate_pv_target_ratio,
                        fix_params, 
                        pulse_masses,
                        pulse_widths,
                        pulse_mass_sdscales,
                        pulse_width_sdscales,
                        pulse_locations)

  patient_chain <- as.data.frame(fit$patient_chain)
  pulse_chain   <- as.data.frame(do.call(rbind, fit$pulse_chains))

  # Convert doubles to ints -- not strictly necessary -- consider.
  pulse_chain$iteration        <- as.integer(pulse_chain$iteration)
  pulse_chain$total_num_pulses <- as.integer(pulse_chain$total_num_pulses)
  pulse_chain$pulse_num        <- as.integer(pulse_chain$pulse_num)
  patient_chain$iteration       <- as.integer(patient_chain$iteration)
  patient_chain$num_pulses      <- as.integer(patient_chain$num_pulses)

  if (use_tibble) {
    patient_chain <- tibble::as_data_frame(patient_chain)
    pulse_chain  <- tibble::as_data_frame(pulse_chain)
  }

  rtn_obj <- 
    structure(list("model"         = "single-subject",
                   "call"          = Call,
                   "patient_chain" = patient_chain,
                   "pulse_chain"   = pulse_chain,
                   "data"          = data,
                   "time_range"    = fit$time_range,
                   "options"       = list("time"       = time,
                                          "conc"       = conc,
                                          "thinning"   = thin,
                                          "iterations" = iters),
                   "spec"          = spec),
              class = "pulse_fit")

  return(rtn_obj)

}

param_validation <- function(data, time, conc, subject_id, spec, iters, thin,
                             burnin, use_tibble, verbose, verbose_iter,
                             fix_params,pulse_masses, pulse_widths,
                             pulse_mass_sdscales, pulse_width_sdscales,
                             pulse_locations) {
  
  if (burnin >= iters) stop("burnin >= iters")
  
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
     (length(pulse_locations) != spec$priors$pulse_count)) {
    stop("length(pulse_locations) != prior_mean_pulse_count")
  }
  if(fix_params[["pulse_mass"]] & 
     (length(pulse_masses) != spec$priors$pulse_count)) {
    stop("length(pulse_masses) != prior_mean_pulse_count")
  }
  if(fix_params[["pulse_width"]] & 
     (length(pulse_widths) != spec$priors$pulse_count)) {
    stop("length(pulse_widths) != prior_mean_pulse_count")
  }
  if(fix_params[["pulse_mass_sdscale"]] & 
     (length(pulse_mass_sdscales) != spec$priors$pulse_count)) {
    stop("length(pulse_mass_sdscales) != prior_mean_pulse_count")
  }
  if(fix_params[["pulse_width_sdscale"]] & 
     (length(pulse_width_sdscales) != spec$priors$pulse_count)) {
    stop("length(pulse_width_sdscales) != prior_mean_pulse_count")
  }
}
