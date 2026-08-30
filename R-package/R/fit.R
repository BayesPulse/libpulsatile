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
#' @importFrom tibble as_tibble
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
                      #model_type = c("single-series"), #, "population", 
                      #               "single-series associational", 
                      #               "population associational"),
                      thin       = 50,
                      burnin     = as.integer(0.1 * iters),
                      use_tibble = TRUE,
                      verbose    = FALSE
                      ) {

  
  if (inherits(data, "pulse_sim")) {
    indata <- data$data
  } else {
    indata <- data.frame("time" = data[[time]], "concentration" = data[[conc]])
  }

  stopifnot(is.numeric(indata[[time]]), is.numeric(indata[[conc]]),
            is.logical(use_tibble), is.logical(verbose))
  if (burnin >= iters) stop("burnin >= iters")

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

  if (!inherits(spec, "pulse_spec")) {
    stop("spec is invalid -- see the fit_pulse() and pulse_spec()
         documentation.")
  }

  fit <- singlesubject_(indata$concentration,
                        indata$time,
                        spec$location_prior,
                        priors,
                        proposalvariances,
                        startingvalues,
                        iters, thin, burnin, verbose,
                        pv_adjust_iter, pv_adjust_max_iter,
                        bivariate_pv_target_ratio, univariate_pv_target_ratio)

  patient_chain <- as.data.frame(fit$patient_chain)
  pulse_chain   <- as.data.frame(do.call(rbind, fit$pulse_chains))

  # Convert doubles to ints -- not strictly necessary -- consider.
  pulse_chain$iteration        <- as.integer(pulse_chain$iteration)
  pulse_chain$total_num_pulses <- as.integer(pulse_chain$total_num_pulses)
  pulse_chain$pulse_num        <- as.integer(pulse_chain$pulse_num)
  patient_chain$iteration       <- as.integer(patient_chain$iteration)
  patient_chain$num_pulses      <- as.integer(patient_chain$num_pulses)

  if (use_tibble) {
    patient_chain <- tibble::as_tibble(patient_chain)
    pulse_chain  <- tibble::as_tibble(pulse_chain)
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
