#-------------------------------------------------------------------------------
# fit.R - Wrappers for estimating features of a pulsatile time-series dataset
#-------------------------------------------------------------------------------


#' fit_pulse
#'
#' Primary function for fitting deconvolution model for time series of pulsatile
#' hormone data.  Expects a tidy, long-format dataset (columns for time,
#' concentration, and subject ID; row for every concentration measurement).
#'
#' @useDynLib pulsatile decon_r_interface
#' @param .data A dataset containing \code{time} and \code{conc}.  Can
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
#' @import tibble
#' @keywords pulse fit
#' @examples
#' this_pulse <- simulate_pulse()
#' this_spec  <- pulse_spec()
#' this_fit   <- fit_pulse(.data = this_pulse, iters = 1000, thin = 10, 
#'                         burnin = 100, spec = this_spec)
#' @export
fit_pulse <- function(.data,
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

  
  if (all(class(.data) == "pulse_sim")) {
    indata <- .data$data
  } else {
    indata <- data.frame("time" = .data[[time]], "concentration" = .data[[conc]])
  }

  stopifnot(is.numeric(indata[[time]]), is.numeric(indata[[conc]]),
            is.logical(use_tibble), is.logical(verbose))
  if (burnin >= iters) stop("burnin >= iters")


  #model_type <- match.arg(model_type)
  # ideas via survival::coxph 
  Call  <- match.call()
  arg_indx <- match(c(".data", "time", "conc", "iters", "thin",
                      "spec"), 
                    names(Call), nomatch = 0)

  if (class(spec) != "pulse_spec") {
    stop("spec is invalid -- see the fit_pulse() and pulse_spec()
         documentation.")
  }

  # 
  #this_env <- environment()
  #list2env(spec, envir = this_env)

  rtn <- .Call(decon_r_interface, PACKAGE = "pulsatile",
               indata,
               "single-series",
               as.integer(thin),
               as.integer(burnin),
               as.integer(iters),
               as.integer(verbose),
               spec$strauss_location_prior,
               spec$priors$pulse_mass$mean,
               spec$priors$pulse_mass$var,
               spec$priors$pulse_width$mean,
               spec$priors$pulse_width$var,
               spec$priors$pulse_location$gamma,
               spec$priors$pulse_location$range,
               spec$priors$pulse_location$count,
               spec$priors$max_sd$mass,
               spec$priors$max_sd$width,
               spec$prior$baseline$mean,
               spec$prior$baseline$var,
               spec$prior$halflife$mean,
               spec$prior$halflife$var,
               spec$prior$error$alpha,
               spec$prior$error$beta,
               spec$starting_values$pulse_mass$mean,
               spec$starting_values$pulse_mass$sd,
               spec$starting_values$pulse_width$mean,
               spec$starting_values$pulse_width$sd,
               spec$starting_values$baseline$mean,
               spec$starting_values$halflife$mean,
               spec$starting_values$error$var,
               spec$proposal_variances$mean_pulse_mass,
               spec$proposal_variances$mean_pulse_width,
               spec$proposal_variances$indiv_pulse_mass,
               spec$proposal_variances$indiv_pulse_width,
               spec$proposal_variances$sd_pulse_mass,
               spec$proposal_variances$sd_pulse_width,
               spec$proposal_variances$pulse_location,
               spec$proposal_variances$baseline,
               spec$proposal_variances$halflife,
               spec$proposal_variances$etamass,
               spec$proposal_variances$etawidth)

  common_chain <- as.data.frame(rtn[[1]])
  pulse_chain  <- as.data.frame(do.call(rbind, rtn[[2]]))

  # 
  pulse_chain$iteration <- as.integer(pulse_chain$iteration)
  pulse_chain$total_num_pulses <- as.integer(pulse_chain$total_num_pulses)
  pulse_chain$pulse_num <- as.integer(pulse_chain$pulse_num)
  common_chain$iteration <- as.integer(common_chain$iteration)
  common_chain$num_pulses <- as.integer(common_chain$num_pulses)

  if (use_tibble) {
    common_chain = tibble::as_data_frame(common_chain)
    pulse_chain = tibble::as_data_frame(pulse_chain)
  }

  rtn_obj <- 
    structure(list("model"        = "single-subject",
                   "call"         = Call,
                   "common_chain" = common_chain,
                   "pulse_chain"  = pulse_chain,
                   "data"         = .data,
                   "options"      = list("time" = time,
                                         "conc" = conc,
                                         "thinning"   = thin, 
                                         "iterations" = iters),
                   "spec"         = spec),
              class = "pulse_fit")

  return(rtn_obj)

}
