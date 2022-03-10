#' pulse_spec
#'
#' Generates a pulse_spec object -- the specification object required for
#' fitting a fit_pulse model.
#' 
#' @param location_prior_type Takes on two values: "order-statistic" and
#' "strauss". "order-statistic" uses every third order statistic of a Uniform
#' distribution for the pulse location prior and requires specification of the
#' prior parameter for mean pulse count ("prior_mean_pulse_count").
#' "strauss" uses the Strauss interacting point-process as a prior and requires
#' specification of "prior_mean_pulse_count", "prior_location_gamma", and
#' "prior_location_range".
#' @param prior_mass_mean Mean of Gaussian prior on patient mean pulse mass
#' @param prior_mass_var Variance of the Gaussian prior on patient mean pulse mass
#' @param prior_mass_sd Scale of the Half-Cauchy prior on the standard deviation of patient mean pulse mass
#' @param prior_width_mean Mean of Gaussian prior on patient mean pulse width
#' @param prior_width_var Variance of the Gaussian prior on patient mean pulse width
#' @param prior_width_sd Scale of the Half-Cauchy prior on the standard deviation of patient mean pulse width
#' @param prior_baseline_mean Mean of Gaussian prior on patient baseline concentration
#' @param prior_baseline_var Variance of the Gaussian prior on baseline concentration
#' @param prior_halflife_mean Mean of Gaussian prior on patient mean elimination half-life
#' @param prior_halflife_var Variance of the Gaussian prior on patient mean elimination half-life
#' @param prior_error_alpha Shape parameter of Gamma prior on patient error variance
#' @param prior_error_beta Scale parameter of Gamma prior on patient error variance
#' @param prior_location_gamma Repulsion parameter for Strauss prior on pulse count
#' @param prior_location_range Rate parameter for Strauss prior on pulse count
#' @param prior_mean_pulse_count Rate parameter for Strauss prior on pulse count
#' @param sv_mass_mean Starting value for patient mean pulse mass
#' @param sv_mass_sd Starting value for standard deviation of pulse-to-pulse variation of masses
#' @param sv_width_mean Starting value for patient mean pulse width
#' @param sv_width_sd Starting value for standard deviation of pulse-to-pulse variation of widths
#' @param sv_baseline_mean Starting value for patient mean baseline
#' @param sv_halflife_mean Starting value for patient mean half-life
#' @param sv_error_var Starting value for patient error variance
#' @param pv_pat_baseline Proposal variance for patient baseline concentration
#' @param pv_pat_halflife Proposal variance for patient elimination half-life
#' @param pv_pat_mass_mean Proposal variance for patient mean pulse mass
#' @param pv_pat_width_mean Proposal variance for patient mean pulse width
#' @param pv_ind_pulse_mass Proposal variance for individual pulse masses
#' @param pv_ind_pulse_width Proposal variance for individual pulse widths
#' @param pv_sd_pulse_mass Proposal variance for individual pulse mass standard deviation
#' @param pv_sd_pulse_width Proposal variance for individual pulse width standard deviation
#' @param pv_sdscale_pulse_mass Proposal variance for individual pulse mass distribution standard deviation scale
#' @param pv_sdscale_pulse_width Proposal variance for individual pulse width distribution standard deviation scale
#' @param pv_pulse_location Proposal variance for individual pulse locations
#' @export
#' @keywords pulse simulation
pulse_spec <-
  function(location_prior_type = c("strauss", "order-statistic"),
           prior_mass_mean        = 3.5,
           prior_mass_var         = 100,
           prior_width_mean       = 25,
           prior_width_var        = 1000,
           prior_baseline_mean    = 3,
           prior_baseline_var     = 100,
           prior_halflife_mean    = 45,
           prior_halflife_var     = 1000,
           prior_error_alpha      = 0.0001,
           prior_error_beta       = 0.0001,
           prior_location_gamma   = 0,
           prior_location_range   = 40,
           prior_mass_sd          = 1,
           prior_width_sd         = 1,
           prior_mean_pulse_count = 8,
           sv_mass_mean           = 3.5,
           sv_width_mean          = 25,
           sv_baseline_mean       = 3.0,
           sv_halflife_mean       = 45,
           sv_error_var           = 0.015,
           sv_mass_sd             = 0.5,
           sv_width_sd            = 5,
           pv_pat_baseline        = 0.015,
           pv_pat_halflife        = 25,
           pv_pat_mass_mean       = 0.5,
           pv_pat_width_mean      = 21,
           pv_ind_pulse_mass      = 1,
           pv_ind_pulse_width     = 50,
           pv_sd_pulse_mass       = 0.25,
           pv_sd_pulse_width      = 10,
           pv_sdscale_pulse_mass  = 4,
           pv_sdscale_pulse_width = 4,
           pv_pulse_location      = 40)
  {

    # TODO: Research better ways to do this range/valid-value checking.  Pretty
    # much all of the args need it.
    location_prior_type <- match.arg(location_prior_type)
    if (length(location_prior_type) > 1L) 
      stop(paste("location_prior_type is a required argument -- choose",
                 "'order-statistic' or 'strauss'"))
    if (prior_mean_pulse_count <= 0)
      stop(paste("prior_mean_pulse_count must be > 0."))

    if (location_prior_type == "strauss") {

      if (is.null(prior_location_gamma) | is.null(prior_location_range)) 
        stop(paste("prior_location_gamma and prior_location_range are required",
                   "arguments when location_prior_type == 'strauss'"))
      if (prior_location_gamma < 0 | prior_location_gamma > 1) 
        stop(paste("Invalid value for argument 'prior_location_gamma'; should",
                   "be in [0,1]"))  
      if (prior_location_range < 0)
        stop(paste("Invalid value for argument 'prior_location_range'; should",
                   "be >= 0"))  

    } else {

      if (!is.null(prior_location_gamma) | !is.null(prior_location_range))
        message(paste("When location_prior_type is set to 'order-statistic'",
                      "prior_location_gamma and prior_location_range are not used."))  

    }

    # Structure for single-subject, strauss prior model.
    # NOTE: sv's use std dev, while priors use variances (want this consistent
    # for any reason?)
    # NOTE: need more clear label for max_sd's 
    ps_obj <- 
      structure(
        list(location_prior = location_prior_type,
             priors = list(baseline_mean           = prior_baseline_mean,
                           baseline_variance       = prior_baseline_var,
                           halflife_mean           = prior_halflife_mean,
                           halflife_variance       = prior_halflife_var,
                           mass_mean               = prior_mass_mean,
                           mass_variance           = prior_mass_var,
                           width_mean              = prior_width_mean,
                           width_variance          = prior_width_var,
                           mass_sd_param           = prior_mass_sd,
                           width_sd_param          = prior_width_sd,
                           error_alpha             = prior_error_alpha,
                           error_beta              = prior_error_beta,
                           pulse_count             = prior_mean_pulse_count,
                           strauss_repulsion       = prior_location_gamma,
                           strauss_repulsion_range = prior_location_range),
             proposal_variances = list(mass_mean   = pv_pat_mass_mean,
                                       width_mean  = pv_pat_width_mean,
                                       mass_sd     = pv_sd_pulse_mass,
                                       width_sd    = pv_sd_pulse_width,
                                       baseline    = pv_pat_baseline,
                                       halflife    = pv_pat_halflife,
                                       location    = pv_pulse_location,
                                       pulse_mass  = pv_ind_pulse_mass,
                                       pulse_width = pv_ind_pulse_width,
                                       sdscale_pulse_mass  = pv_sdscale_pulse_mass ,
                                       sdscale_pulse_width = pv_sdscale_pulse_width),
             starting_values = list(baseline       = sv_baseline_mean,
                                    halflife       = sv_halflife_mean,
                                    errorsq        = sv_error_var,
                                    mass_mean      = sv_mass_mean,
                                    width_mean     = sv_width_mean,
                                    mass_sd        = sv_mass_sd,
                                    width_sd       = sv_width_sd)),
                class = "pulse_spec")

    return(ps_obj)

  }


#' @export 
print.pulse_spec <- function(x, ...) {

  cat("\nBayesian time-series analysis of pulsatile hormone data: 
      Model Specification Object\n\n")
  cat("Model type: Single Patient")
  cat("\n")
  cat("Pulse mass:\n")
  cat("  Fixed effect (mean)\n")
  cat("    prior mean =", x$priors$mass_mean, "\n") 
  cat("    prior variance =", x$priors$mass_variance, "\n") 
  cat("    starting value =", x$starting_values$mass_mean, "\n") 
  cat("    proposal variance =", x$proposal_variances$mass_mean, "\n")
  cat("  Fixed effect (SD)\n")
  cat("    prior parameter =", x$priors$mass_sd_param, "\n")
  cat("    starting value =", x$starting_values$mass_sd, "\n") 
  cat("    proposal variance =", x$proposal_variances$mass_sd, "\n")
  cat("  Random effects (individual pulses)\n")
  cat("    proposal variance =", x$proposal_variances$pulse_mass, "\n")
  cat("\n")
  cat("Pulse width:\n")
  cat("  Fixed effect (mean)\n")
  cat("    prior mean =", x$priors$width_mean, "\n") 
  cat("    prior variance =", x$priors$width_variance, "\n") 
  cat("    starting value =", x$starting_values$width_mean, "\n") 
  cat("    proposal variance =", x$proposal_variances$width_mean, "\n")
  cat("  Fixed effect (SD)\n")
  cat("    prior parameter =", x$priors$width_sd_param, "\n")
  cat("    starting value =", x$starting_values$width_sd, "\n") 
  cat("    proposal variance =", x$proposal_variances$width_sd, "\n")
  cat("  Random effects (individual pulses)\n")
  cat("    proposal variance =", x$proposal_variances$pulse_width, "\n")

}


#------------------------------------------------------------------------------#
#    End of file # End of file  # End of file  # End of file  # End of file    #
#------------------------------------------------------------------------------#
