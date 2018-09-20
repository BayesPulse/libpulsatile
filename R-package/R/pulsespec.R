
#-------------------------------------------------------------------------------
# Functions for creating a pulse model specification
#-------------------------------------------------------------------------------

#' pulse_spec
#'
#' Generates a pulse_spec object -- the specification object required for
#' fitting a fit_pulse model.
#'   
#' 
#' @param location_prior_type Takes on two values: "order-statistic" and
#' "strauss". "order-statistic" uses every third order statistic of a Uniform
#' distribution for the pulse location prior and requires specification of the
#' prior parameter for mean pulse count ("prior_mean_pulse_count").
#' "strauss" uses the Strauss interacting point-process as a prior and requires
#' specification of "prior_mean_pulse_count", "prior_location_gamma", and
#' "prior_location_range".
#' @param prior_mass_mean mass mean hyperparm
#' @param prior_mass_var mass variance hyperparm
#' @param prior_width_mean width mean hyperparm (on variance scale)
#' @param prior_width_var width variance hyperparm (on variance scale)
#' @param prior_baseline_mean mean of prior on baseline
#' @param prior_baseline_var variance of prior on baseline
#' @param prior_halflife_mean mean of prior on half-life
#' @param prior_halflife_var variance of prior on half-life
#' @param prior_error_alpha Gamma shape parameter
#' @param prior_error_beta Gamma rate parameter
#' @param prior_location_gamma placeholder
#' @param prior_location_range placeholder
#' @param prior_max_sd_mass placeholder
#' @param prior_max_sd_width placeholder
#' @param prior_mean_pulse_count placeholder
#' @param sv_mass_mean placeholder
#' @param sv_width_mean placeholder
#' @param sv_baseline_mean placeholder
#' @param sv_halflife_mean placeholder
#' @param sv_error_var placeholder
#' @param sv_mass_sd placeholder
#' @param sv_width_sd placeholder
#' @param pv_baseline placeholder
#' @param pv_halflife placeholder
#' @param pv_mean_pulse_mass placeholder
#' @param pv_mean_pulse_width placeholder
#' @param pv_indiv_pulse_mass placeholder
#' @param pv_indiv_pulse_width placeholder
#' @param pv_sd_pulse_mass placeholder
#' @param pv_sd_pulse_width Proposal variance of the SD of the pulse widths (pulse widths are on variance scale)
#' @param pv_sdscale_pulse_mass placeholder
#' @param pv_sdscale_pulse_width placeholder
#' @param pv_pulse_location placeholder
#' @export
#' @keywords pulse simulation
pulse_spec <-
  function(location_prior_type = c("order-statistic", "strauss"),
           prior_mass_mean        = 3.5,
           prior_mass_var         = 100,
           prior_width_mean       = 42,
           prior_width_var        = 1000,
           prior_baseline_mean    = 2.6,
           prior_baseline_var     = 100,
           prior_halflife_mean    = 45,
           prior_halflife_var     = 100,
           prior_error_alpha      = 0.0001,
           prior_error_beta       = 0.0001,
           prior_location_gamma   = NULL,
           prior_location_range   = NULL,
           prior_max_sd_mass      = 100,
           prior_max_sd_width     = 1000,
           prior_mean_pulse_count = 12,
           sv_mass_mean           = 3.5,
           sv_width_mean          = 42,
           sv_baseline_mean       = 2.6,
           sv_halflife_mean       = 45,
           sv_error_var           = 0.005,
           sv_mass_sd             = 1.6,
           sv_width_sd            = 35,
           pv_baseline            = 0.5,
           pv_halflife            = 45,
           pv_mean_pulse_mass     = 2,
           pv_mean_pulse_width    = 100,
           pv_indiv_pulse_mass    = 2,
           pv_indiv_pulse_width   = 100,
           pv_sd_pulse_mass       = 2,
           pv_sd_pulse_width      = 100,
           pv_sdscale_pulse_mass  = 1,
           pv_sdscale_pulse_width = 1,
           pv_pulse_location      = 10) 
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
                           mass_sdmax            = prior_max_sd_mass,
                           width_sdmax           = prior_max_sd_width,
                           error_alpha             = prior_error_alpha,
                           error_beta              = prior_error_beta,
                           pulse_count             = prior_mean_pulse_count,
                           strauss_repulsion       = prior_location_gamma,
                           strauss_repulsion_range = prior_location_range),
             proposal_variances = list(mass_mean   = pv_mean_pulse_mass,
                                       width_mean  = pv_mean_pulse_width,
                                       mass_sd     = pv_sd_pulse_mass,
                                       width_sd    = pv_sd_pulse_width,
                                       baseline    = pv_baseline,
                                       halflife    = pv_halflife,
                                       location    = pv_pulse_location,
                                       pulse_mass  = pv_indiv_pulse_mass,
                                       pulse_width = pv_indiv_pulse_width,
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
  cat("Model type:", paste0(x$model$model, "\n"))
  cat("Number of iterations:", 
      formatC(x$model$iterations, format = "d", big.mark = ","), "\n")
  cat("\n")
  cat("Pulse mass:\n")
  cat("  Fixed effect (mean)\n")
  cat("    prior mean =", x$priors$mass_mean, "\n") 
  cat("    prior variance =", x$priors$mass_variance, "\n") 
  cat("    starting value =", x$starting_values$mass_mean, "\n") 
  cat("    proposal variance =", x$proposal_variances$mass_mean, "\n")
  cat("  Fixed effect (SD)\n")
  cat("    prior maximum =", x$priors$mass_sdmax, "\n")
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
  cat("    prior maximum =", x$priors$width_sdmax, "\n")
  cat("    starting value =", x$starting_values$width_sd, "\n") 
  cat("    proposal variance =", x$proposal_variances$width_sd, "\n")
  cat("  Random effects (individual pulses)\n")
  cat("    proposal variance =", x$proposal_variances$pulse_width, "\n")

}


#------------------------------------------------------------------------------#
#    End of file # End of file  # End of file  # End of file  # End of file    #
#------------------------------------------------------------------------------#
