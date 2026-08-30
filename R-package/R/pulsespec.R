
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
#' @param prior_mass_mean Pulse mass prior mean. Under \code{pulse_distribution =
#'   "lognormal"} (default) this is the mean of log-mass; under \code{"truncnorm"}
#'   it is the natural-scale mass mean.
#' @param prior_mass_var Pulse mass prior variance (on the same scale as
#'   \code{prior_mass_mean}).
#' @param prior_width_mean Pulse width prior mean. Under \code{"lognormal"}
#'   (default) this is the mean of log-width; under \code{"truncnorm"} it is the
#'   natural-scale width mean (on the variance scale).
#' @param prior_width_var Pulse width prior variance (on the same scale as
#'   \code{prior_width_mean}).
#' @param prior_baseline_mean mean of prior on baseline
#' @param prior_baseline_var variance of prior on baseline
#' @param prior_halflife_mean mean of prior on half-life
#' @param prior_halflife_var variance of prior on half-life
#' @param prior_error_alpha Gamma shape parameter
#' @param prior_error_beta Gamma rate parameter
#' @param prior_location_gamma placeholder
#' @param prior_location_range placeholder
#' @param prior_sd_mass Pulse-to-pulse SD-of-mass prior parameter. Under
#'   \code{sd_prior = "uniform"} (default) it is the upper bound of the
#'   Uniform(0, .) prior; under \code{sd_prior = "half_cauchy"} it is the
#'   half-Cauchy scale. Defaults to 10 (lognormal) or 5 (truncnorm).
#' @param prior_sd_width Pulse-to-pulse SD-of-width prior parameter, interpreted
#'   like \code{prior_sd_mass}. Defaults to 10 (lognormal) or 5 (truncnorm).
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
#' @param pulse_distribution Character, one of \code{"lognormal"} (default) or
#'   \code{"truncnorm"}. \code{"lognormal"} places the pulse mass/width random
#'   effects on the log scale (\eqn{\log\theta \sim N(\mu, \sigma^2/\kappa)}) --
#'   the parameterization used by the Mulvahill thesis and Horton et al. (2017).
#'   \code{"truncnorm"} keeps the legacy natural-scale truncated-normal random
#'   effects. This choice also selects the paper-aligned (lognormal) vs. the
#'   legacy natural-scale (truncnorm) numeric defaults for any prior, starting
#'   value, or proposal variance left unspecified.
#' @param sd_prior Character, one of \code{"uniform"} (default) or
#'   \code{"half_cauchy"}. \code{"uniform"} gives the pulse-to-pulse SDs a
#'   Uniform(0, max) prior (papers' default); \code{"half_cauchy"} gives them the
#'   legacy half-Cauchy prior. Controls how \code{prior_sd_mass}/
#'   \code{prior_sd_width} are interpreted (uniform upper bound vs. Cauchy scale).
#' @param student_t_pulses Logical. If \code{TRUE}, pulse mass and width random
#'   effects follow a Student-t distribution via a per-pulse t-scale
#'   (\code{tvarscale}, kappa) scale-mixture. If \code{FALSE} (default), the
#'   t-scale is fixed at 1 for every pulse and never sampled, giving Gaussian
#'   pulse random effects. Setting this to \code{FALSE} removes the weak
#'   identifiability ridge between the pulse-to-pulse SD and the per-pulse
#'   t-scales, which can improve mixing of the SD parameters (notably the SD of
#'   pulse width) at the cost of no longer accommodating heavy-tailed pulses.
#'
#' @section Breaking change (default parameterization):
#' As of this version the out-of-the-box defaults are
#' \code{{lognormal, uniform, gaussian}} -- i.e. \code{pulse_distribution =
#' "lognormal"}, \code{sd_prior = "uniform"}, \code{student_t_pulses = FALSE} --
#' which reproduces the published single-subject (Mulvahill thesis) model.
#' Previously the defaults were \code{{truncnorm, half_cauchy, student-t}}. This
#' flips fitted results AND changes the meaning of \code{prior_width_mean} /
#' \code{sv_width_sd} etc.: under \code{"lognormal"} these live on the LOG scale
#' (e.g. mean log-width \eqn{\approx} 3.5), whereas under \code{"truncnorm"} they
#' are natural-scale (e.g. width mean 42 on the variance scale). To reproduce the
#' exact pre-change behavior, call \code{pulse_spec(pulse_distribution =
#' "truncnorm", sd_prior = "half_cauchy", student_t_pulses = TRUE)}. Any numeric
#' argument you supply explicitly always overrides the parameterization default.
#' @export
#' @keywords pulse simulation
pulse_spec <-
  function(location_prior_type = c("strauss", "order-statistic"),
           prior_mass_mean        = NULL,
           prior_mass_var         = NULL,
           prior_width_mean       = NULL,
           prior_width_var        = NULL,
           prior_baseline_mean    = 2.6,
           prior_baseline_var     = 100,
           prior_halflife_mean    = 45,
           prior_halflife_var     = 100,
           prior_error_alpha      = 0.0001,
           prior_error_beta       = 0.0001,
           prior_location_gamma   = 0,
           prior_location_range   = 40,
           prior_sd_mass          = NULL,
           prior_sd_width         = NULL,
           prior_mean_pulse_count = 12,
           sv_mass_mean           = NULL,
           sv_width_mean          = NULL,
           sv_baseline_mean       = 2.6,
           sv_halflife_mean       = 45,
           sv_error_var           = 0.005,
           sv_mass_sd             = NULL,
           sv_width_sd            = NULL,
           pv_baseline            = 0.02,
           pv_halflife            = 1.5,
           pv_mean_pulse_mass     = NULL,
           pv_mean_pulse_width    = NULL,
           pv_indiv_pulse_mass    = NULL,
           pv_indiv_pulse_width   = NULL,
           pv_sd_pulse_mass       = NULL,
           pv_sd_pulse_width      = NULL,
           pv_sdscale_pulse_mass  = 4,
           pv_sdscale_pulse_width = 4,
           pv_pulse_location      = 65,
           pulse_distribution     = c("lognormal", "truncnorm"),
           sd_prior               = c("uniform", "half_cauchy"),
           student_t_pulses       = FALSE)
  {

    # TODO: Research better ways to do this range/valid-value checking.  Pretty
    # much all of the args need it.
    location_prior_type <- match.arg(location_prior_type)
    pulse_distribution  <- match.arg(pulse_distribution)
    sd_prior            <- match.arg(sd_prior)
    if (length(location_prior_type) > 1L)
      stop(paste("location_prior_type is a required argument -- choose",
                 "'order-statistic' or 'strauss'"))
    if (!is.logical(student_t_pulses) || length(student_t_pulses) != 1L ||
        is.na(student_t_pulses))
      stop("student_t_pulses must be a single logical (TRUE or FALSE)")
    if (prior_mean_pulse_count <= 0)
      stop(paste("prior_mean_pulse_count must be > 0."))

    # Parameterization flags carried through the priors list into C++.
    lognormal_pulses <- identical(pulse_distribution, "lognormal")
    uniform_sd_prior <- identical(sd_prior, "uniform")

    # Conditional numeric defaults. Any argument left NULL is resolved from the
    # parameterization-specific default table; a user-supplied value overrides.
    # "truncnorm" reproduces the exact pre-change (natural-scale) spec object;
    # "lognormal" gives the paper log-scale defaults (thesis-matched). The
    # fixed-effect mean/SD proposal variances (pv_mean_*, pv_sd_*) are LOG-scale,
    # since those parameters live on the log scale. The INDIVIDUAL-pulse proposals
    # (pv_indiv_pulse_mass/width), however, are random walks on the NATURAL-scale
    # pulse value (log-normal draw), so they must be scaled to the natural pulse
    # magnitude exp(mean): mass ~ exp(1.2) ~ 3.3 -> pv 0.5; width ~ exp(3) ~ 20
    # -> pv 9. Shrinking these to log-scale sizes (~0.25) leaves individual widths
    # barely moving, so the pulse-to-pulse width SD mixes very slowly (a sim-study
    # confirmed width_sd ESS jumps ~3.6x, width_mean ESS ~5.7x, going 0.25 -> 9).
    ln_defaults <- list(
      prior_mass_mean  = 1.2,  prior_mass_var  = 25,
      prior_width_mean = 3.5,  prior_width_var = 25,
      prior_sd_mass    = 10,   prior_sd_width  = 10,
      sv_mass_mean     = 1.0,  sv_width_mean   = 3.0,
      sv_mass_sd       = 0.5,  sv_width_sd     = 0.7,
      pv_mean_pulse_mass  = 0.5,  pv_mean_pulse_width  = 0.5,
      pv_indiv_pulse_mass = 0.5,  pv_indiv_pulse_width = 9,
      pv_sd_pulse_mass    = 0.1,  pv_sd_pulse_width    = 0.1)
    tn_defaults <- list(
      prior_mass_mean  = 3.5,  prior_mass_var  = 100,
      prior_width_mean = 42,   prior_width_var = 1000,
      prior_sd_mass    = 5,    prior_sd_width  = 5,
      sv_mass_mean     = 3.5,  sv_width_mean   = 42,
      sv_mass_sd       = 1.6,  sv_width_sd     = 35,
      pv_mean_pulse_mass  = 6,    pv_mean_pulse_width  = 3700,
      pv_indiv_pulse_mass = 1,    pv_indiv_pulse_width = 15000,
      pv_sd_pulse_mass    = 4.5,  pv_sd_pulse_width    = 4000)
    defaults <- if (lognormal_pulses) ln_defaults else tn_defaults
    resolve <- function(value, name) if (is.null(value)) defaults[[name]] else value
    prior_mass_mean     <- resolve(prior_mass_mean,     "prior_mass_mean")
    prior_mass_var      <- resolve(prior_mass_var,      "prior_mass_var")
    prior_width_mean    <- resolve(prior_width_mean,    "prior_width_mean")
    prior_width_var     <- resolve(prior_width_var,     "prior_width_var")
    prior_sd_mass       <- resolve(prior_sd_mass,       "prior_sd_mass")
    prior_sd_width      <- resolve(prior_sd_width,      "prior_sd_width")
    sv_mass_mean        <- resolve(sv_mass_mean,        "sv_mass_mean")
    sv_width_mean       <- resolve(sv_width_mean,       "sv_width_mean")
    sv_mass_sd          <- resolve(sv_mass_sd,          "sv_mass_sd")
    sv_width_sd         <- resolve(sv_width_sd,         "sv_width_sd")
    pv_mean_pulse_mass  <- resolve(pv_mean_pulse_mass,  "pv_mean_pulse_mass")
    pv_mean_pulse_width <- resolve(pv_mean_pulse_width, "pv_mean_pulse_width")
    pv_indiv_pulse_mass <- resolve(pv_indiv_pulse_mass, "pv_indiv_pulse_mass")
    pv_indiv_pulse_width <- resolve(pv_indiv_pulse_width, "pv_indiv_pulse_width")
    pv_sd_pulse_mass    <- resolve(pv_sd_pulse_mass,    "pv_sd_pulse_mass")
    pv_sd_pulse_width   <- resolve(pv_sd_pulse_width,   "pv_sd_pulse_width")

    # Under a Uniform(0, max) SD prior the starting SD must lie strictly inside
    # the support, or the chain begins outside its own prior. (No hard bound
    # under the half-Cauchy prior, so only checked for the uniform case.)
    if (uniform_sd_prior) {
      if (sv_mass_sd >= prior_sd_mass)
        stop(sprintf(paste("sv_mass_sd = %g is outside the Uniform(0, %g) prior",
                           "support; lower it or raise prior_sd_mass."),
                     sv_mass_sd, prior_sd_mass))
      if (sv_width_sd >= prior_sd_width)
        stop(sprintf(paste("sv_width_sd = %g is outside the Uniform(0, %g) prior",
                           "support; lower it or raise prior_sd_width."),
                     sv_width_sd, prior_sd_width))
    }

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
                           mass_sd_param           = prior_sd_mass,
                           width_sd_param          = prior_sd_width,
                           mass_sd_max             = prior_sd_mass,
                           width_sd_max            = prior_sd_width,
                           error_alpha             = prior_error_alpha,
                           error_beta              = prior_error_beta,
                           pulse_count             = prior_mean_pulse_count,
                           strauss_repulsion       = prior_location_gamma,
                           strauss_repulsion_range = prior_location_range,
                           lognormal_pulses        = lognormal_pulses,
                           uniform_sd_prior        = uniform_sd_prior,
                           student_t_pulses        = student_t_pulses),
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
