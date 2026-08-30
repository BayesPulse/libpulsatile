#-------------------------------------------------------------------------------
# joint_spec.R - Functions for creating a joint hormone model specification
#-------------------------------------------------------------------------------

#' joint_spec
#'
#' Generates a joint_spec object -- the specification object required for
#' fitting a joint hormone model with fit_pulse_joint().
#'
#' The joint model fits two coupled pulsatile hormones where driver hormone
#' pulses influence the intensity of response hormone pulse occurrence.
#'
#' @param location_prior_type Takes on value "strauss" (only option currently).
#'   "strauss" uses the Strauss interacting point-process as a prior and requires
#'   specification of "prior_mean_pulse_count", "prior_location_gamma", and
#'   "prior_location_range".
#'
#' @section Driver Hormone Priors:
#' The \code{prior_driver_*} arguments set the prior distributions on the
#' driver hormone's pulse mass and width, baseline, half-life, model error,
#' pulse count, and pulse locations.
#'
#' @param prior_driver_mass_mean Prior mean for driver pulse mass
#' @param prior_driver_mass_var Prior variance for driver pulse mass
#' @param prior_driver_width_mean Prior mean for driver pulse width (on variance scale)
#' @param prior_driver_width_var Prior variance for driver pulse width (on variance scale)
#' @param prior_driver_baseline_mean Prior mean for driver baseline
#' @param prior_driver_baseline_var Prior variance for driver baseline
#' @param prior_driver_halflife_mean Prior mean for driver half-life
#' @param prior_driver_halflife_var Prior variance for driver half-life
#' @param prior_driver_error_alpha Gamma shape parameter for driver error
#' @param prior_driver_error_beta Gamma rate parameter for driver error
#' @param prior_driver_sd_mass Scale parameter of the half-Cauchy prior on the driver pulse-to-pulse SD of mass
#' @param prior_driver_sd_width Scale parameter of the half-Cauchy prior on the driver pulse-to-pulse SD of width
#' @param prior_driver_pulse_count Mean of Poisson prior on driver pulse count
#' @param prior_driver_location_gamma Strauss repulsion parameter for driver (0-1, 0=no repulsion)
#' @param prior_driver_location_range Strauss interaction range for driver (in time units)
#'
#' @section Response Hormone Priors:
#' The \code{prior_response_*} arguments mirror the driver priors for the
#' response hormone.
#'
#' @param prior_response_mass_mean Prior mean for response pulse mass
#' @param prior_response_mass_var Prior variance for response pulse mass
#' @param prior_response_width_mean Prior mean for response pulse width (on variance scale)
#' @param prior_response_width_var Prior variance for response pulse width (on variance scale)
#' @param prior_response_baseline_mean Prior mean for response baseline
#' @param prior_response_baseline_var Prior variance for response baseline
#' @param prior_response_halflife_mean Prior mean for response half-life
#' @param prior_response_halflife_var Prior variance for response half-life
#' @param prior_response_error_alpha Gamma shape parameter for response error
#' @param prior_response_error_beta Gamma rate parameter for response error
#' @param prior_response_sd_mass Scale parameter of the half-Cauchy prior on the response pulse-to-pulse SD of mass
#' @param prior_response_sd_width Scale parameter of the half-Cauchy prior on the response pulse-to-pulse SD of width
#' @param prior_response_pulse_count Mean of Poisson prior on response pulse count
#' @param prior_response_location_gamma Strauss repulsion parameter for response (0-1, 0=no repulsion)
#' @param prior_response_location_range Strauss interaction range for response (in time units)
#'
#' @section Association Priors:
#' The \code{prior_log_rho_*} and \code{prior_log_nu_*} arguments set normal
#' priors on the log of the coupling strength (\code{rho}) and coupling
#' temporal spread (\code{nu}).
#'
#' @param prior_log_rho_mean Prior mean for log(rho) - coupling strength parameter
#' @param prior_log_rho_var Prior variance for log(rho)
#' @param prior_log_nu_mean Prior mean for log(nu) - coupling temporal spread parameter
#' @param prior_log_nu_var Prior variance for log(nu)
#'
#' @section Driver Starting Values:
#' The \code{sv_driver_*} arguments set the MCMC starting values for the
#' driver hormone's common parameters.
#'
#' @param sv_driver_mass_mean Starting value for driver mean pulse mass
#' @param sv_driver_width_mean Starting value for driver mean pulse width
#' @param sv_driver_baseline_mean Starting value for driver baseline
#' @param sv_driver_halflife_mean Starting value for driver half-life
#' @param sv_driver_error_var Starting value for driver error variance
#' @param sv_driver_mass_sd Starting value for driver pulse-to-pulse SD of mass
#' @param sv_driver_width_sd Starting value for driver pulse-to-pulse SD of width
#'
#' @section Response Starting Values:
#' The \code{sv_response_*} arguments mirror the driver starting values for
#' the response hormone.
#'
#' @param sv_response_mass_mean Starting value for response mean pulse mass
#' @param sv_response_width_mean Starting value for response mean pulse width
#' @param sv_response_baseline_mean Starting value for response baseline
#' @param sv_response_halflife_mean Starting value for response half-life
#' @param sv_response_error_var Starting value for response error variance
#' @param sv_response_mass_sd Starting value for response pulse-to-pulse SD of mass
#' @param sv_response_width_sd Starting value for response pulse-to-pulse SD of width
#'
#' @section Association Starting Values:
#' \code{sv_rho} and \code{sv_nu} set the MCMC starting values for the
#' coupling parameters.
#'
#' @param sv_rho Starting value for rho (coupling strength)
#' @param sv_nu Starting value for nu (coupling temporal spread, is variance)
#'
#' @section Driver Proposal Variances:
#' The \code{pv_driver_*} arguments set the Metropolis-Hastings proposal
#' variances for the driver hormone's parameters.
#'
#' @param pv_driver_baseline Proposal variance for driver baseline
#' @param pv_driver_halflife Proposal variance for driver half-life
#' @param pv_driver_mean_pulse_mass Proposal variance for driver mean pulse mass
#' @param pv_driver_mean_pulse_width Proposal variance for driver mean pulse width
#' @param pv_driver_indiv_pulse_mass Proposal variance for driver individual pulse masses
#' @param pv_driver_indiv_pulse_width Proposal variance for driver individual pulse widths
#' @param pv_driver_sd_pulse_mass Proposal variance for driver pulse mass SD
#' @param pv_driver_sd_pulse_width Proposal variance for driver pulse width SD
#' @param pv_driver_sdscale_pulse_mass Proposal variance for driver pulse mass t-distribution scales
#' @param pv_driver_sdscale_pulse_width Proposal variance for driver pulse width t-distribution scales
#' @param pv_driver_pulse_location Proposal variance for driver pulse locations
#'
#' @section Response Proposal Variances:
#' The \code{pv_response_*} arguments mirror the driver proposal variances
#' for the response hormone.
#'
#' @param pv_response_baseline Proposal variance for response baseline
#' @param pv_response_halflife Proposal variance for response half-life
#' @param pv_response_mean_pulse_mass Proposal variance for response mean pulse mass
#' @param pv_response_mean_pulse_width Proposal variance for response mean pulse width
#' @param pv_response_indiv_pulse_mass Proposal variance for response individual pulse masses
#' @param pv_response_indiv_pulse_width Proposal variance for response individual pulse widths
#' @param pv_response_sd_pulse_mass Proposal variance for response pulse mass SD
#' @param pv_response_sd_pulse_width Proposal variance for response pulse width SD
#' @param pv_response_sdscale_pulse_mass Proposal variance for response pulse mass t-distribution scales
#' @param pv_response_sdscale_pulse_width Proposal variance for response pulse width t-distribution scales
#' @param pv_response_pulse_location Proposal variance for response pulse locations
#'
#' @section Association Proposal Variances:
#' \code{pv_log_rho} and \code{pv_log_nu} set the proposal variances for the
#' coupling parameters, on the log scale.
#'
#' @param pv_log_rho Proposal variance for log(rho)
#' @param pv_log_nu Proposal variance for log(nu)
#'
#' @param pulse_distribution Character, one of \code{"lognormal"} (default) or
#'   \code{"truncnorm"}, applied to both hormones. \code{"lognormal"} places the
#'   pulse mass/width random effects on the log scale and selects the log-scale
#'   numeric defaults for any argument left unspecified; \code{"truncnorm"} keeps
#'   the legacy natural-scale truncated-normal random effects and defaults.
#' @param sd_prior Character, one of \code{"uniform"} (default) or
#'   \code{"half_cauchy"}, applied to both hormones. Selects the prior on the
#'   pulse-to-pulse SDs; controls whether \code{prior_*_sd_mass}/
#'   \code{prior_*_sd_width} are the Uniform upper bound or the half-Cauchy scale.
#' @param student_t_pulses Logical. If \code{TRUE}, pulse mass and width random
#'   effects follow a Student-t distribution via a per-pulse t-scale
#'   (\code{tvarscale}, kappa) scale-mixture. If \code{FALSE} (default), the
#'   t-scale is fixed at 1 for every pulse and never sampled, giving Gaussian
#'   pulse random effects. Applies to both the driver and response hormones.
#'   Setting this to \code{FALSE} removes the weak-identifiability ridge between
#'   the pulse-to-pulse SD and the per-pulse t-scales, which can improve mixing
#'   of the SD parameters (notably the SD of pulse width).
#'
#' @section Breaking change (default parameterization):
#' The out-of-the-box defaults are now \code{{lognormal, uniform, gaussian}}
#' (both hormones), mirroring the single-subject paper model. Previously the
#' defaults were \code{{truncnorm, half_cauchy, student-t}}. Under
#' \code{"lognormal"} the mass/width prior means, starting values, and proposal
#' variances live on the LOG scale. Call \code{joint_spec(pulse_distribution =
#' "truncnorm", sd_prior = "half_cauchy", student_t_pulses = TRUE)} for the exact
#' pre-change behavior.
#'
#' @return A list of class \code{joint_spec} containing:
#'   \item{location_prior}{Type of location prior ("strauss")}
#'   \item{driver_priors}{List of prior parameters for driver hormone}
#'   \item{response_priors}{List of prior parameters for response hormone}
#'   \item{association_priors}{List of prior parameters for coupling}
#'   \item{proposal_variances}{List of proposal variances for all parameters}
#'   \item{driver_starting_values}{List of starting values for driver parameters}
#'   \item{response_starting_values}{List of starting values for response parameters}
#'   \item{association_starting_values}{List of starting values for coupling parameters}
#'
#' @export
#' @keywords pulse joint
#' @examples
#' # Create a joint specification with default values
#' joint_spec <- joint_spec()
#'
#' # Customize for specific hormones (e.g., LH driving FSH)
#' custom_spec <- joint_spec(
#'   prior_driver_mass_mean = 4.0,
#'   prior_response_mass_mean = 2.5,
#'   prior_log_rho_mean = 0.0,  # log(1) = neutral coupling
#'   prior_log_nu_mean = 3.0    # log(20) for temporal spread
#' )
joint_spec <- function(
  # Location prior
  location_prior_type = "strauss",

  # Driver priors
  prior_driver_mass_mean = NULL,
  prior_driver_mass_var = NULL,
  prior_driver_width_mean = NULL,
  prior_driver_width_var = NULL,
  prior_driver_baseline_mean = 2.6,
  prior_driver_baseline_var = 100,
  prior_driver_halflife_mean = 45,
  prior_driver_halflife_var = 100,
  prior_driver_error_alpha = 0.0001,
  prior_driver_error_beta = 0.0001,
  prior_driver_sd_mass = NULL,
  prior_driver_sd_width = NULL,
  prior_driver_pulse_count = 12,
  prior_driver_location_gamma = 0,
  prior_driver_location_range = 40,

  # Response priors
  prior_response_mass_mean = NULL,
  prior_response_mass_var = NULL,
  prior_response_width_mean = NULL,
  prior_response_width_var = NULL,
  prior_response_baseline_mean = 2.6,
  prior_response_baseline_var = 100,
  prior_response_halflife_mean = 45,
  prior_response_halflife_var = 100,
  prior_response_error_alpha = 0.0001,
  prior_response_error_beta = 0.0001,
  prior_response_sd_mass = NULL,
  prior_response_sd_width = NULL,
  prior_response_pulse_count = 12,
  prior_response_location_gamma = 0,
  prior_response_location_range = 40,

  # Association priors (on log scale)
  prior_log_rho_mean = 0.0,    # log(1) = neutral coupling
  prior_log_rho_var = 1.0,
  prior_log_nu_mean = 3.0,     # log(20) for temporal spread
  prior_log_nu_var = 1.0,

  # Driver starting values
  sv_driver_mass_mean = NULL,
  sv_driver_width_mean = NULL,
  sv_driver_baseline_mean = 2.6,
  sv_driver_halflife_mean = 45,
  sv_driver_error_var = 0.005,
  sv_driver_mass_sd = NULL,
  sv_driver_width_sd = NULL,

  # Response starting values
  sv_response_mass_mean = NULL,
  sv_response_width_mean = NULL,
  sv_response_baseline_mean = 2.6,
  sv_response_halflife_mean = 45,
  sv_response_error_var = 0.005,
  sv_response_mass_sd = NULL,
  sv_response_width_sd = NULL,

  # Association starting values
  sv_rho = 1.0,
  sv_nu = 20.0,

  # Driver proposal variances
  pv_driver_baseline = 0.02,
  pv_driver_halflife = 1.5,
  pv_driver_mean_pulse_mass = NULL,
  pv_driver_mean_pulse_width = NULL,
  pv_driver_indiv_pulse_mass = NULL,
  pv_driver_indiv_pulse_width = NULL,
  pv_driver_sd_pulse_mass = NULL,
  pv_driver_sd_pulse_width = NULL,
  pv_driver_sdscale_pulse_mass = 4,
  pv_driver_sdscale_pulse_width = 4,
  pv_driver_pulse_location = 65,

  # Response proposal variances
  pv_response_baseline = 0.02,
  pv_response_halflife = 1.5,
  pv_response_mean_pulse_mass = NULL,
  pv_response_mean_pulse_width = NULL,
  pv_response_indiv_pulse_mass = NULL,
  pv_response_indiv_pulse_width = NULL,
  pv_response_sd_pulse_mass = NULL,
  pv_response_sd_pulse_width = NULL,
  pv_response_sdscale_pulse_mass = 4,
  pv_response_sdscale_pulse_width = 4,
  pv_response_pulse_location = 65,

  # Association proposal variances (on log scale)
  pv_log_rho = 0.1,
  pv_log_nu = 0.1,

  # Random-effects parameterization (applies to both hormones)
  pulse_distribution = c("lognormal", "truncnorm"),
  sd_prior = c("uniform", "half_cauchy"),
  student_t_pulses = FALSE
) {

  # Input validation
  if (location_prior_type != "strauss") {
    stop("location_prior_type must be 'strauss' (only option currently supported)")
  }

  pulse_distribution <- match.arg(pulse_distribution)
  sd_prior           <- match.arg(sd_prior)

  if (!is.logical(student_t_pulses) || length(student_t_pulses) != 1L ||
      is.na(student_t_pulses)) {
    stop("student_t_pulses must be a single logical (TRUE or FALSE)")
  }

  if (prior_driver_pulse_count <= 0 || prior_response_pulse_count <= 0) {
    stop("prior pulse counts must be > 0")
  }

  # Parameterization flags carried through both priors lists into C++.
  lognormal_pulses <- identical(pulse_distribution, "lognormal")
  uniform_sd_prior <- identical(sd_prior, "uniform")

  # Conditional numeric defaults: driver and response each mirror the
  # single-subject (thesis) parameterization. Any argument left NULL is resolved
  # from the parameterization-specific table (log-scale for "lognormal", legacy
  # natural-scale for "truncnorm"); a user-supplied value overrides. Under the
  # log-scale default the pulse-to-pulse SD has a Uniform(0, 10) prior (the C++
  # default bound), so the starting SDs must sit inside it -- hence the converted
  # sv_*_width_sd = 0.7 rather than the natural-scale 35.
  ln_defaults <- list(
    mass_mean = 1.2,  mass_var = 25,  width_mean = 3.5,  width_var = 25,
    sd_mass = 10,     sd_width = 10,
    sv_mass_mean = 1.0, sv_width_mean = 3.0, sv_mass_sd = 0.5, sv_width_sd = 0.7,
    pv_mean_pulse_mass = 0.5,  pv_mean_pulse_width = 0.5,
    # Individual-pulse proposals are natural-scale random walks on the log-normal
    # pulse value, so they scale to the natural magnitude exp(mean) (width ~
    # exp(3.5) ~ 33 -> 9, mass -> 0.5), NOT log-scale sizes -- otherwise pulse
    # widths barely move and the width SD mixes slowly (see
    # benchmarks/sim_study_lognormal_report.md).
    pv_indiv_pulse_mass = 0.5,  pv_indiv_pulse_width = 9,
    pv_sd_pulse_mass = 0.1,     pv_sd_pulse_width = 0.1)
  tn_defaults <- list(
    mass_mean = 3.5,  mass_var = 100, width_mean = 42,   width_var = 1000,
    sd_mass = 5,      sd_width = 5,
    sv_mass_mean = 3.5, sv_width_mean = 42, sv_mass_sd = 1.6, sv_width_sd = 35,
    pv_mean_pulse_mass = 6,    pv_mean_pulse_width = 3700,
    pv_indiv_pulse_mass = 1,   pv_indiv_pulse_width = 15000,
    pv_sd_pulse_mass = 4.5,    pv_sd_pulse_width = 4000)
  d <- if (lognormal_pulses) ln_defaults else tn_defaults
  r <- function(value, name) if (is.null(value)) d[[name]] else value

  prior_driver_mass_mean   <- r(prior_driver_mass_mean,   "mass_mean")
  prior_driver_mass_var    <- r(prior_driver_mass_var,    "mass_var")
  prior_driver_width_mean  <- r(prior_driver_width_mean,  "width_mean")
  prior_driver_width_var   <- r(prior_driver_width_var,   "width_var")
  prior_driver_sd_mass     <- r(prior_driver_sd_mass,     "sd_mass")
  prior_driver_sd_width    <- r(prior_driver_sd_width,    "sd_width")
  prior_response_mass_mean  <- r(prior_response_mass_mean,  "mass_mean")
  prior_response_mass_var   <- r(prior_response_mass_var,   "mass_var")
  prior_response_width_mean <- r(prior_response_width_mean, "width_mean")
  prior_response_width_var  <- r(prior_response_width_var,  "width_var")
  prior_response_sd_mass    <- r(prior_response_sd_mass,    "sd_mass")
  prior_response_sd_width   <- r(prior_response_sd_width,   "sd_width")
  sv_driver_mass_mean   <- r(sv_driver_mass_mean,   "sv_mass_mean")
  sv_driver_width_mean  <- r(sv_driver_width_mean,  "sv_width_mean")
  sv_driver_mass_sd     <- r(sv_driver_mass_sd,     "sv_mass_sd")
  sv_driver_width_sd    <- r(sv_driver_width_sd,    "sv_width_sd")
  sv_response_mass_mean  <- r(sv_response_mass_mean,  "sv_mass_mean")
  sv_response_width_mean <- r(sv_response_width_mean, "sv_width_mean")
  sv_response_mass_sd    <- r(sv_response_mass_sd,    "sv_mass_sd")
  sv_response_width_sd   <- r(sv_response_width_sd,   "sv_width_sd")
  pv_driver_mean_pulse_mass   <- r(pv_driver_mean_pulse_mass,   "pv_mean_pulse_mass")
  pv_driver_mean_pulse_width  <- r(pv_driver_mean_pulse_width,  "pv_mean_pulse_width")
  pv_driver_indiv_pulse_mass  <- r(pv_driver_indiv_pulse_mass,  "pv_indiv_pulse_mass")
  pv_driver_indiv_pulse_width <- r(pv_driver_indiv_pulse_width, "pv_indiv_pulse_width")
  pv_driver_sd_pulse_mass     <- r(pv_driver_sd_pulse_mass,     "pv_sd_pulse_mass")
  pv_driver_sd_pulse_width    <- r(pv_driver_sd_pulse_width,    "pv_sd_pulse_width")
  pv_response_mean_pulse_mass   <- r(pv_response_mean_pulse_mass,   "pv_mean_pulse_mass")
  pv_response_mean_pulse_width  <- r(pv_response_mean_pulse_width,  "pv_mean_pulse_width")
  pv_response_indiv_pulse_mass  <- r(pv_response_indiv_pulse_mass,  "pv_indiv_pulse_mass")
  pv_response_indiv_pulse_width <- r(pv_response_indiv_pulse_width, "pv_indiv_pulse_width")
  pv_response_sd_pulse_mass     <- r(pv_response_sd_pulse_mass,     "pv_sd_pulse_mass")
  pv_response_sd_pulse_width    <- r(pv_response_sd_pulse_width,    "pv_sd_pulse_width")

  # Under a Uniform(0, max) SD prior the starting SD must lie strictly inside the
  # support, or the chain begins outside its own prior (see
  # ss_draw_sdrandomeffects.h::parameter_support). Checked for the driver and
  # response, only under the uniform prior (no hard bound under the half-Cauchy).
  if (uniform_sd_prior) {
    if (sv_driver_mass_sd >= prior_driver_sd_mass)
      stop(sprintf(paste("sv_driver_mass_sd = %g is outside the Uniform(0, %g)",
                         "prior support; lower it or raise prior_driver_sd_mass."),
                   sv_driver_mass_sd, prior_driver_sd_mass))
    if (sv_driver_width_sd >= prior_driver_sd_width)
      stop(sprintf(paste("sv_driver_width_sd = %g is outside the Uniform(0, %g)",
                         "prior support; lower it or raise prior_driver_sd_width."),
                   sv_driver_width_sd, prior_driver_sd_width))
    if (sv_response_mass_sd >= prior_response_sd_mass)
      stop(sprintf(paste("sv_response_mass_sd = %g is outside the Uniform(0, %g)",
                         "prior support; lower it or raise prior_response_sd_mass."),
                   sv_response_mass_sd, prior_response_sd_mass))
    if (sv_response_width_sd >= prior_response_sd_width)
      stop(sprintf(paste("sv_response_width_sd = %g is outside the Uniform(0, %g)",
                         "prior support; lower it or raise prior_response_sd_width."),
                   sv_response_width_sd, prior_response_sd_width))
  }

  # Validate Strauss prior parameters for driver
  if (prior_driver_location_gamma < 0 || prior_driver_location_gamma > 1) {
    stop("prior_driver_location_gamma must be in [0,1]")
  }
  if (prior_driver_location_range < 0) {
    stop("prior_driver_location_range must be >= 0")
  }

  # Validate Strauss prior parameters for response
  if (prior_response_location_gamma < 0 || prior_response_location_gamma > 1) {
    stop("prior_response_location_gamma must be in [0,1]")
  }
  if (prior_response_location_range < 0) {
    stop("prior_response_location_range must be >= 0")
  }

  # Validate positive parameters
  if (any(c(sv_driver_mass_sd, sv_driver_width_sd, sv_driver_error_var,
            sv_response_mass_sd, sv_response_width_sd, sv_response_error_var,
            sv_rho, sv_nu) <= 0)) {
    stop("All starting value SD, variance, and association parameters must be > 0")
  }

  # Create specification object
  spec_obj <- structure(
    list(
      location_prior = location_prior_type,

      # Driver priors
      driver_priors = list(
        baseline_mean = prior_driver_baseline_mean,
        baseline_variance = prior_driver_baseline_var,
        halflife_mean = prior_driver_halflife_mean,
        halflife_variance = prior_driver_halflife_var,
        mass_mean = prior_driver_mass_mean,
        mass_variance = prior_driver_mass_var,
        width_mean = prior_driver_width_mean,
        width_variance = prior_driver_width_var,
        mass_sd_param = prior_driver_sd_mass,
        width_sd_param = prior_driver_sd_width,
        mass_sd_max = prior_driver_sd_mass,
        width_sd_max = prior_driver_sd_width,
        error_alpha = prior_driver_error_alpha,
        error_beta = prior_driver_error_beta,
        pulse_count = prior_driver_pulse_count,
        strauss_repulsion = prior_driver_location_gamma,
        strauss_repulsion_range = prior_driver_location_range,
        lognormal_pulses = lognormal_pulses,
        uniform_sd_prior = uniform_sd_prior,
        student_t_pulses = student_t_pulses
      ),

      # Response priors
      response_priors = list(
        baseline_mean = prior_response_baseline_mean,
        baseline_variance = prior_response_baseline_var,
        halflife_mean = prior_response_halflife_mean,
        halflife_variance = prior_response_halflife_var,
        mass_mean = prior_response_mass_mean,
        mass_variance = prior_response_mass_var,
        width_mean = prior_response_width_mean,
        width_variance = prior_response_width_var,
        mass_sd_param = prior_response_sd_mass,
        width_sd_param = prior_response_sd_width,
        mass_sd_max = prior_response_sd_mass,
        width_sd_max = prior_response_sd_width,
        error_alpha = prior_response_error_alpha,
        error_beta = prior_response_error_beta,
        pulse_count = prior_response_pulse_count,
        strauss_repulsion = prior_response_location_gamma,
        strauss_repulsion_range = prior_response_location_range,
        lognormal_pulses = lognormal_pulses,
        uniform_sd_prior = uniform_sd_prior,
        student_t_pulses = student_t_pulses
      ),

      # Association priors
      association_priors = list(
        log_rho_mean = prior_log_rho_mean,
        log_rho_var = prior_log_rho_var,
        log_nu_mean = prior_log_nu_mean,
        log_nu_var = prior_log_nu_var
      ),

      # Proposal variances
      proposal_variances = list(
        # Driver
        driver_mass_mean = pv_driver_mean_pulse_mass,
        driver_width_mean = pv_driver_mean_pulse_width,
        driver_mass_sd = pv_driver_sd_pulse_mass,
        driver_width_sd = pv_driver_sd_pulse_width,
        driver_baseline = pv_driver_baseline,
        driver_halflife = pv_driver_halflife,
        driver_location = pv_driver_pulse_location,
        driver_pulse_mass = pv_driver_indiv_pulse_mass,
        driver_pulse_width = pv_driver_indiv_pulse_width,
        driver_sdscale_pulse_mass = pv_driver_sdscale_pulse_mass,
        driver_sdscale_pulse_width = pv_driver_sdscale_pulse_width,
        # Response
        response_mass_mean = pv_response_mean_pulse_mass,
        response_width_mean = pv_response_mean_pulse_width,
        response_mass_sd = pv_response_sd_pulse_mass,
        response_width_sd = pv_response_sd_pulse_width,
        response_baseline = pv_response_baseline,
        response_halflife = pv_response_halflife,
        response_location = pv_response_pulse_location,
        response_pulse_mass = pv_response_indiv_pulse_mass,
        response_pulse_width = pv_response_indiv_pulse_width,
        response_sdscale_pulse_mass = pv_response_sdscale_pulse_mass,
        response_sdscale_pulse_width = pv_response_sdscale_pulse_width,
        # Association (C++ code expects "rho" and "nu", not "log_rho" and "log_nu")
        rho = pv_log_rho,
        nu = pv_log_nu
      ),

      # Driver starting values
      driver_starting_values = list(
        mass_mean = sv_driver_mass_mean,
        width_mean = sv_driver_width_mean,
        baseline = sv_driver_baseline_mean,
        halflife = sv_driver_halflife_mean,
        mass_sd = sv_driver_mass_sd,
        width_sd = sv_driver_width_sd,
        errorsq = sv_driver_error_var
      ),

      # Response starting values
      response_starting_values = list(
        mass_mean = sv_response_mass_mean,
        width_mean = sv_response_width_mean,
        baseline = sv_response_baseline_mean,
        halflife = sv_response_halflife_mean,
        mass_sd = sv_response_mass_sd,
        width_sd = sv_response_width_sd,
        errorsq = sv_response_error_var
      ),

      # Association starting values
      association_starting_values = list(
        rho = sv_rho,
        nu = sv_nu
      )
    ),
    class = "joint_spec"
  )

  return(spec_obj)
}


#' @export
print.joint_spec <- function(x, ...) {
  cat("\nBayesian Joint Hormone Model Specification\n\n")
  cat("Location prior:", x$location_prior, "\n\n")

  cat("Driver Hormone:\n")
  cat("  Mean pulse count prior:", x$driver_priors$pulse_count, "\n")
  cat("  Pulse mass mean prior:", x$driver_priors$mass_mean,
      "(variance:", x$driver_priors$mass_variance, ")\n")
  cat("  Pulse width mean prior:", x$driver_priors$width_mean,
      "(variance:", x$driver_priors$width_variance, ")\n")
  cat("  Starting values: baseline =", x$driver_starting_values$baseline,
      ", halflife =", x$driver_starting_values$halflife, "\n\n")

  cat("Response Hormone:\n")
  cat("  Mean pulse count prior:", x$response_priors$pulse_count, "\n")
  cat("  Pulse mass mean prior:", x$response_priors$mass_mean,
      "(variance:", x$response_priors$mass_variance, ")\n")
  cat("  Pulse width mean prior:", x$response_priors$width_mean,
      "(variance:", x$response_priors$width_variance, ")\n")
  cat("  Starting values: baseline =", x$response_starting_values$baseline,
      ", halflife =", x$response_starting_values$halflife, "\n\n")

  cat("Association (Driver -> Response):\n")
  cat("  log(rho) prior: mean =", x$association_priors$log_rho_mean,
      ", variance =", x$association_priors$log_rho_var, "\n")
  cat("  log(nu) prior: mean =", x$association_priors$log_nu_mean,
      ", variance =", x$association_priors$log_nu_var, "\n")
  cat("  Starting values: rho =", x$association_starting_values$rho,
      ", nu =", x$association_starting_values$nu, "\n\n")

  invisible(x)
}


#------------------------------------------------------------------------------#
#    End of file                                                               #
#------------------------------------------------------------------------------#
