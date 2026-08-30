#-------------------------------------------------------------------------------
# population_spec.R - Functions for creating a population model specification
#-------------------------------------------------------------------------------

#' population_spec
#'
#' Generates a population_spec object -- the specification object required for
#' fitting a population model with fit_pulse_population().
#'
#' This function specifies priors, starting values, and proposal variances for
#' hierarchical Bayesian modeling of pulsatile hormone data across multiple
#' subjects. The population model includes:
#' - Population-level means for pulse mass, width, baseline, and halflife
#' - Subject-to-subject variation parameters
#' - Pulse-to-pulse variation parameters
#'
#' @param location_prior_type Takes on two values: "strauss" (default) and
#'   "order-statistic". "strauss" uses the Strauss interacting point-process
#'   as a prior and requires specification of "prior_mean_pulse_count",
#'   "prior_location_gamma", and "prior_location_range".
#'   "order-statistic" uses every third order statistic of a Uniform
#'   distribution for the pulse location prior.
#' @param prior_mass_mean_mean Prior mean for population mean pulse mass
#' @param prior_mass_mean_var Prior variance for population mean pulse mass
#' @param prior_width_mean_mean Prior mean for population mean pulse width
#' @param prior_width_mean_var Prior variance for population mean pulse width
#' @param prior_baseline_mean_mean Prior mean for population baseline
#' @param prior_baseline_mean_var Prior variance for population baseline
#' @param prior_halflife_mean_mean Prior mean for population halflife
#' @param prior_halflife_mean_var Prior variance for population halflife
#' @param prior_mass_mean_sd_max Upper bound (uniform prior) for subject-to-subject SD of mean mass
#' @param prior_width_mean_sd_max Upper bound (uniform prior) for subject-to-subject SD of mean width
#' @param prior_baseline_sd_max Upper bound (uniform prior) for subject-to-subject SD of baseline
#' @param prior_halflife_sd_max Upper bound (uniform prior) for subject-to-subject SD of halflife
#' @param prior_mass_sd_max Upper bound (uniform prior) for pulse-to-pulse SD of mass
#' @param prior_width_sd_max Upper bound (uniform prior) for pulse-to-pulse SD of width
#' @param prior_error_alpha Gamma shape parameter for model error variance
#' @param prior_error_beta Gamma rate parameter for model error variance
#' @param prior_mean_pulse_count Mean of Poisson prior on pulse count
#' @param prior_location_gamma Strauss repulsion parameter (0-1, 0=no repulsion)
#' @param prior_location_range Strauss interaction range (in time units)
#' @param sv_mass_mean Starting value for population mean pulse mass
#' @param sv_width_mean Starting value for population mean pulse width
#' @param sv_baseline_mean Starting value for population baseline
#' @param sv_halflife_mean Starting value for population halflife
#' @param sv_mass_mean_sd Starting value for subject-to-subject SD of mean mass
#' @param sv_width_mean_sd Starting value for subject-to-subject SD of mean width
#' @param sv_baseline_sd Starting value for subject-to-subject SD of baseline
#' @param sv_halflife_sd Starting value for subject-to-subject SD of halflife
#' @param sv_mass_sd Starting value for pulse-to-pulse SD of mass
#' @param sv_width_sd Starting value for pulse-to-pulse SD of width
#' @param sv_error_var Starting value for model error variance
#' @param pv_mass_mean Proposal variance for population mean mass
#' @param pv_width_mean Proposal variance for population mean width
#' @param pv_baseline_mean Proposal variance for population baseline
#' @param pv_halflife_mean Proposal variance for population halflife
#' @param pv_mass_mean_sd Proposal variance for subject-to-subject SD of mean mass
#' @param pv_width_mean_sd Proposal variance for subject-to-subject SD of mean width
#' @param pv_baseline_sd Proposal variance for subject-to-subject SD of baseline
#' @param pv_halflife_sd Proposal variance for subject-to-subject SD of halflife
#' @param pv_mass_sd Proposal variance for pulse-to-pulse SD of mass
#' @param pv_width_sd Proposal variance for pulse-to-pulse SD of width
#' @param pv_subject_baseline Proposal variance for subject baselines
#' @param pv_subject_halflife Proposal variance for subject halflives
#' @param pv_subject_mean_pulse_mass Proposal variance for subject mean pulse mass
#' @param pv_subject_mean_pulse_width Proposal variance for subject mean pulse width
#' @param pv_indiv_pulse_mass Proposal variance for individual pulse masses
#' @param pv_indiv_pulse_width Proposal variance for individual pulse widths
#' @param pv_sdscale_pulse_mass Proposal variance for pulse mass t-distribution scales
#' @param pv_sdscale_pulse_width Proposal variance for pulse width t-distribution scales
#' @param pv_pulse_location Proposal variance for pulse locations
#' @param pulse_distribution Character, one of \code{"lognormal"} (default) or
#'   \code{"truncnorm"}. \code{"lognormal"} places the pulse mass/width random
#'   effects on the log scale -- the Horton et al. (2017) parameterization --
#'   and selects the paper-aligned (log-scale) numeric defaults for any argument
#'   left unspecified. \code{"truncnorm"} keeps the legacy natural-scale
#'   truncated-normal random effects and the legacy natural-scale defaults.
#' @param sd_prior Character, one of \code{"uniform"} (default) or
#'   \code{"half_cauchy"}. \strong{The population model only implements the
#'   Uniform SD prior}; \code{"half_cauchy"} is not implemented and, if supplied,
#'   is ignored with a warning (a Uniform prior is always used). Horton uses
#'   \code{"uniform"}.
#' @param student_t_pulses Logical. If \code{TRUE}, pulse mass and width random
#'   effects follow a Student-t distribution via a per-pulse t-scale
#'   (\code{tvarscale}, kappa) scale-mixture. If \code{FALSE} (default), the
#'   t-scale is fixed at 1 for every pulse and never sampled, giving Gaussian
#'   pulse random effects. Applies to every subject. Setting this to
#'   \code{FALSE} removes the weak-identifiability ridge between the
#'   pulse-to-pulse SD and the per-pulse t-scales, which can improve mixing of
#'   the SD parameters (notably the SD of pulse width).
#'
#' @section Breaking change (default parameterization):
#' The out-of-the-box defaults are now \code{{lognormal, uniform, gaussian}},
#' reproducing the published population model (Horton et al. 2017). Previously
#' the defaults were \code{{truncnorm, half_cauchy, student-t}}. This flips
#' fitted results and changes the meaning of \code{prior_width_mean_mean} /
#' \code{sv_width_sd} etc.: under \code{"lognormal"} these are on the LOG scale.
#' Call \code{population_spec(pulse_distribution = "truncnorm", sd_prior =
#' "half_cauchy", student_t_pulses = TRUE)} for the exact pre-change behavior.
#'
#' @return A list of class \code{population_spec} containing:
#'   \item{location_prior}{Type of location prior ("strauss" or "order-statistic")}
#'   \item{population_priors}{List of population-level prior parameters}
#'   \item{proposal_variances}{List of proposal variances for all parameters}
#'   \item{population_starting_values}{List of starting values for population parameters}
#'
#' @export
#' @keywords pulse population
#' @examples
#' # Create a population specification with default values
#' pop_spec <- population_spec()
#'
#' # Customize priors for a specific study
#' custom_spec <- population_spec(
#'   prior_mass_mean_mean = 4.0,
#'   prior_width_mean_mean = 50,
#'   prior_baseline_mean_mean = 3.0,
#'   prior_mean_pulse_count = 15
#' )
population_spec <- function(
  # Location prior
  location_prior_type = c("strauss", "order-statistic"),

  # Population-level priors (means)
  prior_mass_mean_mean = NULL,
  prior_mass_mean_var = NULL,
  prior_width_mean_mean = NULL,
  prior_width_mean_var = NULL,
  prior_baseline_mean_mean = 2.6,
  prior_baseline_mean_var = 100,
  prior_halflife_mean_mean = 45,
  # Baseline/half-life priors are axis-INDEPENDENT: they do not switch with
  # pulse_distribution (only the pulse mass/width priors do). Fixed defaults
  # matching pulse_spec()/joint_spec().
  prior_halflife_mean_var = 100,

  # Subject-to-subject variation priors (uniform upper bounds)
  prior_mass_mean_sd_max = NULL,
  prior_width_mean_sd_max = NULL,
  prior_baseline_sd_max = 3,
  prior_halflife_sd_max = 20,

  # Pulse-to-pulse variation priors (uniform upper bounds)
  prior_mass_sd_max = NULL,
  prior_width_sd_max = NULL,

  # Error and location priors
  prior_error_alpha = 0.0001,
  prior_error_beta = 0.0001,
  prior_mean_pulse_count = 12,
  prior_location_gamma = 0,
  prior_location_range = 40,

  # Population starting values
  sv_mass_mean = NULL,
  sv_width_mean = NULL,
  sv_baseline_mean = 2.6,
  sv_halflife_mean = 45,
  sv_mass_mean_sd = NULL,
  sv_width_mean_sd = NULL,
  sv_baseline_sd = 0.5,
  sv_halflife_sd = 5,
  sv_mass_sd = NULL,
  sv_width_sd = NULL,
  sv_error_var = 0.005,

  # Proposal variances - population parameters
  pv_mass_mean = NULL,
  pv_width_mean = NULL,
  pv_baseline_mean = 0.02,
  pv_halflife_mean = 1.5,
  pv_mass_mean_sd = NULL,
  pv_width_mean_sd = NULL,
  pv_baseline_sd = 0.5,
  pv_halflife_sd = 3,
  pv_mass_sd = NULL,
  pv_width_sd = NULL,

  # Proposal variances - subject-level parameters
  pv_subject_baseline = 0.02,
  pv_subject_halflife = 1.5,
  pv_subject_mean_pulse_mass = NULL,
  pv_subject_mean_pulse_width = NULL,
  pv_indiv_pulse_mass = NULL,
  pv_indiv_pulse_width = NULL,
  pv_sdscale_pulse_mass = 4,
  pv_sdscale_pulse_width = 4,
  pv_pulse_location = 65,

  # Random-effects parameterization
  pulse_distribution = c("lognormal", "truncnorm"),
  sd_prior = c("uniform", "half_cauchy"),
  student_t_pulses = FALSE
) {

  # Input validation
  location_prior_type <- match.arg(location_prior_type)
  pulse_distribution  <- match.arg(pulse_distribution)
  sd_prior            <- match.arg(sd_prior)

  if (!is.logical(student_t_pulses) || length(student_t_pulses) != 1L ||
      is.na(student_t_pulses)) {
    stop("student_t_pulses must be a single logical (TRUE or FALSE)")
  }

  if (length(location_prior_type) > 1L) {
    stop(paste("location_prior_type is a required argument -- choose",
               "'order-statistic' or 'strauss'"))
  }

  if (prior_mean_pulse_count <= 0) {
    stop("prior_mean_pulse_count must be > 0")
  }

  # The population model's SD draws are Uniform-only (see population.h); the
  # half-Cauchy prior is not implemented for it. Accept the argument (so existing
  # calls don't break) but warn that it is ignored and a Uniform prior is used.
  if (identical(sd_prior, "half_cauchy")) {
    warning("The population model only supports a Uniform SD prior; ",
            "sd_prior = \"half_cauchy\" is ignored and a Uniform prior is used.")
    # Reset so the uniform_sd_prior flag threaded into C++ matches the actual
    # (always Uniform) sampler behavior instead of reading FALSE.
    sd_prior <- "uniform"
  }

  # Parameterization flags carried through the priors list into C++.
  lognormal_pulses <- identical(pulse_distribution, "lognormal")
  uniform_sd_prior <- identical(sd_prior, "uniform")

  # Conditional numeric defaults (Horton-matched log-scale vs. legacy natural
  # scale). Any argument left NULL is resolved from the parameterization-specific
  # table; a user-supplied value overrides. "truncnorm" reproduces the exact
  # pre-change spec; "lognormal" uses log-scale priors, starting values within
  # the Uniform(0, 10) SD support, and log-scale proposal variances that mix.
  ln_defaults <- list(
    prior_mass_mean_mean = 1.0,  prior_mass_mean_var = 100,
    prior_width_mean_mean = 3.0, prior_width_mean_var = 100,
    prior_mass_mean_sd_max = 10, prior_width_mean_sd_max = 10,
    prior_mass_sd_max = 10,      prior_width_sd_max = 10,
    sv_mass_mean = 1.0,          sv_width_mean = 3.0,
    sv_mass_mean_sd = 0.5,       sv_width_mean_sd = 0.5,
    sv_mass_sd = 0.5,            sv_width_sd = 0.7,
    pv_mass_mean = 0.5,          pv_width_mean = 0.5,
    pv_mass_mean_sd = 0.1,       pv_width_mean_sd = 0.1,
    pv_mass_sd = 0.1,            pv_width_sd = 0.1,
    pv_subject_mean_pulse_mass = 0.5, pv_subject_mean_pulse_width = 0.5,
    # Individual-pulse proposals are natural-scale random walks on the log-normal
    # pulse value, so they are scaled to the natural magnitude exp(mean)
    # (mass ~ exp(1) ~ 2.7 -> 0.5; width ~ exp(3) ~ 20 -> 9), NOT to log-scale
    # sizes -- otherwise pulse widths barely move and the pulse-to-pulse width SD
    # mixes very slowly (see benchmarks/sim_study_lognormal_report.md).
    pv_indiv_pulse_mass = 0.5,   pv_indiv_pulse_width = 9)
  tn_defaults <- list(
    prior_mass_mean_mean = 3.5,  prior_mass_mean_var = 100,
    prior_width_mean_mean = 42,  prior_width_mean_var = 1000,
    prior_mass_mean_sd_max = 5,  prior_width_mean_sd_max = 30,
    prior_mass_sd_max = 5,       prior_width_sd_max = 30,
    sv_mass_mean = 3.5,          sv_width_mean = 42,
    sv_mass_mean_sd = 1.0,       sv_width_mean_sd = 10,
    sv_mass_sd = 1.6,            sv_width_sd = 15,
    pv_mass_mean = 6,            pv_width_mean = 3700,
    pv_mass_mean_sd = 2,         pv_width_mean_sd = 500,
    pv_mass_sd = 4.5,            pv_width_sd = 4000,
    pv_subject_mean_pulse_mass = 6, pv_subject_mean_pulse_width = 3700,
    pv_indiv_pulse_mass = 1,     pv_indiv_pulse_width = 15000)
  defaults <- if (lognormal_pulses) ln_defaults else tn_defaults
  resolve <- function(value, name) if (is.null(value)) defaults[[name]] else value
  prior_mass_mean_mean    <- resolve(prior_mass_mean_mean,    "prior_mass_mean_mean")
  prior_mass_mean_var     <- resolve(prior_mass_mean_var,     "prior_mass_mean_var")
  prior_width_mean_mean   <- resolve(prior_width_mean_mean,   "prior_width_mean_mean")
  prior_width_mean_var    <- resolve(prior_width_mean_var,    "prior_width_mean_var")
  prior_mass_mean_sd_max  <- resolve(prior_mass_mean_sd_max,  "prior_mass_mean_sd_max")
  prior_width_mean_sd_max <- resolve(prior_width_mean_sd_max, "prior_width_mean_sd_max")
  prior_mass_sd_max       <- resolve(prior_mass_sd_max,       "prior_mass_sd_max")
  prior_width_sd_max      <- resolve(prior_width_sd_max,      "prior_width_sd_max")
  sv_mass_mean            <- resolve(sv_mass_mean,            "sv_mass_mean")
  sv_width_mean           <- resolve(sv_width_mean,           "sv_width_mean")
  sv_mass_mean_sd         <- resolve(sv_mass_mean_sd,         "sv_mass_mean_sd")
  sv_width_mean_sd        <- resolve(sv_width_mean_sd,        "sv_width_mean_sd")
  sv_mass_sd              <- resolve(sv_mass_sd,              "sv_mass_sd")
  sv_width_sd             <- resolve(sv_width_sd,             "sv_width_sd")
  pv_mass_mean            <- resolve(pv_mass_mean,            "pv_mass_mean")
  pv_width_mean           <- resolve(pv_width_mean,           "pv_width_mean")
  pv_mass_mean_sd         <- resolve(pv_mass_mean_sd,         "pv_mass_mean_sd")
  pv_width_mean_sd        <- resolve(pv_width_mean_sd,        "pv_width_mean_sd")
  pv_mass_sd              <- resolve(pv_mass_sd,              "pv_mass_sd")
  pv_width_sd             <- resolve(pv_width_sd,             "pv_width_sd")
  pv_subject_mean_pulse_mass  <- resolve(pv_subject_mean_pulse_mass,  "pv_subject_mean_pulse_mass")
  pv_subject_mean_pulse_width <- resolve(pv_subject_mean_pulse_width, "pv_subject_mean_pulse_width")
  pv_indiv_pulse_mass     <- resolve(pv_indiv_pulse_mass,     "pv_indiv_pulse_mass")
  pv_indiv_pulse_width    <- resolve(pv_indiv_pulse_width,    "pv_indiv_pulse_width")
  
  # Validate Strauss prior parameters
  if (location_prior_type == "strauss") {
    if (is.null(prior_location_gamma) | is.null(prior_location_range)) {
      stop(paste("prior_location_gamma and prior_location_range are required",
                 "arguments when location_prior_type == 'strauss'"))
    }
    if (prior_location_gamma < 0 | prior_location_gamma > 1) {
      stop(paste("Invalid value for argument 'prior_location_gamma'; should",
                 "be in [0,1]"))
    }
    if (prior_location_range < 0) {
      stop(paste("Invalid value for argument 'prior_location_range'; should",
                 "be >= 0"))
    }
  } else {
    if (!is.null(prior_location_gamma) | !is.null(prior_location_range)) {
      message(paste("When location_prior_type is set to 'order-statistic',",
                    "prior_location_gamma and prior_location_range are not used."))
    }
  }
  
  # Validate positive parameters
  if (any(c(prior_mass_mean_sd_max, prior_width_mean_sd_max, 
            prior_baseline_sd_max, prior_halflife_sd_max,
            prior_mass_sd_max, prior_width_sd_max) <= 0)) {
    stop("All SD max parameters must be > 0")
  }
  
  if (any(c(sv_mass_mean_sd, sv_width_mean_sd, sv_baseline_sd,
            sv_halflife_sd, sv_mass_sd, sv_width_sd, sv_error_var) <= 0)) {
    stop("All starting value SD and variance parameters must be > 0")
  }

  # Each SD parameter has a Uniform(0, max) prior, so its starting value must lie
  # strictly inside that support -- otherwise the chain begins outside its own
  # prior and every proposal is rejected until it happens to land in range.
  sd_starts <- c(mass_mean_sd = sv_mass_mean_sd,
                 width_mean_sd = sv_width_mean_sd,
                 baseline_sd = sv_baseline_sd,
                 halflife_sd = sv_halflife_sd,
                 mass_sd = sv_mass_sd,
                 width_sd = sv_width_sd)
  sd_maxes  <- c(mass_mean_sd = prior_mass_mean_sd_max,
                 width_mean_sd = prior_width_mean_sd_max,
                 baseline_sd = prior_baseline_sd_max,
                 halflife_sd = prior_halflife_sd_max,
                 mass_sd = prior_mass_sd_max,
                 width_sd = prior_width_sd_max)
  bad <- which(sd_starts >= sd_maxes)
  if (length(bad) > 0) {
    stop(paste0("Starting value(s) outside the Uniform(0, max) prior support: ",
                paste(sprintf("sv_%s = %g >= max %g",
                              names(sd_starts)[bad], sd_starts[bad], sd_maxes[bad]),
                      collapse = "; "),
                ". Lower the starting value or raise the corresponding *_sd_max."))
  }

  # Create specification object
  spec_obj <- structure(
    list(
      location_prior = location_prior_type,
      
      # Population priors
      population_priors = list(
        mass_mean_mean = prior_mass_mean_mean,
        mass_mean_var = prior_mass_mean_var,
        width_mean_mean = prior_width_mean_mean,
        width_mean_var = prior_width_mean_var,
        baseline_mean_mean = prior_baseline_mean_mean,
        baseline_mean_var = prior_baseline_mean_var,
        halflife_mean_mean = prior_halflife_mean_mean,
        halflife_mean_var = prior_halflife_mean_var,
        mass_mean_sd_max = prior_mass_mean_sd_max,
        width_mean_sd_max = prior_width_mean_sd_max,
        baseline_sd_max = prior_baseline_sd_max,
        halflife_sd_max = prior_halflife_sd_max,
        mass_sd_max = prior_mass_sd_max,
        width_sd_max = prior_width_sd_max,
        error_alpha = prior_error_alpha,
        error_beta = prior_error_beta,
        pulse_count = prior_mean_pulse_count,
        strauss_repulsion = prior_location_gamma,
        strauss_repulsion_range = prior_location_range,
        lognormal_pulses = lognormal_pulses,
        uniform_sd_prior = uniform_sd_prior,
        student_t_pulses = student_t_pulses
      ),
      
      # Proposal variances
      proposal_variances = list(
        # Population-level (with pop_ prefix as expected by C++)
        pop_mass_mean_sd = pv_mass_mean_sd,
        pop_width_mean_sd = pv_width_mean_sd,
        pop_baseline_sd = pv_baseline_sd,
        pop_halflife_sd = pv_halflife_sd,
        pop_mass_sd = pv_mass_sd,
        pop_width_sd = pv_width_sd,
        # Subject-level
        baseline = pv_subject_baseline,
        halflife = pv_subject_halflife,
        mass_mean = pv_subject_mean_pulse_mass,
        width_mean = pv_subject_mean_pulse_width,
        pulse_mass = pv_indiv_pulse_mass,
        pulse_width = pv_indiv_pulse_width,
        sdscale_pulse_mass = pv_sdscale_pulse_mass,
        sdscale_pulse_width = pv_sdscale_pulse_width,
        location = pv_pulse_location
      ),
      
      # Population starting values
      population_starting_values = list(
        mass_mean = sv_mass_mean,
        width_mean = sv_width_mean,
        baseline_mean = sv_baseline_mean,
        halflife_mean = sv_halflife_mean,
        mass_mean_sd = sv_mass_mean_sd,
        width_mean_sd = sv_width_mean_sd,
        baseline_sd = sv_baseline_sd,
        halflife_sd = sv_halflife_sd,
        mass_sd = sv_mass_sd,
        width_sd = sv_width_sd,
        errorsq = sv_error_var
      )
    ),
    class = "population_spec"
  )
  
  return(spec_obj)
}


#' @export
print.population_spec <- function(x, ...) {
  cat("\nBayesian Population Model Specification\n\n")
  cat("Location prior:", x$location_prior, "\n")
  cat("Mean pulse count prior:", x$population_priors$pulse_count, "\n\n")
  
  cat("Population Parameters:\n")
  cat("  Pulse mass (population mean):\n")
  cat("    prior mean =", x$population_priors$mass_mean_mean, "\n")
  cat("    prior variance =", x$population_priors$mass_mean_var, "\n")
  cat("    starting value =", x$population_starting_values$mass_mean, "\n\n")
  
  cat("  Pulse width (population mean):\n")
  cat("    prior mean =", x$population_priors$width_mean_mean, "\n")
  cat("    prior variance =", x$population_priors$width_mean_var, "\n")
  cat("    starting value =", x$population_starting_values$width_mean, "\n\n")
  
  cat("Subject-to-Subject Variation:\n")
  cat("  Mean mass SD:\n")
  cat("    prior max =", x$population_priors$mass_mean_sd_max, "\n")
  cat("    starting value =", x$population_starting_values$mass_mean_sd, "\n\n")
  
  cat("Pulse-to-Pulse Variation:\n")
  cat("  Mass SD:\n")
  cat("    prior max =", x$population_priors$mass_sd_max, "\n")
  cat("    starting value =", x$population_starting_values$mass_sd, "\n\n")
}


#------------------------------------------------------------------------------#
#    End of file                                                               #
#------------------------------------------------------------------------------#
