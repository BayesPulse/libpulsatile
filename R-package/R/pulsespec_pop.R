#-------------------------------------------------------------------------------
# Functions for creating a population pulse model specification
#-------------------------------------------------------------------------------

#' pop_spec
#' 
#' @param prior_mass_mean Mean of the Gaussian prior on population mean pulse mass
#' @param prior_mass_var Variance of the Gaussian prior on population mean pulse mass
#' @param prior_mass_s2s_sd Scale of the Half-Cauchy prior on the standard deviation of patient mean pulse mass
#' @param prior_mass_p2p_sd  Scale of the Half-Cauchy prior on the standard deviation of pulse-to-pulse variation of masses for a patient
#' @param prior_width_mean Mean of the Gaussian prior on population mean pulse width
#' @param prior_width_var  Variance of the Gaussian prior on population mean pulse width
#' @param prior_width_s2s_sd  Scale of the Half-Cauchy prior on the standard deviation of patient mean pulse width
#' @param prior_width_p2p_sd Scale of the Half-Cauchy prior on the standard deviation of pulse-to-pulse variation of widths for a patient
#' @param prior_baseline_mean Mean of the Gaussian prior on population mean baseline
#' @param prior_baseline_var Variance of the Gaussian prior on population mean baseline
#' @param prior_baseline_s2s_sd  Scale of the Half-Cauchy prior on the standard deviation of patient baselines
#' @param prior_halflife_mean Mean of the Gaussian prior on population mean half-life
#' @param prior_halflife_var Variance of the Gaussian prior on population mean half-life
#' @param prior_halflife_s2s_sd Scale of the Half-Cauchy prior on the standard deviation of patient half-lives
#' @param prior_error_alpha Shape parameter of Gamma prior on patient error variance
#' @param prior_error_beta Scale parameter of Gamma prior on patient error variance
#' @param prior_mean_pulse_count Rate parameter for Strauss prior on pulse count
#' @param prior_strauss_repulsion Repulsion parameter for Strauss prior on pulse count
#' @param prior_strauss_repulsion_range Rate parameter for Strauss prior on pulse count
#' @param sv_mass_mean Starting value for population mean pulse mass 
#' @param sv_mass_s2s_sd Starting value for standard deviation of patient mean pulse mass
#' @param sv_mass_p2p_sd Starting value for standard deviation of pulse-to-pulse variation of masses
#' @param sv_width_mean Starting value for population mean pulse width
#' @param sv_width_s2s_sd Starting value for standard deviation of patient mean pulse width
#' @param sv_width_p2p_sd Starting value for standard deviation of pulse-to-pulse variation of widths 
#' @param sv_baseline_mean Starting value for population mean baseline
#' @param sv_baseline_s2s_sd Starting value for standard deviation of patient mean baseline
#' @param sv_halflife_mean Starting value for population mean half-life
#' @param sv_halflife_s2s_sd Starting value for standard deviation of patient mean half-life
#' @param sv_error_var Starting value for patient error variance
#' @param pv_pat_mass_mean Proposal variance for patient mean pulse mass
#' @param pv_pop_mass_mean Proposal variance for population mean pulse mass
#' @param pv_ind_mass_sd Proposal variance for mass pulse-to-pulse standard deviation
#' @param pv_pat_mass_sd Proposal variance for standard deviation of patient mean pulse mass 
#' @param pv_pat_width_mean Proposal variance for patient mean pulse width
#' @param pv_pop_width_mean Proposal variance for population mean pulse width 
#' @param pv_ind_width_sd Proposal variance for width pulse-to-pulse standard deviation 
#' @param pv_pat_width_sd Proposal variance for standard deviation of patient mean pulse width 
#' @param pv_pop_baseline_mean Proposal variance for population mean baseline 
#' @param pv_pat_baseline Proposal variance for patient baselines 
#' @param pv_pat_baseline_sd Proposal variance for standard deviation of patient baseline 
#' @param pv_pop_halflife_mean Proposal variance for population mean half-life 
#' @param pv_pat_halflife Proposal variance for standard deviation of patient half-life
#' @param pv_pat_halflife_sd Proposal variance for standard deviation of patient half-life 
#' @param pv_pulse_location Proposal variance for pulse locations
#' @param pv_sdscale_pulse_mass Proposal variance for individual pulse mass distribution standard deviation scale 
#' @param pv_sdscale_pulse_width Proposal variance for individual pulse width distribution standard deviation scale 
#' 
#' 
#' @export

pop_spec <- function(prior_mass_mean = 3.5,
                     prior_mass_var = 100,
                     prior_mass_s2s_sd = 1,
                     prior_mass_p2p_sd = 1,
                     prior_width_mean = 25,
                     prior_width_var = 1000,
                     prior_width_s2s_sd = 1,
                     prior_width_p2p_sd = 1,
                     prior_baseline_mean = 3,
                     prior_baseline_var = 100,
                     prior_baseline_s2s_sd = 1,
                     prior_halflife_mean = 45,
                     prior_halflife_var = 1000,
                     prior_halflife_s2s_sd = 1,
                     prior_error_alpha = 0.0001,
                     prior_error_beta = 0.0001,
                     prior_mean_pulse_count = 8,
                     prior_strauss_repulsion = 0,
                     prior_strauss_repulsion_range = 40,
                     sv_mass_mean = 3.5,
                     sv_mass_s2s_sd = 0.5,
                     sv_mass_p2p_sd = 0.5,
                     sv_width_mean = 25,
                     sv_width_s2s_sd = 5,
                     sv_width_p2p_sd = 5,
                     sv_baseline_mean = 3.0,
                     sv_baseline_s2s_sd = 0.5,
                     sv_halflife_mean = 45,
                     sv_halflife_s2s_sd = 5,
                     sv_error_var = 0.015,
                     pv_pat_mass_mean = 0.5,
                     pv_pop_mass_mean = 0.5,
                     pv_ind_mass_sd = 1,
                     pv_pat_mass_sd = 1,
                     pv_pat_width_mean = 21,
                     pv_pop_width_mean = 21,
                     pv_ind_width_sd = 50,
                     pv_pat_width_sd = 25,
                     pv_pop_baseline_mean = 0.015,
                     pv_pat_baseline = 0.015,
                     pv_pat_baseline_sd = 1,
                     pv_pop_halflife_mean = 25,
                     pv_pat_halflife = 25,
                     pv_pat_halflife_sd = 1,
                     pv_pulse_location = 40,
                     pv_sdscale_pulse_mass = 4,
                     pv_sdscale_pulse_width = 4
) {
  
  pop_ps_obj <- 
    structure(
      list(#location_prior = location_prior_type,
           priors = list(mass_mean               = prior_mass_mean,
                         mass_var                = prior_mass_var,
                         mass_p2p_sd             = prior_mass_p2p_sd,
                         mass_s2s_sd             = prior_mass_s2s_sd,
                         width_mean              = prior_width_mean,
                         width_var               = prior_width_var,
                         width_p2p_sd            = prior_width_p2p_sd,
                         width_s2s_sd            = prior_width_s2s_sd,
                         baseline_mean           = prior_baseline_mean,
                         baseline_var            = prior_baseline_var,
                         baseline_s2s_sd         = prior_baseline_s2s_sd,
                         halflife_mean           = prior_halflife_mean,
                         halflife_var            = prior_halflife_var,
                         halflife_s2s_sd         = prior_halflife_s2s_sd,
                         error_alpha             = prior_error_alpha,
                         error_beta              = prior_error_beta,
                         mean_pulse_count        = prior_mean_pulse_count,
                         strauss_repulsion       = prior_strauss_repulsion,
                         strauss_repulsion_range = prior_strauss_repulsion_range),
           proposal_variances = list(pat_mass_mean       = pv_pat_mass_mean,
                                     pop_mass_mean       = pv_pop_mass_mean,
                                     ind_mass_sd         = pv_ind_mass_sd,
                                     pat_mass_sd         = pv_pat_mass_sd,
                                     pat_width_mean      = pv_pat_width_mean,
                                     pop_width_mean      = pv_pat_width_mean,
                                     ind_width_sd        = pv_ind_width_sd,
                                     pat_width_sd        = pv_pat_width_sd,
                                     pop_baseline_mean   = pv_pop_baseline_mean,
                                     pat_baseline        = pv_pat_baseline,
                                     pat_baseline_sd     = pv_pat_baseline_sd,
                                     pop_halflife_mean   = pv_pop_halflife_mean,
                                     pat_halflife        = pv_pat_halflife,
                                     pat_halflife_sd     = pv_pat_halflife_sd,
                                     pulse_location      = pv_pulse_location,
                                     sdscale_pulse_mass  =
                                       pv_sdscale_pulse_mass,
                                     sdscale_pulse_width =
                                       pv_sdscale_pulse_width),
           starting_values = list(mass_mean           = sv_mass_mean,
                                  mass_p2p_sd         = sv_mass_p2p_sd,
                                  mass_s2s_sd         = sv_mass_s2s_sd,
                                  width_mean          = sv_width_mean,
                                  width_p2p_sd        = sv_width_p2p_sd,
                                  width_s2s_sd        = sv_width_s2s_sd,
                                  baseline_mean       = sv_baseline_mean,
                                  baseline_s2s_sd     = sv_baseline_s2s_sd,
                                  halflife_mean       = sv_halflife_mean,
                                  halflife_s2s_sd     = sv_halflife_s2s_sd,
                                  error_var           = sv_error_var)),
           class = "pop_pulse_spec")
      
      return(pop_ps_obj)
}
