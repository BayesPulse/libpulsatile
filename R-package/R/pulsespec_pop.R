#-------------------------------------------------------------------------------
# Functions for creating a population pulse model specification
#-------------------------------------------------------------------------------

#' pop_spec
#' 
#' @param prior_mass_mean Placeholder
#' @param prior_mass_variance Placeholder
#' 
#' 
#' @export

pop_spec <- function(prior_mass_mean, #4
                     prior_mass_var,
                     prior_mass_p2p_sd_var_max,
                     prior_mass_s2s_sd_var_max,
                     prior_width_mean, #5
                     prior_width_var,
                     prior_width_p2p_sd_var_max,
                     prior_width_s2s_sd_var_max,
                     prior_baseline_mean, #6
                     prior_baseline_var,
                     prior_baseline_s2s_sd_var_max,
                     prior_halflife_mean, #7
                     prior_halflife_var,
                     prior_halflife_s2s_sd_var_max,
                     prior_error_alpha, #8
                     prior_error_beta,
                     prior_mean_pulse_count, #9
                     sv_mass_mean, #10
                     sv_mass_p2p_sd_var,
                     sv_mass_s2s_sd_var,
                     sv_width_mean, #11
                     sv_width_p2p_sd_var,
                     sv_width_s2s_sd_var,
                     sv_baseline_mean, #12
                     sv_baseline_s2s_sd_var,
                     sv_halflife_mean, #13
                     sv_halflife_s2s_sd_var,
                     sv_error_var, #14
                     pv_sub_mass_mean, #15
                     pv_pop_mass_mean,
                     pv_ind_mass,
                     pv_sub_width_mean, #16
                     pv_pop_width_mean,
                     pv_ind_width,
                     pv_pop_baseline_mean, #17
                     pv_sub_baseline,
                     pv_pop_halflife_mean, #18
                     pv_sub_halflife,
                     pv_pulse_location #19 
) {
  
  pop_ps_obj <- 
    structure(
      list(#location_prior = location_prior_type,
           priors = list(mass_mean               = prior_mass_mean,
                         mass_var                = prior_mass_var,
                         mass_p2p_sd_var_max     = prior_mass_p2p_sd_var_max,
                         mass_s2s_sd_var_max     = prior_mass_s2s_sd_var_max,
                         width_mean              = prior_width_mean,
                         width_var               = prior_width_var,
                         width_p2p_sd_var_max    = prior_width_p2p_sd_var_max,
                         width_s2s_sd_var_max    = prior_width_s2s_sd_var_max,
                         baseline_mean           = prior_baseline_mean,
                         baseline_var            = prior_baseline_var,
                         baseline_s2s_sd_var_max = prior_baseline_s2s_sd_var_max,
                         halflife_mean           = prior_halflife_mean,
                         halflife_var            = prior_halflife_var,
                         halflife_s2s_sd_var_max = prior_halflife_s2s_sd_var_max,
                         error_alpha             = prior_error_alpha,
                         error_beta              = prior_error_beta,
                         mean_pulse_count        = prior_mean_pulse_count),
           proposal_variances = list(sub_mass_mean     = pv_sub_mass_mean,
                                     pop_mass_mean     = pv_pop_mass_mean,
                                     ind_mass          = pv_ind_mass,
                                     sub_width_mean    = pv_sub_width_mean,
                                     pop_width_mean    = pv_pop_width_mean,
                                     ind_width         = pv_ind_width,
                                     pop_baseline_mean = pv_pop_baseline_mean,
                                     sub_baseline      = pv_sub_baseline,
                                     pop_halflife_mean = pv_pop_halflife_mean,
                                     sub_halflife      = pv_sub_halflife,
                                     pulse_location    = pv_pulse_location),
           starting_values = list(mass_mean           = sv_mass_mean,
                                  mass_p2p_sd_var     = sv_mass_p2p_sd_var,
                                  mass_p2p_sd_var     = sv_mass_s2s_sd_var,
                                  width_mean          = sv_width_mean,
                                  width_p2p_sd_var    = sv_width_p2p_sd_var,
                                  width_s2s_sd_var    = sv_width_s2s_sd_var,
                                  baseline_mean       = sv_baseline_mean,
                                  baseline_s2s_sd_var = sv_baseline_s2s_sd_var,
                                  halflife_mean       = sv_halflife_mean,
                                  halflife_s2s_sd_var = sv_halflife_s2s_sd_var,
                                  error_var           = sv_error_var)),
           class = "pop_pulse_spec")
      
      return(pop_ps_obj)
}