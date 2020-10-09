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

pop_spec <- function(prior_mass_mean,
                     prior_mass_var,
                     prior_mass_p2p_sd_var_max,
                     prior_mass_s2s_sd_var_max,
                     prior_width_mean,
                     prior_width_var,
                     prior_width_p2p_sd_var_max,
                     prior_width_s2s_sd_var_max,
                     prior_baseline_mean,
                     prior_baseline_var,
                     prior_baseline_s2s_sd_var_max,
                     prior_halflife_mean,
                     prior_halflife_var,
                     prior_halflife_s2s_sd_var_max,
                     sv_mass_mean,
                     sv_mass_p2p_sd_var,
                     sv_mass_s2s_sd_var,
                     sv_width_mean,
                     sv_width_p2p_sd_var,
                     sv_width_s2s_sd_var,
                     sv_baseline_mean,
                     sv_baseline_s2s_sd_var,
                     sv_halflife_mean,
                     sv_halflife_s2s_sd_var,
                     sv_error_var,
                     sv_error_alpha,
                     sv_error_beta,
                     sv_mean_pulse_count,
                     sv_strauss_repulsion,
                     sv_strauss_repulsion_range,
                     pv_sub_mass_mean,
                     pv_pop_mass_mean,
                     pv_ind_mass,
                     pv_sub_mass,
                     pv_sub_width_mean,
                     pv_pop_width_mean,
                     pv_ind_width,
                     pv_sub_width,
                     pv_pop_baseline_mean,
                     pv_sub_baseline,
                     pv_pop_halflife_mean,
                     pv_sub_halflife,
                     pv_pulse_location
) {
  
  pop_ps_obj <- 
    structure(
      list(#location_prior = location_prior_type,
           priors = list(mass_mean               = prior_mass_mean,
                         mass_var                = prior_mass_var,
                         mass_p2p_sd_var         = prior_mass_p2p_sd_var_max,
                         mass_s2s_sd_var         = prior_mass_s2s_sd_var_max,
                         width_mean              = prior_width_mean,
                         width_var               = prior_width_var,
                         width_p2p_sd_var        = prior_width_p2p_sd_var_max,
                         width_s2s_sd_var        = prior_width_s2s_sd_var_max,
                         baseline_mean           = prior_baseline_mean,
                         baseline_var            = prior_baseline_var,
                         baseline_s2s_sd_var     = prior_baseline_s2s_sd_var_max,
                         halflife_mean           = prior_halflife_mean,
                         halflife_var            = prior_halflife_var,
                         halflife_s2s_sd_var     = prior_halflife_s2s_sd_var_max),
           proposal_variances = list(sub_mass_mean     = pv_sub_mass_mean,
                                     pop_mass_mean     = pv_pop_mass_mean,
                                     ind_mass          = pv_ind_mass,
                                     sub_mass          = pv_sub_mass,
                                     sub_width_mean    = pv_sub_width_mean,
                                     pop_width_mean    = pv_pop_width_mean,
                                     ind_width         = pv_ind_width,
                                     sub_width         = pv_sub_width,
                                     pop_baseline_mean = pv_pop_baseline_mean,
                                     sub_baseline      = pv_sub_baseline,
                                     pop_halflife_mean = pv_pop_halflife_mean,
                                     sub_halflife      = pv_sub_halflife,
                                     pulse_location    = pv_pulse_location),
           starting_values = list(mass_mean           = sv_mass_mean,
                                  mass_p2p_sd_var     = sv_mass_p2p_sd_var,
                                  mass_s2s_sd_var     = sv_mass_s2s_sd_var,
                                  width_mean          = sv_width_mean,
                                  width_p2p_sd_var    = sv_width_p2p_sd_var,
                                  width_s2s_sd_var    = sv_width_s2s_sd_var,
                                  baseline_mean       = sv_baseline_mean,
                                  baseline_s2s_sd_var = sv_baseline_s2s_sd_var,
                                  halflife_mean       = sv_halflife_mean,
                                  halflife_s2s_sd_var = sv_halflife_s2s_sd_var,
                                  error_var           = sv_error_var,
                                  error_alpha         = sv_error_alpha,
                                  error_beta          = sv_error_beta,
                                  mean_pulse_count    = sv_mean_pulse_count,
                                  strauss_repulsion   = sv_strauss_repulsion,
                                  strauss_repulsion_range = sv_strauss_repulsion_range)),
           class = "pop_pulse_spec")
      
      return(pop_ps_obj)
}
