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

pop_spec <- function(prior_mass_mean = 3.5,
                     prior_mass_var = 100,
                     prior_mass_p2p_sd_var_max = 1,
                     prior_mass_s2s_sd_var_max = 1,
                     prior_width_mean = 25,
                     prior_width_var = 1000,
                     prior_width_p2p_sd_var_max = 1,
                     prior_width_s2s_sd_var_max = 1,
                     prior_baseline_mean = 3,
                     prior_baseline_var = 100,
                     prior_baseline_s2s_sd_var_max = 1,
                     prior_halflife_mean = 45,
                     prior_halflife_var = 1000,
                     prior_halflife_s2s_sd_var_max = 1,
                     prior_error_alpha = 0.0001,
                     prior_error_beta = 0.0001,
                     sv_mass_mean = 3.5,
                     sv_mass_p2p_sd_var = 0.5,
                     sv_mass_s2s_sd_var = 0.5,
                     sv_width_mean = 25,
                     sv_width_p2p_sd_var = 5,
                     sv_width_s2s_sd_var = 5,
                     sv_baseline_mean = 3.0,
                     sv_baseline_s2s_sd_var = 0.5,
                     sv_halflife_mean = 45,
                     sv_halflife_s2s_sd_var = 5,
                     sv_error_var = 0.015,
                     sv_mean_pulse_count = 8,
                     sv_strauss_repulsion = 0,
                     sv_strauss_repulsion_range = 40,
                     pv_sub_mass_mean = 0.5,
                     pv_pop_mass_mean = 0.5,
                     pv_ind_mass_sd = 1,
                     pv_sub_mass_sd = 1,
                     pv_sub_width_mean = 21,
                     pv_pop_width_mean = 21,
                     pv_ind_width_sd = 50,
                     pv_sub_width_sd = 25,
                     pv_pop_baseline_mean = 0.015,
                     pv_sub_baseline = 0.015,
                     pv_sub_baseline_sd = 1,
                     pv_pop_halflife_mean = 25,
                     pv_sub_halflife = 25,
                     pv_sub_halflife_sd = 1,
                     pv_pulse_location = 40,
                     pv_sdscale_pulse_mass = 4,
                     pv_sdscale_pulse_width = 4
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
                         halflife_s2s_sd_var     = prior_halflife_s2s_sd_var_max,
                         error_alpha             = prior_error_alpha,
                         error_beta              = prior_error_beta),
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
                                     pulse_location    = pv_pulse_location,
                                     sdscale_pulse_mass =
                                       pv_sdscale_pulse_mass,
                                     sdscale_pulse_width =
                                       pv_sdscale_pulse_width),
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
                                  mean_pulse_count    = sv_mean_pulse_count,
                                  strauss_repulsion   = sv_strauss_repulsion,
                                  strauss_repulsion_range = sv_strauss_repulsion_range)),
           class = "pop_pulse_spec")
      
      return(pop_ps_obj)
}
