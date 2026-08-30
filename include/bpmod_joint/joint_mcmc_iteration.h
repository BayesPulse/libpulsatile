#ifndef GUARD_bpmod_joint_mcmc_iteration_h
#define GUARD_bpmod_joint_mcmc_iteration_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/patient.h>
#include <bp_datastructures/associationparameters.h>
#include <bpmod_singlesubject/bpmod_singlesubject.h>
#include <bpmod_joint/joint_birthdeath.h>
#include <bpmod_joint/joint_draw_association.h>
#include <bpmod_joint/joint_update_lambda.h>

//
// joint_mcmc_iteration.h
//   Orchestrates a single MCMC iteration for the joint hormone model
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// The joint MCMC iteration proceeds as follows:
//
// 1. Driver hormone updates:
//    a. Birth-death process for driver pulses
//    b. Driver baseline and halflife
//    c. Driver mass_mean and width_mean
//    d. Driver pulse parameters (mass, width, location, tvarscale)
//    e. Driver error variance
//    f. **Update all response pulse lambda values**
//
// 2. Response hormone updates:
//    a. Birth-death process for response pulses (coupling-aware)
//    b. Response baseline and halflife
//    c. Response mass_mean and width_mean
//    d. Response pulse parameters (mass, width, location, tvarscale)
//    e. Response error variance
//
// 3. Association parameter updates:
//    a. ρ (cluster size) - MH sampler
//    b. ν (cluster width) - MH sampler
//    c. **Update all response pulse lambda values**
//

using namespace Rcpp;


//
// JointSamplers struct
//   Container for all MCMC samplers needed for joint model
//
struct JointSamplers {

  // Birth-death samplers
  JointBirthDeathProcess joint_bd;

  // Driver hormone samplers (standard single-subject)
  SS_DrawBaselineHalflife *driver_draw_blhl;
  SS_DrawFixedEffects *driver_draw_fixeff_mass;
  SS_DrawFixedEffects *driver_draw_fixeff_width;
  SS_DrawLocations *driver_draw_locations;
  SS_DrawRandomEffects *driver_draw_masses;
  SS_DrawRandomEffects *driver_draw_widths;
  SS_DrawSDRandomEffects *driver_draw_sd_mass;
  SS_DrawSDRandomEffects *driver_draw_sd_width;
  SS_DrawTVarScale *driver_draw_tvarscale_mass;
  SS_DrawTVarScale *driver_draw_tvarscale_width;
  SS_DrawError driver_draw_error;

  // Response hormone samplers (standard single-subject)
  SS_DrawBaselineHalflife *response_draw_blhl;
  SS_DrawFixedEffects *response_draw_fixeff_mass;
  SS_DrawFixedEffects *response_draw_fixeff_width;
  SS_DrawLocations *response_draw_locations;
  SS_DrawRandomEffects *response_draw_masses;
  SS_DrawRandomEffects *response_draw_widths;
  SS_DrawSDRandomEffects *response_draw_sd_mass;
  SS_DrawSDRandomEffects *response_draw_sd_width;
  SS_DrawTVarScale *response_draw_tvarscale_mass;
  SS_DrawTVarScale *response_draw_tvarscale_width;
  SS_DrawError response_draw_error;

  // Association parameter samplers
  Joint_DrawRho *draw_rho;
  Joint_DrawNu *draw_nu;

  // Constructor
  JointSamplers(Rcpp::List proposalvars,
                int adj_iter,
                int adj_max,
                double biv_target,
                double univ_target,
                bool verbose,
                int verbose_iter,
                std::string loc_prior,
                bool driver_lognormal = false,
                bool response_lognormal = false);

  // Destructor
  ~JointSamplers();
};


// Constructor implementation
inline JointSamplers::JointSamplers(Rcpp::List proposalvars,
                                    int adj_iter,
                                    int adj_max,
                                    double biv_target,
                                    double univ_target,
                                    bool verbose,
                                    int verbose_iter,
                                    std::string loc_prior,
                                    bool driver_lognormal,
                                    bool response_lognormal) {

  // Driver samplers
  arma::vec driver_bhl_pv = { Rcpp::as<double>(proposalvars["driver_baseline"]),
                              Rcpp::as<double>(proposalvars["driver_halflife"]) };
  driver_draw_blhl = new SS_DrawBaselineHalflife(driver_bhl_pv, adj_iter, adj_max,
                                                  biv_target, verbose, verbose_iter);

  // Thread the driver's log-normal flag into its mass/width MEAN full
  // conditionals (the MMH container type is `bool`, so the flag is passed at
  // construction rather than read from the Patient).
  driver_draw_fixeff_mass = new SS_DrawFixedEffects(
    Rcpp::as<double>(proposalvars["driver_mass_mean"]),
    adj_iter, adj_max, univ_target, false, verbose, verbose_iter, driver_lognormal);

  driver_draw_fixeff_width = new SS_DrawFixedEffects(
    Rcpp::as<double>(proposalvars["driver_width_mean"]),
    adj_iter, adj_max, univ_target, true, verbose, verbose_iter, driver_lognormal);

  if (loc_prior == "strauss") {
    driver_draw_locations = new SS_DrawLocationsStrauss(
      Rcpp::as<double>(proposalvars["driver_location"]),
      adj_iter, adj_max, univ_target, verbose, verbose_iter);
  } else {
    Rcpp::stop("Order statistic prior not yet implemented");
  }

  driver_draw_masses = new SS_DrawRandomEffects(
    Rcpp::as<double>(proposalvars["driver_pulse_mass"]),
    adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

  driver_draw_widths = new SS_DrawRandomEffects(
    Rcpp::as<double>(proposalvars["driver_pulse_width"]),
    adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

  driver_draw_sd_mass = new SS_DrawSDRandomEffects(
    Rcpp::as<double>(proposalvars["driver_mass_sd"]),
    adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

  driver_draw_sd_width = new SS_DrawSDRandomEffects(
    Rcpp::as<double>(proposalvars["driver_width_sd"]),
    adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

  driver_draw_tvarscale_mass = new SS_DrawTVarScale(
    Rcpp::as<double>(proposalvars["driver_sdscale_pulse_mass"]),
    adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

  driver_draw_tvarscale_width = new SS_DrawTVarScale(
    Rcpp::as<double>(proposalvars["driver_sdscale_pulse_width"]),
    adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

  // Response samplers
  arma::vec response_bhl_pv = { Rcpp::as<double>(proposalvars["response_baseline"]),
                                Rcpp::as<double>(proposalvars["response_halflife"]) };
  response_draw_blhl = new SS_DrawBaselineHalflife(response_bhl_pv, adj_iter, adj_max,
                                                    biv_target, verbose, verbose_iter);

  // Thread the response's log-normal flag into its mass/width MEAN full
  // conditionals (see driver note above).
  response_draw_fixeff_mass = new SS_DrawFixedEffects(
    Rcpp::as<double>(proposalvars["response_mass_mean"]),
    adj_iter, adj_max, univ_target, false, verbose, verbose_iter, response_lognormal);

  response_draw_fixeff_width = new SS_DrawFixedEffects(
    Rcpp::as<double>(proposalvars["response_width_mean"]),
    adj_iter, adj_max, univ_target, true, verbose, verbose_iter, response_lognormal);

  if (loc_prior == "strauss") {
    response_draw_locations = new SS_DrawLocationsStrauss(
      Rcpp::as<double>(proposalvars["response_location"]),
      adj_iter, adj_max, univ_target, verbose, verbose_iter);
  } else {
    Rcpp::stop("Order statistic prior not yet implemented");
  }

  response_draw_masses = new SS_DrawRandomEffects(
    Rcpp::as<double>(proposalvars["response_pulse_mass"]),
    adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

  response_draw_widths = new SS_DrawRandomEffects(
    Rcpp::as<double>(proposalvars["response_pulse_width"]),
    adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

  response_draw_sd_mass = new SS_DrawSDRandomEffects(
    Rcpp::as<double>(proposalvars["response_mass_sd"]),
    adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

  response_draw_sd_width = new SS_DrawSDRandomEffects(
    Rcpp::as<double>(proposalvars["response_width_sd"]),
    adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

  response_draw_tvarscale_mass = new SS_DrawTVarScale(
    Rcpp::as<double>(proposalvars["response_sdscale_pulse_mass"]),
    adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

  response_draw_tvarscale_width = new SS_DrawTVarScale(
    Rcpp::as<double>(proposalvars["response_sdscale_pulse_width"]),
    adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

  // Association parameter samplers
  draw_rho = new Joint_DrawRho(
    Rcpp::as<double>(proposalvars["rho"]),
    adj_iter, adj_max, univ_target, verbose, verbose_iter);

  draw_nu = new Joint_DrawNu(
    Rcpp::as<double>(proposalvars["nu"]),
    adj_iter, adj_max, univ_target, verbose, verbose_iter);
}


// Destructor implementation
inline JointSamplers::~JointSamplers() {
  delete driver_draw_blhl;
  delete driver_draw_fixeff_mass;
  delete driver_draw_fixeff_width;
  delete driver_draw_locations;
  delete driver_draw_masses;
  delete driver_draw_widths;
  delete driver_draw_sd_mass;
  delete driver_draw_sd_width;
  delete driver_draw_tvarscale_mass;
  delete driver_draw_tvarscale_width;

  delete response_draw_blhl;
  delete response_draw_fixeff_mass;
  delete response_draw_fixeff_width;
  delete response_draw_locations;
  delete response_draw_masses;
  delete response_draw_widths;
  delete response_draw_sd_mass;
  delete response_draw_sd_width;
  delete response_draw_tvarscale_mass;
  delete response_draw_tvarscale_width;

  delete draw_rho;
  delete draw_nu;
}


//
// joint_mcmc_iteration()
//   Execute one complete MCMC iteration for the joint model
//
inline void joint_mcmc_iteration(Patient *driver_patient,
                                 Patient *response_patient,
                                 AssociationEstimates *assoc_est,
                                 AssociationPriors *assoc_priors,
                                 JointSamplers &samplers,
                                 int iteration) {

  //
  // STAGE 1: Driver hormone updates
  //

  // Birth-death for driver pulses
  samplers.joint_bd.sample_driver(driver_patient, iteration);

  // Driver parameters
  samplers.driver_draw_blhl->sample(driver_patient,
                                   &driver_patient->estimates.baseline_halflife,
                                   iteration);
  samplers.driver_draw_fixeff_mass->sample(driver_patient,
                                           &driver_patient->estimates.mass_mean,
                                           iteration);
  samplers.driver_draw_fixeff_width->sample(driver_patient,
                                            &driver_patient->estimates.width_mean,
                                            iteration);

  // Driver pulse-to-pulse SDs of mass and width (patient-level). Previously the
  // joint model never sampled these, leaving mass_sd/width_sd frozen at their
  // starting values for the whole chain.
  samplers.driver_draw_sd_mass->sample(driver_patient,
                                       &driver_patient->estimates.mass_sd,
                                       driver_patient, iteration);
  samplers.driver_draw_sd_width->sample(driver_patient,
                                        &driver_patient->estimates.width_sd,
                                        driver_patient, iteration);

  // Driver pulse-level parameters
  samplers.driver_draw_locations->sample_pulses(driver_patient, iteration);
  samplers.driver_draw_masses->sample_pulses(driver_patient, iteration);
  samplers.driver_draw_widths->sample_pulses(driver_patient, iteration);
  // Skip the t-scale (kappa) draws under Gaussian random effects.
  if (!driver_patient->gaussian_random_effects) {
    samplers.driver_draw_tvarscale_mass->sample_pulses(driver_patient, iteration);
    samplers.driver_draw_tvarscale_width->sample_pulses(driver_patient, iteration);
  }

  // Driver error variance
  samplers.driver_draw_error.sample(driver_patient);

  // **CRITICAL**: Update lambda values for all response pulses after driver changes
  update_lambda_values(driver_patient->pulses, response_patient->pulses, *assoc_est);


  //
  // STAGE 2: Response hormone updates
  //

  // Birth-death for response pulses (coupling-aware)
  samplers.joint_bd.sample_response(driver_patient, response_patient,
                                    *assoc_est, iteration);

  // Response parameters
  samplers.response_draw_blhl->sample(response_patient,
                                     &response_patient->estimates.baseline_halflife,
                                     iteration);
  samplers.response_draw_fixeff_mass->sample(response_patient,
                                             &response_patient->estimates.mass_mean,
                                             iteration);
  samplers.response_draw_fixeff_width->sample(response_patient,
                                              &response_patient->estimates.width_mean,
                                              iteration);

  // Response pulse-to-pulse SDs of mass and width (patient-level).
  samplers.response_draw_sd_mass->sample(response_patient,
                                         &response_patient->estimates.mass_sd,
                                         response_patient, iteration);
  samplers.response_draw_sd_width->sample(response_patient,
                                          &response_patient->estimates.width_sd,
                                          response_patient, iteration);

  // Response pulse-level parameters
  samplers.response_draw_locations->sample_pulses(response_patient, iteration);
  samplers.response_draw_masses->sample_pulses(response_patient, iteration);
  samplers.response_draw_widths->sample_pulses(response_patient, iteration);
  // Skip the t-scale (kappa) draws under Gaussian random effects.
  if (!response_patient->gaussian_random_effects) {
    samplers.response_draw_tvarscale_mass->sample_pulses(response_patient, iteration);
    samplers.response_draw_tvarscale_width->sample_pulses(response_patient, iteration);
  }

  // Update lambda values after response pulse locations may have changed
  // This ensures diagnostic output reflects correct lambda for each pulse position
  update_lambda_values(driver_patient->pulses, response_patient->pulses, *assoc_est);

  // Response error variance
  samplers.response_draw_error.sample(response_patient);


  //
  // STAGE 3: Association parameter updates
  //

  // Create JointData for samplers
  JointData joint_data(driver_patient, response_patient, assoc_priors);

  // Sample ρ and ν
  // NOTE: assoc_est is passed as BOTH sampling_unit's companion and the container.
  // Joint_DrawRho/Joint_DrawNu::posterior_function read the CURRENT rho/nu from the
  // container argument; the 3-arg sample() overload would default-construct an empty
  // AssociationEstimates (rho=1, nu=1, log_rho=0), so the MH acceptance ratio would be
  // computed against those phantom constants instead of the real current values,
  // breaking detailed balance. Passing the real assoc_est as the container fixes this.
  samplers.draw_rho->sample(&joint_data, &assoc_est->rho, assoc_est, iteration);
  samplers.draw_nu->sample(&joint_data, &assoc_est->nu, assoc_est, iteration);

  // Update log-scale values to maintain consistency
  assoc_est->set_rho(assoc_est->rho);
  assoc_est->set_nu(assoc_est->nu);

  // **CRITICAL**: Update lambda values after association parameter changes
  update_lambda_values(driver_patient->pulses, response_patient->pulses, *assoc_est);
}


#endif
