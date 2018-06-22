#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <bpmod_singlesubject/bpmod_singlesubject.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif

using namespace Rcpp;

//
// singlesubject.cpp
//   Implementation of the single-subject model.  Creates Patient object
//   (containing data, priors, estimates), the samplers (Metropolis Hastings,
//   birth death process), and output routines (chains, verbose
//   output/diagnostics).  Intended to be called by a function in the R package.
//
// Author: Matt Mulvahill
// Notes:
//


// [[Rcpp::export]]
Rcpp::List singlesubject_(Rcpp::NumericVector concentration,
                          Rcpp::NumericVector time,
                          Rcpp::List priors,
                          Rcpp::List proposalvars,
                          Rcpp::List startingvals,
                          int mcmc_iterations,
                          int thin,
                          int burnin,
                          bool verbose,
                          int pv_adjust_iter,
                          int pv_adjust_max_iter,
                          double bivariate_pv_target_ratio,
                          double univariate_pv_target_ratio)
{

  // Check for valid input
  if ( !priors.inherits("bp_priors") ) stop("priors argument must be a bp_priors object");
  if ( !proposalvars.inherits("bp_proposalvariance") ) stop("proposalvars argument must be a bp_proposalvariance object");
  if ( !startingvals.inherits("bp_startingvals") ) stop("startingvals argument must be a bp_startingvals object");
  if ( concentration.size() != time.size() ) stop("Time and concentration vectors must be the same size");

  // Create shorter names just for cleaner code appearance
  int adj_iter       = pv_adjust_iter;
  int adj_max        = pv_adjust_max_iter;
  double biv_target  = bivariate_pv_target_ratio;
  double univ_target = univariate_pv_target_ratio;


  // Create patient data object
  PatientData pdone(time, concentration);

  //Create priors object
  PatientPriors ppsingle(priors["baseline_mean"],
                         priors["baseline_variance"],
                         priors["halflife_mean"],
                         priors["halflife_variance"],
                         priors["mass_mean"],
                         priors["mass_variance"],
                         priors["width_mean"],
                         priors["width_variance"],
                         priors["mass_sdmax"],
                         priors["width_sdmax"],
                         priors["error_alpha"],
                         priors["error_beta"],
                         priors["pulse_count"],
                         priors["strauss_repulsion"],
                         priors["strauss_repulsion_range"]);

  // Create estimates object (w/ starting vals)
  PatientEstimates pesingle(startingvals["baseline"],
                            startingvals["halflife"],
                            startingvals["errorsq"],
                            startingvals["mass_mean"],
                            startingvals["width_mean"],
                            startingvals["mass_sd"],
                            startingvals["width_sd"]);

  // Create pointers
  PatientData * data_obj = &pdone;
  PatientPriors * priors_obj = &ppsingle;
  PatientEstimates * estimates_obj = &pesingle;

  // Now take all of this and create a Patient object
  Patient pat(data_obj, priors_obj, estimates_obj);
  Patient * patient = &pat;


  //----------------------------------------
  // Create sampler objects
  //----------------------------------------

  // Birth-death process
  BirthDeathProcess birth_death;

  // Modified Metropolis Hastings for fixed effects (mean mass & mean width)
  SS_DrawFixedEffects draw_fixeff_mass(proposalvars["mass_mean"], adj_iter,
                                       adj_max, univ_target, false);
  SS_DrawFixedEffects draw_fixeff_width(proposalvars["width_mean"], adj_iter,
                                        adj_max, univ_target, true);

  // Modified Metropolis Hastings for the standard deviation of the random
  // effects (sd mass & sd width) (patient level estimate)
  //SS_DrawSDRandomEffects draw_sd_masses(proposalvars["mass_sd"], adj_iter,
  //                                      adj_max, univ_target, false);
  //SS_DrawSDRandomEffects draw_sd_widths(proposalvars["width_sd"], adj_iter,
  //                                      adj_max, univ_target, true);

  // Bivariate Modified Metropolis Hastings for the baseline and half-life
  arma::vec bhl_pv = { proposalvars["baseline"], proposalvars["halflife"] };
  SS_DrawBaselineHalflife draw_blhl(bhl_pv, adj_iter, adj_max, biv_target);

  // Modified Metropolis Hastings for pulse locations (pulse level)
  //if ( priors["location_prior_type"] == "strauss" ) {
    SS_DrawLocationsStrauss draw_locations(proposalvars["location"], adj_iter,
                                           adj_max, univ_target);
  //} else {
  //  SS_DrawLocationsOS draw_locations(proposalvars["location"], adj_iter,
  //                                    adj_max, univ_target);
  //}
  SS_DrawRandomEffects draw_masses(proposalvars["pulse_mass"], adj_iter,
                                   adj_max, univ_target, false);
  SS_DrawRandomEffects draw_widths(proposalvars["pulse_width"], adj_iter,
                                   adj_max, univ_target, true);
  //SS_DrawTVarScale draw_tvarscale_mass(proposalvars["sdscale_pulse_mass"],
  //                                     adj_iter, adj_max, univ_target, false);
  //SS_DrawTVarScale draw_tvarscale_width(proposalvars["sdscale_pulse_width"],
  //                                      adj_iter, adj_max, univ_target, true);


  // Create output objects (chains)
  Rcpp::Rcout << "mcmc iterations = " << mcmc_iterations << std::endl;
  Rcpp::Rcout << "thin = " << thin << std::endl;
  Rcpp::Rcout << "burnin = " << burnin << std::endl;
  Chains chains(mcmc_iterations, thin, burnin, false);


  //----------------------------------------
  // Sample MMH objects
  //----------------------------------------
  for (int i = 0; i < mcmc_iterations; i++) {

    checkUserInterrupt();
    birth_death.sample(patient, false);
    draw_fixeff_mass.sample(patient, &patient->estimates->mass_mean);
    draw_fixeff_width.sample(patient, &patient->estimates->width_mean);
    //draw_sd_masses.sample(patient, &patient->estimates->mass_sd, patient);
    //draw_sd_widths.sample(patient, &patient->estimates->mass_sd, patient);
    draw_blhl.sample(patient, &patient->estimates->baseline_halflife);
    draw_locations.sample_pulses(patient);
    draw_masses.sample_pulses(patient);
    draw_widths.sample_pulses(patient);
    //draw_tvarscale_mass.sample_pulses(patient);
    //draw_tvarscale_width.sample_pulses(patient);

    //print_diagnostic_output(verbose);
    ////std::cout << "Iteration " << i << " Number of pulses = " << patient->pulses.size() << std::endl;
    chains.save_sample(patient, i);

  }

  return chains.output();

}
