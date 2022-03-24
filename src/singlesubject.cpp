#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <bpmod_singlesubject/bpmod_singlesubject.h>
#include <iostream>
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
// Authors: Matt Mulvahill and Max McGrath
// Notes:
//


// [[Rcpp::export]]
Rcpp::List singlesubject_(Rcpp::NumericVector concentration,
                          Rcpp::NumericVector time,
                          Rcpp::CharacterVector location_prior,
                          Rcpp::List inpriors,
                          Rcpp::List proposalvars,
                          Rcpp::List startingvals,
                          int mcmc_iterations,
                          int thin,
                          int burnin,
                          bool verbose,
                          int verbose_iter,
                          int pv_adjust_iter,
                          int pv_adjust_max_iter,
                          double bivariate_pv_target_ratio,
                          double univariate_pv_target_ratio,
                          Rcpp::List fix_params,
                          Rcpp::NumericVector masses_vec,
                          Rcpp::NumericVector width_vec,
                          Rcpp::NumericVector mass_tvarscale_vec,
                          Rcpp::NumericVector width_tvarscale_vec,
                          Rcpp::NumericVector location_vec)
{

  // every nth iteration for printing verbose screen output
  //int verbose_iter = 5000;
 
  // Check for valid input
  if ( !inpriors.inherits("bp_priors") ) stop("priors argument must be a bp_priors object");
  if ( !proposalvars.inherits("bp_proposalvariance") ) stop("proposalvars argument must be a bp_proposalvariance object");
  if ( !startingvals.inherits("bp_startingvals") ) stop("startingvals argument must be a bp_startingvals object");
  if ( concentration.size() != time.size() ) stop("Time and concentration vectors must be the same size");
  if ( location_prior.size() != 1 ) stop("Location prior vector must be length 1");

  // Create shorter names just for cleaner code appearance
  int adj_iter       = pv_adjust_iter;
  int adj_max        = pv_adjust_max_iter;
  double biv_target  = bivariate_pv_target_ratio;
  double univ_target = univariate_pv_target_ratio;
  std::string loc_prior = Rcpp::as<std::string>(location_prior);

  // Create patient data object
  PatientData data(time, concentration);

  Rcpp::Rcout << "Location prior is: " << loc_prior << std::endl;

  //Create priors object
  PatientPriors priors(inpriors["baseline_mean"],
                       inpriors["baseline_variance"],
                       inpriors["halflife_mean"],
                       inpriors["halflife_variance"],
                       inpriors["mass_mean"],
                       inpriors["mass_variance"],
                       inpriors["width_mean"],
                       inpriors["width_variance"],
                       inpriors["mass_prec_param"],
                       inpriors["mass_prec_param_rate"],
                       inpriors["width_prec_param"],
                       inpriors["width_prec_param_rate"],
                       inpriors["error_alpha"],
                       inpriors["error_beta"],
                       inpriors["pulse_count"],
                       inpriors["strauss_repulsion"],
                       inpriors["strauss_repulsion_range"]);

  // Create estimates object (w/ starting vals)
  PatientEstimates estimates(startingvals["baseline"],
                             startingvals["halflife"],
                             startingvals["errorsq"],
                             startingvals["mass_mean"],
                             startingvals["width_mean"],
                             startingvals["mass_prec"],
                             startingvals["width_prec"]);

  // TODO: SHOULD WE USE THE INTERNAL TEST DATASET/PATIENT OBJ?
//  bool TEST = false;
  Patient * patient;

//  if (TEST) {
//    // NOTE: For debugging, add pulses w/ true parms and data
//    DataStructuresUtils utils;
//    Patient pat = utils.create_new_test_patient_obj();
//    patient = &pat;
//    patient = utils.add_default_pulses(patient);
//
//  } else {
    // Now take all of this and create a Patient object
    Patient pat(data, priors, estimates);
    patient = &pat;
//  }

  //double like = patient->likelihood(false);
  //Rcpp::Rcout << "initial likelihood is: " << like << std::endl;
  //Rcpp::Rcout << "Mean conc:\n" << patient->mean_concentration(false) << std::endl;
  //std::cout << "pulse count is: " << patient->get_pulsecount() << std::endl;
  
  

  //----------------------------------------
  // Create sampler objects
  //----------------------------------------

  // Birth-death process
  BirthDeathProcess birth_death;

  // Modified Metropolis Hastings for fixed effects (mean mass & mean width)
  SS_DrawFixedEffects draw_fixeff_mass(proposalvars["mass_mean"], adj_iter,
                                       adj_max, univ_target, false, verbose,
                                       verbose_iter);
  SS_DrawFixedEffects draw_fixeff_width(proposalvars["width_mean"], adj_iter,
                                        adj_max, univ_target, true, verbose,
                                        verbose_iter);

  // Modified Metropolis Hastings for the standard deviation of the random
  // effects (sd mass & sd width) (patient level estimate)
  SS_DrawPrecRandomEffects draw_prec_masses(proposalvars["mass_prec"], adj_iter,
                                        adj_max, univ_target, false, verbose,
                                        verbose_iter);
  SS_DrawPrecRandomEffects draw_prec_widths(proposalvars["width_prec"], adj_iter,
                                       adj_max, univ_target, true, verbose,
                                    verbose_iter);

  // Bivariate Modified Metropolis Hastings for the baseline and half-life
  arma::vec bhl_pv = { proposalvars["baseline"], proposalvars["halflife"] };
  SS_DrawBaselineHalflife draw_blhl(bhl_pv, adj_iter, adj_max, biv_target,
                                    verbose, verbose_iter);


  // Modified Metropolis Hastings for pulse locations (pulse level)
  SS_DrawLocations * draw_locations;
  if ( loc_prior =="strauss" ) {
    //Rcpp::Rcout << "USING STRAUSS LOCATION MH" << std::endl;
    draw_locations = new SS_DrawLocationsStrauss(proposalvars["location"],
                                                 adj_iter, adj_max, univ_target,
                                                 verbose, verbose_iter);
  } else {
    Rcpp::stop("Order statistic prior is not yet implemented in the birth death process");
    //Rcpp::Rcout << "USING ORDERSTAT LOCATION MH" << std::endl;
    //draw_locations = new SS_DrawLocationsOS(proposalvars["location"], adj_iter,
    //                                        adj_max, univ_target, verbose,
    //                                        verbose_iter);
  }

  SS_DrawRandomEffects draw_masses(proposalvars["pulse_mass"], adj_iter,
                                   adj_max, univ_target, false, verbose,
                                   verbose_iter);
  SS_DrawRandomEffects draw_widths(proposalvars["pulse_width"], adj_iter,
                                   adj_max, univ_target, true, verbose,
                                   verbose_iter);
  SS_DrawTVarScale draw_tvarscale_mass(proposalvars["sdscale_pulse_mass"],
                                       adj_iter, adj_max, univ_target, false,
                                       verbose, verbose_iter);
  SS_DrawTVarScale draw_tvarscale_width(proposalvars["sdscale_pulse_width"],
                                        adj_iter, adj_max, univ_target, true,
                                        verbose, verbose_iter);
  SS_DrawError draw_error;


  // Create output objects (chains)
  Rcpp::Rcout << "mcmc iterations = " << mcmc_iterations << std::endl;
  Rcpp::Rcout << "thin = " << thin << std::endl;
  Rcpp::Rcout << "burnin = " << burnin << std::endl;
  Chains chains(mcmc_iterations, thin, burnin, false, verbose, verbose_iter);

  //---------------------------------------
  // Manually fix values for testing
  //---------------------------------------
  patient->fix_estimates(fix_params,
                    masses_vec,
                    width_vec,
                    mass_tvarscale_vec,
                    width_tvarscale_vec,
                    location_vec);

  //----------------------------------------
  // Sample MMH objects
  //----------------------------------------
  for (int iteration = 0; iteration < mcmc_iterations; iteration++) {

    checkUserInterrupt();
    chains.print_diagnostic_output(patient, iteration);

    if (!fix_params["pulse_count"]) birth_death.sample(patient, false, iteration);
    if (!fix_params["mass_mean"]) draw_fixeff_mass.sample(patient, &patient->estimates.mass_mean, iteration);
    if (!fix_params["width_mean"]) draw_fixeff_width.sample(patient, &patient->estimates.width_mean, iteration);
    if (!fix_params["mass_prec"]) draw_prec_masses.sample(patient, &patient->estimates.mass_prec, patient, iteration);
    if (!fix_params["width_prec"]) draw_prec_widths.sample(patient, &patient->estimates.width_prec, patient, iteration);
    if (!fix_params["bl_hl"]) draw_blhl.sample(patient, &patient->estimates.baseline_halflife, iteration);
    if (!fix_params["pulse_location"]) draw_locations->sample_pulses(patient, iteration);
    if (!fix_params["pulse_mass"]) draw_masses.sample_pulses(patient, iteration);
    if (!fix_params["pulse_width"]) draw_widths.sample_pulses(patient, iteration);
    if (!fix_params["pulse_mass_sdscale"]) draw_tvarscale_mass.sample_pulses(patient, iteration);
    if (!fix_params["pulse_width_sdscale"]) draw_tvarscale_width.sample_pulses(patient, iteration);
    if (!fix_params["error_var"]) draw_error.sample(patient);
    chains.save_sample(patient, iteration);

    //arma::vec locations(patient->get_pulsecount());
    //int i = 0;
    //for (auto &pulse : patient->pulses) {
    //  locations[i] = pulse.time;
    //  ++i;
    //}

    //std::cout.precision(11);
    //std::cout.setf(std::ios::fixed);
    //std::cout << "Iteration " << iteration <<
    //  " Number of pulses = " << patient->get_pulsecount() <<
    //  " pulse locations = " <<  locations <<
    //  "; FE Mass = " << patient->estimates.mass_mean <<
    //  "; FE Width = " << patient->estimates.width_mean <<
    //  "; FE Mass SD = " << patient->estimates.mass_sd <<
    //  "; FE Width SD = " << patient->estimates.width_sd <<
    //  //"; Baseline = " << patient->estimates.baseline_halflife(0) <<
    //  //"; Halflife = " << patient->estimates.baseline_halflife(1) <<
    //  std::endl;

  }

  // Any objects created with new must be deleted
  delete draw_locations;

  // Return results object
  return chains.output(patient);

}
