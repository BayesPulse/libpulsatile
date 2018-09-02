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
// Author: Matt Mulvahill
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
                          int pv_adjust_iter,
                          int pv_adjust_max_iter,
                          double bivariate_pv_target_ratio,
                          double univariate_pv_target_ratio)
{

  // every nth iteration for printing verbose screen output
  int verbose_iter = 5000;

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
                       inpriors["mass_sdmax"],
                       inpriors["width_sdmax"],
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
                             startingvals["mass_sd"],
                             startingvals["width_sd"]);

  // Now take all of this and create a Patient object
  Patient pat(data, priors, estimates);
  Patient * patient = &pat;

  //
  // ****************************************
  // TODO: NOTE: For debugging, add 11 pulses w/ true parms
  // ****************************************
  //
  DataStructuresUtils utils;
  patient = utils.add_default_pulses(patient);

  //std::cout << "pulse count is: " << patient->get_pulsecount() << std::endl;

  //----------------------------------------
  // Create sampler objects
  //----------------------------------------

  // Birth-death process
  //BirthDeathProcess birth_death;

  // Modified Metropolis Hastings for fixed effects (mean mass & mean width)
  SS_DrawFixedEffects draw_fixeff_mass(proposalvars["mass_mean"], adj_iter,
                                       adj_max, univ_target, false, verbose,
                                       verbose_iter);
  SS_DrawFixedEffects draw_fixeff_width(proposalvars["width_mean"], adj_iter,
                                        adj_max, univ_target, true, verbose,
                                        verbose_iter);

  // Modified Metropolis Hastings for the standard deviation of the random
  // effects (sd mass & sd width) (patient level estimate)
  SS_DrawSDRandomEffects draw_sd_masses(proposalvars["mass_sd"], adj_iter,
                                        adj_max, univ_target, false, verbose,
                                        verbose_iter);
  SS_DrawSDRandomEffects draw_sd_widths(proposalvars["width_sd"], adj_iter,
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
    //Rcpp::Rcout << "USING ORDERSTAT LOCATION MH" << std::endl;
    draw_locations = new SS_DrawLocationsOS(proposalvars["location"], adj_iter,
                                            adj_max, univ_target, verbose,
                                            verbose_iter);
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


  //----------------------------------------
  // Sample MMH objects
  //----------------------------------------
  for (int iteration = 0; iteration < mcmc_iterations; iteration++) {

    checkUserInterrupt();
    chains.print_diagnostic_output(patient, iteration);

    //birth_death.sample(patient, false, iteration);
    //draw_fixeff_mass.sample(patient, &patient->estimates.mass_mean, iteration);
    //draw_fixeff_width.sample(patient, &patient->estimates.width_mean, iteration);
    //draw_sd_masses.sample(patient, &patient->estimates.mass_sd, patient, iteration);
    //draw_sd_widths.sample(patient, &patient->estimates.width_sd, patient, iteration);
    //draw_blhl.sample(patient, &patient->estimates.baseline_halflife, iteration);
    draw_locations.sample_pulses(patient, iteration);
    //draw_masses.sample_pulses(patient, iteration);
    //draw_widths.sample_pulses(patient, iteration);
    draw_tvarscale_mass.sample_pulses(patient, iteration);
    draw_tvarscale_width.sample_pulses(patient, iteration);
    //draw_error.sample(patient);
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

  delete draw_locations;

  return chains.output();

}
