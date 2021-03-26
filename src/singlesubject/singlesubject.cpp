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
                          double univariate_pv_target_ratio,
                          bool test_birthdeath,
                          bool test_fixeff_mass,
                          bool test_fixeff_width,
                          bool test_sd_mass,
                          bool test_sd_width,
                          bool test_blhl,
                          bool test_error,
                          bool test_locations,
                          bool test_masses,
                          bool test_widths,
                          bool test_tvarscale_mass,
                          bool test_tvarscale_width,
                          Rcpp::NumericVector testMassVec,
                          Rcpp::NumericVector testWidthVec,
                          Rcpp::NumericVector testMKappaVec,
                          Rcpp::NumericVector testWKappaVec,
                          Rcpp::NumericVector testLocVec)
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
                       inpriors["mass_sd_param"],
                       inpriors["width_sd_param"],
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


  //---------------------------------------
  // Add pulses (for testing)
  //---------------------------------------
  double position, new_mass, new_width, new_tvarscale_mass, new_tvarscale_width,
         new_t_sd_mass, new_t_sd_width;
  int l = 0;

  if (!test_birthdeath) {

    int num_pulses = inpriors["pulse_count"];

    // Create first pulse (can't do it in loop because first pulse is initialized differently)

    Rcpp::RNGScope rng_scope;
    position = (!test_locations) ? testLocVec(l) : Rf_runif(patient->data.fitstart, patient->data.fitend);
    new_tvarscale_mass = (!test_locations) ? testMKappaVec(l) : Rf_rgamma(2, 0.5);
    new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : Rf_rgamma(2, 0.5);


    if(!test_masses) {
      new_mass = testMassVec(l);
    } else {
      new_t_sd_mass = patient->estimates.mass_sd / sqrt(new_tvarscale_mass);
      new_mass = -1.0;
      while (new_mass < 0) {
        //new_mass = Rf_rnorm(patient->estimates.mass_mean, new_t_sd_mass);
        new_mass = Rf_rnorm(10, .1);
      }
    }

    if(!test_widths) {
      new_width = testWidthVec(l);
    } else {
      new_t_sd_width = patient->estimates.width_sd / sqrt(new_tvarscale_width);
      new_width = -1.0;
      while (new_width < 0) {
        //new_width = Rf_rnorm(patient->estimates.width_mean, new_t_sd_width);
        new_width = Rf_rnorm(70, 1);
      }
    }
    
    PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                             new_tvarscale_width, patient->estimates.get_decay(),
                             patient->data.time);
    patient->pulses.front() = new_pulse;
    l++;
    
    // Create remaining pulses in loop
    for (int k = 0; k < num_pulses - 1; k++) {
      
        Rcpp::RNGScope rng_scope;
        position = (!test_locations) ? testLocVec(l) : Rf_runif(patient->data.fitstart, patient->data.fitend);
        new_tvarscale_mass = (!test_tvarscale_mass) ? testMKappaVec(l) : Rf_rgamma(2, 0.5);
        new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : Rf_rgamma(2, 0.5);

        if(!test_masses) {
          new_mass = testMassVec(l);
        } else {
          new_t_sd_mass = patient->estimates.mass_sd / sqrt(new_tvarscale_mass);
          new_mass = -1.0;
          while (new_mass < 0) {
            //new_mass = Rf_rnorm(patient->estimates.mass_mean, new_t_sd_mass);
            new_mass = Rf_rnorm(10, .1);
          }
        }

        if(!test_widths) {
          new_width = testWidthVec(l);
        } else {
          new_t_sd_width = patient->estimates.width_sd / sqrt(new_tvarscale_width);
          new_width = -1.0;
          while (new_width < 0) {
            //new_width = Rf_rnorm(patient->estimates.width_mean, new_t_sd_width);
            new_width = Rf_rnorm(70, 1);
          }
        }
        
        PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                                 new_tvarscale_width, patient->estimates.get_decay(),
                                 patient->data.time);

        patient->pulses.push_back(new_pulse);

        l++;

  }

  Rcpp::Rcout << "Pulses Added Manually\n";

  Rcpp::Rcout << "Details:\n";
  l = 1;
  for(auto pulse : patient->pulses) {
    Rcpp::Rcout << std::setprecision(3) << l << ": "
                << pulse.mass << " "
                << pulse.width << " "
                << pulse.time << " "
                << pulse.tvarscale_mass << " "
                << pulse.tvarscale_width << " " << "\n";
    l++;
  }

  Rcpp::Rcout << "\n";
  }

  //----------------------------------------
  // Sample MMH objects
  //----------------------------------------
  for (int iteration = 0; iteration < mcmc_iterations; iteration++) {

    checkUserInterrupt();
    chains.print_diagnostic_output(patient, iteration);

    if (test_birthdeath) birth_death.sample(patient, false, iteration);
    if (test_fixeff_mass) draw_fixeff_mass.sample(patient, &patient->estimates.mass_mean, iteration);
    if (test_fixeff_width) draw_fixeff_width.sample(patient, &patient->estimates.width_mean, iteration);
    if (test_sd_mass) draw_sd_masses.sample(patient, &patient->estimates.mass_sd, patient, iteration);
    if (test_sd_width) draw_sd_widths.sample(patient, &patient->estimates.width_sd, patient, iteration);
    if (test_blhl) draw_blhl.sample(patient, &patient->estimates.baseline_halflife, iteration);
    if (test_locations) draw_locations->sample_pulses(patient, iteration);
    if (test_masses) draw_masses.sample_pulses(patient, iteration);
    if (test_widths) draw_widths.sample_pulses(patient, iteration);
    if (test_tvarscale_mass) draw_tvarscale_mass.sample_pulses(patient, iteration);
    if (test_tvarscale_width) draw_tvarscale_width.sample_pulses(patient, iteration);
    if (test_error) draw_error.sample(patient);
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
