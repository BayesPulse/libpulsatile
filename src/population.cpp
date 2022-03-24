#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <bp_population/bp_population.h>
//#include <bpmod_singlesubject/bpmod_singlesubject.h>
#include <iostream>
#ifndef NORINSIDE
#include <RInside.h>
#endif

using namespace Rcpp;

// population.cpp
//   Implementation of population model.
//
// Author: Max McGrath
// Notes:
//

// [[Rcpp::export]]
Rcpp::List population_(Rcpp::List concentrations,
                       Rcpp::List times,
                       Rcpp::CharacterVector location_prior,
                       Rcpp::List inpriors,
                       Rcpp::List proposalvars,
                       Rcpp::List startingvals,
                       int mcmc_iterations,
                       int thin,
                       int burnin,
                       bool verbose,
                       bool verbose_patient,
                       int verbose_iter,
                       int pv_adjust_iter,
                       int pv_adjust_max_iter,
                       double bivariate_pv_target_ratio,
                       double univariate_pv_target_ratio,
                       Rcpp::List fix_params,
                       Rcpp::NumericVector pat_mass_mean_vec,
                       Rcpp::NumericVector pat_width_mean_vec,
                       Rcpp::NumericVector pat_baseline_vec,
                       Rcpp::NumericVector pat_halflife_vec,
                       Rcpp::NumericVector pulse_count_vec,
                       Rcpp::NumericVector masses_vec,
                       Rcpp::NumericVector width_vec,
                       Rcpp::NumericVector mass_sdscale_vec,
                       Rcpp::NumericVector width_sdscale_vec,
                       Rcpp::NumericVector location_vec)
{

  // Set to print every nth iteration for verbose screen output  
  //int verbose_iter = 1;

  // Check for valid input
  if ( !inpriors.inherits("bp_priors") ) stop("priors argument must be a bp_priors object");
  if ( !proposalvars.inherits("bp_proposalvariance") ) stop("proposalvars argument must be a bp_proposalvariance object");
  if ( !startingvals.inherits("bp_startingvals") ) stop("startingvals argument must be a bp_startingvals object");
  //if ( concentrations.nrow() != time.size() ) stop("Time and concentration vectors must be the same size");
  if ( location_prior.size() != 1 ) stop("Location prior vector must be length 1");
  
  // Cleanup variable names 
  int adj_iter       = pv_adjust_iter;
  int adj_max        = pv_adjust_max_iter;
  double biv_target  = bivariate_pv_target_ratio;
  double univ_target = univariate_pv_target_ratio;
  std::string loc_prior = Rcpp::as<std::string>(location_prior);

  // Create population priors and estimates objects
  PopulationPriors popPriors(inpriors["mass_mean"],
                             inpriors["mass_var"],                
                             inpriors["mass_p2p_prec"],
                             inpriors["mass_p2p_prec_rate"],
                             inpriors["mass_s2s_sd"],     
                             inpriors["width_mean"],              
                             inpriors["width_var"],               
                             inpriors["width_p2p_prec"],
                             inpriors["width_p2p_prec_rate"],
                             inpriors["width_s2s_sd"],
                             inpriors["baseline_mean"],
                             inpriors["baseline_var"],
                             inpriors["baseline_s2s_sd"], 
                             inpriors["halflife_mean"],           
                             inpriors["halflife_var"],            
                             inpriors["halflife_s2s_sd"],
                             inpriors["error_alpha"],
                             inpriors["error_beta"]);

  PopulationEstimates popEstimates(startingvals["mass_mean"],
                                   startingvals["mass_s2s_sd"],
                                   startingvals["mass_p2p_prec"],
                                   startingvals["width_mean"],
                                   startingvals["width_s2s_sd"],
                                   startingvals["width_p2p_prec"],
                                   startingvals["baseline_mean"],
                                   startingvals["baseline_s2s_sd"],
                                   startingvals["halflife_mean"],
                                   startingvals["halflife_s2s_sd"]);

  PatientPriors patientPriors(startingvals["mass_mean"],
                              startingvals["width_mean"],
                              startingvals["baseline_mean"],
                              startingvals["halflife_mean"],
                              startingvals["mass_s2s_sd"],
                              startingvals["width_s2s_sd"],
                              startingvals["baseline_s2s_sd"],
                              startingvals["halflife_s2s_sd"],
                              inpriors["mean_pulse_count"],
                              inpriors["strauss_repulsion"],
                              inpriors["strauss_repulsion_range"]);

  // Create subject level estimates object
  PatientEstimates patEstimates(startingvals["baseline_mean"],
                                startingvals["halflife_mean"],
                                startingvals["error_var"],
                                startingvals["mass_mean"],
                                startingvals["width_mean"],
                                startingvals["mass_p2p_prec"],
                                startingvals["width_p2p_prec"],
                                false);

  Rcpp::Rcout << "Patient priors and estimates created\n";

  // Get number of patients (should this be a function arg?)
  int numPats = concentrations.length();
  
  // Create vector of patients
  std::vector<Patient> pats;

  // Fill vector of patients
  for(int i = 0; i < numPats; i++) {
    // Do some param validation while we're looping
    //if ( times[i].length() != concentrations[i].length() ) { 
    //  stop("Lengths of times and concentrations must match for respective patients");
    //}
    PatientData tempdata(times[i], concentrations[i]);
    Patient temppatient(tempdata, patientPriors, patEstimates);
    pats.push_back(temppatient);
  }

  int q = 0;

  for(auto pat : pats) {
    Rcpp::Rcout << "Patient " << q << " : start = " << pat.data.fitstart << " end = " << pat.data.fitend << "\n";
    q++;
  }

  // Create population object and pointer to it
  Population pop(pats, popPriors, popEstimates);
  Population * population = &pop;

  // Fix params if requested
  //if (any_sug(fix_params) {
    pop.fix_estimates(fix_params,
                      pat_mass_mean_vec,
                      pat_width_mean_vec,
                      pat_baseline_vec,
                      pat_halflife_vec,
                      pulse_count_vec,
                      masses_vec,
                      width_vec,
                      mass_sdscale_vec,
                      width_sdscale_vec,
                      location_vec);
  //}
    
  //----------------------------------------
  // Construct MMH Objects
  //----------------------------------------
  BirthDeathProcess birth_death;
  
  // Population Level   
  Pop_DrawPrecRandomEffects draw_prec_masses(proposalvars["ind_mass_prec"], adj_iter,
                                         adj_max, univ_target, false,
                                         verbose, verbose_iter);
  Pop_DrawPrecRandomEffects draw_prec_width(proposalvars["ind_width_prec"], adj_iter,
                                         adj_max, univ_target, true,
                                         verbose, verbose_iter);
  Pop_DrawS2S_SD draw_s2s_sd_width(proposalvars["pat_width_sd"], adj_iter, adj_max, univ_target,
                                         true, false, false, verbose, verbose_iter);
  Pop_DrawS2S_SD draw_s2s_sd_mass(proposalvars["pat_mass_sd"], adj_iter, adj_max, univ_target,
                                         false, true, false, verbose, verbose_iter);
  Pop_DrawS2S_SD draw_s2s_sd_baseline(proposalvars["pat_baseline_sd"], adj_iter, adj_max, univ_target,
                                         false, false, true, verbose, verbose_iter);
  Pop_DrawS2S_SD draw_s2s_sd_halflife(proposalvars["pat_halflife_sd"], adj_iter, adj_max, univ_target,
                                         false, false, false, verbose, verbose_iter);
  Pop_DrawPopMeans draw_pop_means_width(proposalvars["pop_width_mean"], adj_iter, adj_max, univ_target,
                                         true, false, false, verbose, verbose_iter);
  Pop_DrawPopMeans draw_pop_means_mass(proposalvars["pop_mass_mean"], adj_iter, adj_max, univ_target,
                                         false, true, false, verbose, verbose_iter);
  Pop_DrawPopMeans draw_pop_means_baseline(proposalvars["pop_baseline_mean"], adj_iter, adj_max, univ_target,
                                         false, false, true, verbose, verbose_iter);
  Pop_DrawPopMeans draw_pop_means_halflife(proposalvars["pop_halflife_mean"], adj_iter, adj_max, univ_target,
                                         false, false, false, verbose, verbose_iter);

  // Subject Level
  SS_DrawFixedEffects draw_fixeff_mass(proposalvars["pat_mass_mean"], adj_iter, adj_max, univ_target,
                                       false, verbose, verbose_iter);
  SS_DrawFixedEffects draw_fixeff_width(proposalvars["pat_width_mean"], adj_iter, adj_max, univ_target,
                                        true, verbose, verbose_iter);

  arma::vec blhl_pv = { proposalvars["pat_baseline"], proposalvars["pat_halflife"] };
  SS_DrawBaselineHalflife draw_blhl(blhl_pv, adj_iter, adj_max, biv_target,
                                    verbose, verbose_iter);
  SS_DrawError draw_error;

  // Pulse Level
  SS_DrawLocations * draw_locations;
  if ( loc_prior == "strauss" ) {
    draw_locations = new SS_DrawLocationsStrauss(proposalvars["pulse_location"], adj_iter, adj_max,
                                                univ_target, verbose, verbose_iter);
  } else {
    Rcpp::stop("Order statistic prior is not yet implemented in the birth death process");
  }

  SS_DrawRandomEffects draw_masses(proposalvars["ind_pulse_mass"], adj_iter, adj_max, univ_target,
                                   false, verbose, verbose_iter);
  SS_DrawRandomEffects draw_widths(proposalvars["ind_pulse_width"], adj_iter, adj_max, univ_target,
                                   true, verbose, verbose_iter);
  SS_DrawTVarScale draw_tvarscale_mass(proposalvars["sdscale_pulse_mass"],
                                       adj_iter, adj_max, univ_target, false,
                                       verbose, verbose_iter);
  SS_DrawTVarScale draw_tvarscale_width(proposalvars["sdscale_pulse_width"],
                                        adj_iter, adj_max, univ_target, true,
                                        verbose, verbose_iter);

  // Initialize object to store chains 
  PopChains chains(mcmc_iterations, thin, burnin, false, verbose, verbose_iter, numPats);

  //----------------------------------------
  // Sample MMH Objects
  //----------------------------------------

  for(int iteration = 0; iteration < mcmc_iterations; iteration++) {

    checkUserInterrupt(); 
    chains.print_diagnostic_output(population, iteration, verbose_patient);

    if (!fix_params["mass_p2p_prec"]) {
      draw_prec_masses.sample(population, &population->estimates.mass_p2p_prec, population, iteration);
    }
    if (!fix_params["width_p2p_prec"]) {
      draw_prec_width.sample(population, &population->estimates.width_p2p_prec, population, iteration);
    }
    if (!fix_params["pop_width_mean"]) {
      draw_pop_means_width.sample(population, &population->estimates.width_mean, iteration);
    }
    if (!fix_params["pop_mass_mean"]) {
      draw_pop_means_mass.sample(population, &population->estimates.mass_mean, iteration);
    }
    if (!fix_params["pop_baseline_mean"]) {
      draw_pop_means_baseline.sample(population, &population->estimates.baseline_mean, iteration);
    }
    if (!fix_params["pop_halflife_mean"]) {
      draw_pop_means_halflife.sample(population, &population->estimates.halflife_mean, iteration);
    }
    if (!fix_params["width_s2s_sd"]) {
      draw_s2s_sd_width.sample(population, &population->estimates.width_s2s_sd, population, iteration);
    }
    if (!fix_params["mass_s2s_sd"]) {
      draw_s2s_sd_mass.sample(population, &population->estimates.mass_s2s_sd, population, iteration);
    }
    if (!fix_params["baseline_s2s_sd"]) {
      draw_s2s_sd_baseline.sample(population, &population->estimates.baseline_sd, population, iteration);
    }
    if (!fix_params["halflife_s2s_sd"]) {
      draw_s2s_sd_halflife.sample(population, &population->estimates.halflife_sd, population, iteration);
    }

    population->matchPatPriorsToPop();

    for(auto &pat : population->patients) {
      Patient * patient = &pat;
      if (!fix_params["pulse_count"]) birth_death.sample(patient, false, iteration);
      if (!fix_params["pat_mass_mean"]) draw_fixeff_mass.sample(patient, &patient->estimates.mass_mean, iteration);
      if (!fix_params["pat_width_mean"]) draw_fixeff_width.sample(patient, &patient->estimates.width_mean, iteration);
      if (!fix_params["pat_bl_hl"]){ 
        draw_blhl.sample(patient, &patient->estimates.baseline_halflife, iteration);
        patient->estimates.matchBLHL();
      }
      if (!fix_params["pulse_location"]) draw_locations->sample_pulses(patient, iteration);
      if (!fix_params["pulse_mass"]) draw_masses.sample_pulses(patient, iteration);
      if (!fix_params["pulse_width"]) draw_widths.sample_pulses(patient, iteration);
      if (!fix_params["pulse_mass_sdscale"]) draw_tvarscale_mass.sample_pulses(patient, iteration);
      if (!fix_params["pulse_width_sdscale"]) draw_tvarscale_width.sample_pulses(patient, iteration);
    }

    if (!fix_params["pat_error"]) draw_error.sample(population);
    chains.save_sample(population, iteration);

  }

  delete draw_locations;
  
  return(chains.output(population));

};
