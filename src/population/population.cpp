#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <bp_population/bp_population.h>
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
Rcpp::List population_(Rcpp::NumericMatrix concentrations,
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
  // Set to print every nth iteration for verbose screen output  
  int verbose_iter = 5000;

  // Check for valid input
  if ( !inpriors.inherits("bp_priors") ) stop("priors argument must be a bp_priors object");
  if ( !proposalvars.inherits("bp_proposalvariance") ) stop("proposalvars argument must be a bp_proposalvariance object");
  if ( !startingvals.inherits("bp_startingvals") ) stop("startingvals argument must be a bp_startingvals object");
  if ( concentrations.nrow() != time.size() ) stop("Time and concentration vectors must be the same size");
         if ( location_prior.size() != 1 ) stop("Location prior vector must be length 1");
  
  // Cleanup variable names 
  int adj_iter       = pv_adjust_iter;
  int adj_max        = pv_adjust_max_iter;
  double biv_target  = bivariate_pv_target_ratio;
  double univ_target = univariate_pv_target_ratio;
  std::string loc_prior = Rcpp::as<std::string>(location_prior);
    

  Rcpp::Rcout << "Location prior is: " << loc_prior << "\n\n";
  
  // Create priors object
    // Q NEC 5/3/20: Does the order matter?
  PopulationPriors priors(inpriors["mass_mean"],
                          inpriors["mass_var"],
                          inpriors["mass_p2p_sd_param"],
                          inpriros["mass_s2s_sd_param"],
                          inpriors["width_mean"],              
                          inpriors["width_var"],               
                          inpriors["width_p2p_sd_param"],
                          inpriors["width_s2s_sd_param"],
                          inpriors["baseline_mean"],
                          inpriors["baseline_var"],
                          inpriors["baseline_sd_param"],
                          inpriors["halflife_mean"],           
                          inpriors["halflife_var"],            
                          inpriors["halflife_sd_param"],
                          inpriors["error_alpha"],             
                          inpriors["error_beta"],
                          inpriors["pulse_count"],
                          inpriors["strauss_repulsion"],
                          inpriors["strauss_repulsion_range"]);

  std::cout << "Priors created\n";

  // Create patient priors object
  PatientPriors patientPriors(startingvals["mass_mean"],
                              startingvals["mass_p2p_sd_var"],
                              startingvals["mass_s2s_sd_var"],
                              startingvals["width_mean"],
                              startingvals["width_p2p_sd_var"],
                              startingvals["width_s2s_sd_var"],
                              startingvals["baseline_mean"],
                              startingvals["baseline_s2s_sd_var"],
                              startingvals["halflife_mean"],
                              startingvals["halflife_s2s_sd_var"], 
                              startingvals["error_alpha"],
                              startingvals["error_beta"],
                              startingvals["error_mean_pulse_count"],
                              startingvals["strauss_repulsion"],
                              startingvals["strauss_repulsion_range"],
                              true);
 
  // Create population estimates object
  PopulationEstimates estimates(startingvals["mass_mean"],
                                startingvals["mass_sd"],
                                startingvals["mass_sd_s2s"],
                                startingvals["width_mean"],
                                startingvals["width_sd"],
                                startingvals["width_sd_s2s"],
                                startingvals["baseline_mean"],
                                startingvals["baseline_sd"],
                                startingvals["halflife_mean"],
                                startingvals["halflife_sd"],
                                startingvals["error_var"]);

  // Create subject level estimates object
  PatientEstimates patEstimates(startingvals["baseline_mean"],
                                startingvals["halflife_mean"],
                                startingvals["error_var"],
                                startingvals["mass_mean"],
                                startingvals["width_mean"]);

  std::cout << "Patient estimates created\n";


  // Get number of patients (should this be a function arg?)
  int numPats = concentrations.ncol();

  std::cout << "Number of patients: " << numPats << "\n";
  std::cout << "Element 5,5: " << concentrations(5, 5) << "\n";
  
  // Create vector of patients
  std::vector<Patient> pats;

  // Fill vector of patients
  for(int i = 0; i < numPats; i++) {
    PatientData tempdata(time, concentrations.column(i));
    std::cout << "Created data\n";
    Patient temppatient(tempdata, patEstimates);
    std::cout << "Created temp patient\n";
    pats.push_back(temppatient);
    std::cout << "Patient " << i << " created\n";
  }

  std::cout << "Patients created\n";

  // Create population object and pointer to it
  Population pop(pats, priors, patientPriors);
  Population * population = &pop;
    
  std::cout << "Population created\n\n";
  int j = 1; 
  // Output diagnostics for each patient
  for(auto &patient : population->patients) {
 
    std::cout << "Patient " << j << ":\n";
    std::cout << "First conc: " << patient.data.concentration[0] << "\n";
    std::cout << "Num obs: "    << patient.data.number_of_obs << "\n";
    std::cout << "Fit start: "  << patient.data.fitstart << "\n";
    std::cout << "Fit end: "    << patient.data.fitend << "\n";
    std::cout << "Avg period: " << patient.data.avg_period_of_obs << "\n";
    std::cout << "Likelihood: " << patient.likelihood(false) << "\n";
    std::cout << "Mass P2P sd: " << population->patPriors.mass_p2p_sd << "\n";
    std::cout << "Width P2P sd: " << population->patPriors.width_p2p_sd << "\n\n";
    j++;

  }

  //----------------------------------------
  // Sample MMH Objects
  //----------------------------------------
  BirthDeathProcess birth_death;
  
  Pop_DrawSDRandomEffects draw_sd_masses(proposalvars["sub_mass_mean"], adj_iter,
                                         adj_max, univ_target, false,
                                         verbose, verbose_iter);

  Pop_DrawSDRandomEffects draw_sd_width(proposalvars["sub_width_mean"], adj_iter,
                                         adj_max, univ_target, false,
                                         verbose, verbose_iter);

  //----------------------------------------
  // Sample MMH Objects
  //----------------------------------------
  for(int iteration = 0; iteration < 50; iteration++) {

    draw_sd_masses.sample(population, &population->patPriors.mass_p2p_sd, population, iteration);
    draw_sd_width.sample(population, &population->patPriors.width_p2p_sd, population, iteration);

  }


  /*
  // Declare BDMCMC objects
  BirthDeathProcess birth_death;
  Pop_DrawPopMeans draw_mass_mean(proposalvars["pop_mass_mean"], adj_iter,
                                  adj_max, univ_target, false, true, false,
                                  verbose, verbose_iter);
  draw_mass_mean.sample(population, &population->popPriors.mass_mean, iteration);
  
  // Temporary placeholder to specify iteration (will be loop)
  int i = 0;

  // Run birth death algorithm for each patient
  for(int j = 0; j < numPats; j++) {

    // Create pointer to patient
    Patient * temp_patient = &population->patients[j];

    std::cout << "Patient " << j << "\n";
    std::cout << "Pulses before: " << population->patients[j].get_pulsecount() << "\n";

    // BD Sample
    birth_death.sample(temp_patient, false, i);
    std::cout << "Likelihood: " << population->patients[j].likelihood(false) << "\n";

    std::cout << "Pulses after: " << population->patients[j].get_pulsecount() << "\n\n";

  }

  */


  j = 1;
  std::cout << "\n";
  
  // Output diagnostics for each patient
  for(auto &patient : population->patients) {

    std::cout << "Patient " << j << ":\n";
    std::cout << "Likelihood: " << patient.likelihood(false) << "\n"; 
    std::cout << "Mass P2P sd: " << population->patPriors.mass_p2p_sd << "\n";
    std::cout << "Width P2P sd: " << population->patPriors.width_p2p_sd << "\n\n";
    j++;

  }

  Rcpp::List test = List::create(Named("Patient1") = 1,
                                 Named("Patient2") = 2);

  return test;
};
