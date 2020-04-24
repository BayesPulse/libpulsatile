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
  PopulationPriors priors(inpriors["mass_mean"],
                          inpriors["mass_var"],                
                          inpriors["mass_p2p_sd_var"],     
                          inpriors["mass_s2s_sd_var"],     
                          inpriors["width_mean"],              
                          inpriors["width_var"],               
                          inpriors["width_p2p_sd_var"],    
                          inpriors["width_s2s_sd_var"],
                          inpriors["baseline_mean"],
                          inpriors["baseline_var"],
                          inpriors["baseline_s2s_sd_var"], 
                          inpriors["halflife_mean"],           
                          inpriors["halflife_var"],            
                          inpriors["halflife_s2s_sd_var"], 
                          inpriors["error_alpha"],             
                          inpriors["error_beta"],              
                          inpriors["error_mean_pulse_count"]);

  std::cout << "Priors created\n";
 
  // Create population estimates object
  PopulationEstimates estimates(startingvals["mass_mean"],
                                startingvals["mass_p2p_sd_var"],
                                startingvals["mass_s2s_sd_var"],
                                startingvals["width_mean"],
                                startingvals["width_p2p_sd_var"],
                                startingvals["width_s2s_sd_var"],
                                startingvals["baseline_mean"],
                                startingvals["baseline_s2s_sd_var"],
                                startingvals["halflife_mean"],
                                startingvals["halflife_s2s_sd_var"],
                                startingvals["error_var"]);

  std::cout << "Population estimates created\n";

  // Create subject level estimates object
  PatientEstimates patEstimates(startingvals["baseline_mean"],
                                startingvals["halflife_mean"],
                                startingvals["error_var"],
                                startingvals["mass_mean"],
                                startingvals["width_mean"]);

  // Get number of patients (should this be a function arg?)
  int numPats = concentrations.ncol();
  
  // Create vector of patients
  std::vector<Patient> pats;

  // Fill vector of patients
  for(int i = 0; i < numPats; i++) {
    PatientData tempdata(time, concentrations.column(i));
    Patient temppatient(tempdata, patEstimates);
    pats.push_back(temppatient);
  }

  std::cout << "Patients created\n";

  // Create population object and pointer to it
  Population pop(pats, priors, estimates);
  Population * population = &pop;
    
  std::cout << "Population created\n\n";

  //std::cout << population->patients[0].data.concentration[0] << "\n\n";

  // Initialize birth death object
  BirthDeathProcess birth_death;

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

  // Output diagnostics for each patient
  for(int j = 0; j < numPats; j++) {

    std::cout << "Patient " << j << "\n";
    std::cout << "First conc: " << population->patients[j].data.concentration[0] << "\n";
    std::cout << "Num obs: " << population->patients[j].data.number_of_obs << "\n";
    std::cout << "Fit start: " << population->patients[j].data.fitstart << "\n";
    std::cout << "Fit end: " << population->patients[j].data.fitend << "\n";
    std::cout << "Avg period: " << population->patients[j].data.avg_period_of_obs << "\n";
    std::cout << "Likelihood: " << population->patients[j].likelihood(false) << "\n\n";

  }

  Rcpp::List test = List::create(Named("Patient1") = 1,
                                 Named("Patient2") = 2);

  return test;
};
