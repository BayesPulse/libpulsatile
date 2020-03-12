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
    

  Rcpp::Rcout << "Location prior is: " << loc_prior << std::endl;
  
  //Create priors object
  PatientPriors priors(inpriors["mass_mean"],
                       inpriors["mass_var"],                
                       inpriors["mass_p2p_sd_var_max"],     
                       inpriors["mass_s2s_sd_var_max"],     
                       inpriors["width_mean"],              
                       inpriors["width_var"],               
                       inpriors["width_p2p_sd_var_max"],    
                       inpriors["width_s2s_sd_var_max"],
                       inpriors["baseline_mean"],
                       inpriors["baseline_var"],
                       inpriors["baseline_s2s_sd_var_max"], 
                       inpriors["halflife_mean"],           
                       inpriors["halflife_var"],            
                       inpriors["halflife_s2s_sd_var_max"], 
                       inpriors["error_alpha"],             
                       inpriors["error_beta"],              
                       inpriors["error_mean_pulse_count"]);

  std::cout << "Priors work \n";  

  PatientEstimates estimates(startingvals["mass_mean"],
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

  std::cout << "Estimates work \n";
  
  // Get number of patients (should this be a function arg?)
  int numPat = concentrations.ncol();

  // Create vector of patients
  std::vector<Patient> patients;

  // Fill vector of patients
  for(int i = 0; i < numPat; i++) {
    // Create temp data object for patient
    PatientData tempdata(time, concentrations.column(i));
    std::cout << "Created data.\n";
    Patient temppatient(tempdata, priors, estimates);
    std::cout << "Created patient.\n"; 
    // Add patient to patients vector
    patients.push_back(temppatient);
    std::cout << "Added patient\n.";
  }

  std::cout << "Patients created";

  Rcpp::List test = List::create(Named("test") = 1,
                                       Named("test2") = 2);
  


  return test;
};
