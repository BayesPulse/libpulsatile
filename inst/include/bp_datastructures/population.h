#ifndef GUARD_population_h
#define GUARD_population_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include "bp_datastructures/patient.h"
#include "bp_datastructures/patientdata.h"
#include "bp_datastructures/patientestimates.h"
#include "bp_datastructures/patientpriors.h"
#include "bp_datastructures/pulseestimates.h"
#include "bp_datastructures/populationestimates.h"

//
// population.h
//   defining the population class and subclasses
//


// In context of population model, patient priors are population estimates.
//typedef struct PatientPriors PopulationEstimates;


// The user defined values for the priors on the population level parameters
// ** This structure contains all the variables that the user sets when setting the population priors in the population model **
struct PopulationPriors {

  double baseline_mean;             //Mean of the prior on the pop mean baseline
  double baseline_variance;         //Variance of the prior on the pop mean baseline
  double baseline_sd_param;         //Subj to Subj SD of the patient baselines
    
  double halflife_mean;             //Mean of the prior on the pop mean halflife
  double halflife_variance;         //Variance of the prior on the pop mean halflife
  double halflife_sd_param;         //Subj-to-subj SD of the patient halflives
  

  double mass_mean;                 //Mean of the prior on the pop mean pulse mass
  double mass_variance;             //Var of the prior on the pop mean pulse mass
  double mass_p2p_sd_param;         //pulse to pulse SD of the ind. pulse masses
  double mass_s2s_sd_param;         //subj to subj SD of the subj mean pulse mass
    
  double width_mean;                //Mean of the prior on the pop mean pulse width
  double width_variance;            //Var of the prior on the pop mean pulse width
  double width_p2p_sd_param;        //pulse to pulse SD of the ind pulse widths
  double width_s2s_sd_param;        //subj to subj SD of the subj mean pulse width

  double error_alpha;               // Alpha for gamma prior on model error
  double error_beta;                // Beta for gamma prior on model error
    

  // Constructor
  PopulationPriors(double prior_mass_mean,
                   double prior_mass_var,                
                   double prior_mass_p2p_sd_param,
                   double prior_mass_s2s_sd_param,
                   double prior_width_mean,              
                   double prior_width_var,               
                   double prior_width_p2p_sd_param,
                   double prior_width_s2s_sd_param,
                   double prior_baseline_mean,
                   double prior_baseline_var,
                   double prior_baseline_sd_param,
                   double prior_halflife_mean,           
                   double prior_halflife_var,            
                   double prior_halflife_sd_param,
                   double prior_error_alpha,
                   double prior_error_beta
                   ){

    mass_mean          = prior_mass_mean;
    mass_variance      = prior_mass_var;
    mass_p2p_sd_param  = prior_mass_p2p_sd_param;
    mass_s2s_sd_param  = prior_mass_s2s_sd_param;

    width_mean         = prior_width_mean;
    width_variance     = prior_width_var;
    width_p2p_sd_param = prior_width_p2p_sd_param;
    width_s2s_sd_param = prior_width_s2s_sd_param;
      
    baseline_mean      = prior_baseline_mean;
    baseline_variance  = prior_baseline_var;
    baseline_sd_param  = prior_baseline_sd_param;
      
    halflife_mean      = prior_halflife_mean;
    halflife_variance  = prior_halflife_var;
    halflife_sd_param  = prior_halflife_sd_param;

    error_alpha        = prior_error_alpha;
    error_beta         = 1 / prior_error_beta;

  }
};


struct PopulationEstimates {
  double mass_mean;
  double mass_s2s_sd;
  double mass_p2p_sd;

  double width_mean;
  double width_s2s_sd;
  double width_p2p_sd;

  double baseline_mean;
  double baseline_sd;

  double halflife_mean;
  double halflife_sd;

  PopulationEstimates(double sv_mass_mean,
                      double sv_mass_s2s_sd,
                      double sv_mass_p2p_sd,
                      double sv_width_mean,
                      double sv_width_s2s_sd,
                      double sv_width_p2p_sd,
                      double sv_baseline_mean,
                      double sv_baseline_sd,
                      double sv_halflife_mean,
                      double sv_halflife_sd) {

  mass_mean = sv_mass_mean;
  mass_s2s_sd = sv_mass_s2s_sd;
  mass_p2p_sd = sv_mass_p2p_sd;
  width_mean = sv_width_mean;
  width_s2s_sd = sv_width_s2s_sd;
  width_p2p_sd = sv_width_p2p_sd;
  baseline_mean = sv_baseline_mean;
  baseline_sd = sv_baseline_sd;
  halflife_mean = sv_halflife_mean;
  halflife_sd = sv_halflife_sd;

  }


};


// Parameters in the cox process for the response hormone intensity
struct AssocEstimates {
  double cluster_size;      // cluster size (rho)
  double log_cluster_size;  // log scale cluster size
  double cluster_width;     // cluster width (nu)
  double log_cluster_width;
};



struct Population {

  // Member objects
  std::vector<Patient> patients;
  PopulationPriors priors;
  PopulationEstimates estimates;
  //PatientPriors patPriors;
  //AssocEstimates associations;

  /*  // Constructor (w/ response hormone)
  Population(std::vector<Patient> in_patients,
             PopulationPriors in_priors,
             PopulationEstimates in_estimates,
             AssocEstimates in_associations) {

    patients     = in_patients;
    priors       = in_priors;
    estimates    = in_estimates;
    associations = in_associations;

  } */
  
  // Constructor (w/o response hormone)
  Population(std::vector<Patient> in_patients,
             PopulationPriors in_priors,
             PopulationEstimates in_estimates) : 
         priors(in_priors),
         estimates(in_estimates) {

    patients  = in_patients;
    priors    = in_priors;
    estimates = in_estimates;


  }

  int get_patientcount() { return patients.size(); };

  void matchPatPriorsToPop() {
    for(auto &pat : patients) {
      // For draw_fixeff
      pat.priors.mass_mean = estimates.mass_mean;
      pat.priors.width_mean = estimates.width_mean;
      pat.priors.baseline_mean = estimates.baseline_mean;
      pat.priors.halflife_mean = estimates.halflife_mean;
      pat.priors.mass_variance = pow(estimates.mass_s2s_sd, 2);
      pat.priors.width_variance = pow(estimates.width_s2s_sd, 2);
      pat.priors.baseline_variance = pow(estimates.baseline_sd, 2);
      pat.priors.halflife_variance = pow(estimates.halflife_sd, 2);

      // For draw_randomeff
      pat.estimates.mass_sd = estimates.mass_p2p_sd;
      pat.estimates.width_sd = estimates.width_p2p_sd;
    }
  };

// This function fixes parameter estimates for parameters the user
//   chooses not to estimate 
void fix_estimates(Rcpp::List fix_params,
                     Rcpp::NumericVector pat_mass_mean_vec,
                     Rcpp::NumericVector pat_width_mean_vec,
                     Rcpp::NumericVector pat_baseline_vec,
                     Rcpp::NumericVector pat_halflife_vec,
                     Rcpp::NumericVector pulse_count_vec,
                     Rcpp::NumericVector masses_vec,
                     Rcpp::NumericVector width_vec,
                     Rcpp::NumericVector mass_tvarscale_vec,
                     Rcpp::NumericVector width_tvarscale_vec,
                     Rcpp::NumericVector location_vec) {

  int i = 0;
  int j = 0;
  int l = 0;

  // Fix individual mass means
  if(fix_params["pat_mass_mean"]) {
    i = 0;
    for(auto &pat : patients) {
      pat.estimates.mass_mean = pat_mass_mean_vec(i);
      i++;
    }
  Rcpp::Rcout << "Patient mass means fixed\n";
  }

  // Fix individual width means
  if(fix_params["pat_width_mean"]) {
    i = 0;
    for(auto &pat : patients) {
      pat.estimates.width_mean = pat_width_mean_vec(i);
      i++;
    }
    Rcpp::Rcout << "Patient width means fixed\n";
  
  }
  
  // Fix individual baseline
  if(fix_params["pat_bl_hl"]) {
    i = 0;
    for(auto &pat : patients) {
      pat.estimates.baseline = pat_baseline_vec(i);
      pat.estimates.baseline_halflife(0) = pat_baseline_vec(i);
      pat.estimates.halflife = pat_halflife_vec(i);
      pat.estimates.baseline_halflife(1) = pat_halflife_vec(i);
      i++;
    }
    Rcpp::Rcout << "Patient baseline/halflife fixed\n";
  
  }

  // Fix pulses
  double position, new_mass, new_width, new_tvarscale_mass, new_tvarscale_width,
         new_t_sd_mass, new_t_sd_width;

  if(fix_params["pulse_count"]) {
    i = 0;

    for(auto &pat : patients) {

      // Get number of pulses for patient
      j = pulse_count_vec[i];

      Rcpp::RNGScope rng_scope;

      // Generate first pulse (must be done outside loop, otherwise push_back()
      //  adds first pulse on top of the one generated by pulse constructor)
      position = (fix_params["pulse_location"]) 
        ? location_vec(l) : Rf_runif(pat.data.fitstart, pat.data.fitend);
      new_tvarscale_mass = (fix_params["pulse_mass_sdscale"]) 
        ? mass_tvarscale_vec(l) : Rf_rgamma(2, 0.5);
      new_tvarscale_width = (fix_params["pulse_width_sdscale"]) 
        ? width_tvarscale_vec(l) : Rf_rgamma(2, 0.5);

      // Can fix values to known value (rather than gammas) for testing purposes
      //new_tvarscale_mass = (fix_params["mass_tvarscale"]) ? testMKappaVec(l) : 10;
      //new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : 10;
  
      if(fix_params["pulse_mass"]) {
        new_mass = masses_vec(l);
      } else {
        new_t_sd_mass = pat.estimates.mass_sd / sqrt(new_tvarscale_mass);
        new_mass = -1.0;
        while (new_mass < 0) {
          new_mass = Rf_rnorm(pat.estimates.mass_mean, new_t_sd_mass);
          //new_mass = Rf_rnorm(10, .1);
        }
      }

      if(fix_params["pulse_width"]) {
        new_width = width_vec(l);
      } else {
        new_t_sd_width = pat.estimates.width_sd / sqrt(new_tvarscale_width);
        new_width = -1.0;
        while (new_width < 0) {
          new_width = Rf_rnorm(pat.estimates.width_mean, new_t_sd_width);
          //new_width = Rf_rnorm(70, 1);
        }
      }
    
      PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                               new_tvarscale_width, pat.estimates.get_decay(),
                               pat.data.time);

      pat.pulses.front() = new_pulse;

      l++;

      // Fix remainder of pulses using loop
      for(int k = 0; k < j-1; k++) {
        position = (fix_params["pulse_location"]) 
           ? location_vec(l) : Rf_runif(pat.data.fitstart, pat.data.fitend);
         new_tvarscale_mass = (fix_params["pulse_mass_sdscale"]) 
           ? mass_tvarscale_vec(l) : Rf_rgamma(2, 0.5);
         new_tvarscale_width = (fix_params["pulse_width_sdscale"]) 
           ? width_tvarscale_vec(l) : Rf_rgamma(2, 0.5);

         // Can fix values to known value (rather than gammas) for testing purposes
         //new_tvarscale_mass = (fix_params["mass_tvarscale"]) ? testMKappaVec(l) : 10;
         //new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : 10;
      
         if(fix_params["pulse_mass"]) {
           new_mass = masses_vec(l);
         } else {
           new_t_sd_mass = pat.estimates.mass_sd / sqrt(new_tvarscale_mass);
           new_mass = -1.0;
           while (new_mass < 0) {
             new_mass = Rf_rnorm(pat.estimates.mass_mean, new_t_sd_mass);
             //new_mass = Rf_rnorm(10, .1);
           }
         }

         if(fix_params["pulse_width"]) {
           new_width = width_vec(l);
         } else {
           new_t_sd_width = pat.estimates.width_sd / sqrt(new_tvarscale_width);
           new_width = -1.0;
           while (new_width < 0) {
             new_width = Rf_rnorm(pat.estimates.width_mean, new_t_sd_width);
             //new_width = Rf_rnorm(70, 1);
           }
         }
        
         PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                                  new_tvarscale_width, pat.estimates.get_decay(),
                                  pat.data.time);

         pat.pulses.push_back(new_pulse);

         l++;
      }

      i++;

    }

    Rcpp::Rcout << "Pulses Added Manually\n";

    Rcpp::Rcout << "Details:\n";
    j = l = 1;
    for(auto pat : patients) {
      for(auto pulse : pat.pulses) {
        Rcpp::Rcout << std::setprecision(3) << j << " "
                    << l << " "
                    << pulse.mass << " "
                    << pulse.width << " "
                    << pulse.time << " "
                    << pulse.tvarscale_mass << " "
                    << pulse.tvarscale_width << " " << "\n";
        l++;
      }
      l = 1;
      j++;
    }

    Rcpp::Rcout << "\n";

    }

  };

};

#endif
