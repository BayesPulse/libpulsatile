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

//
// population.h
//   defining the population class and subclasses
//


// In context of population model, patient priors are population estimates.
//typedef struct PatientPriors PopulationEstimates;


// The user defined values for the priors on the population level parameters
// ** This structure contains all the variables that the user sets when setting the priors **
struct PopulationPriors {

  // All models
  double baseline_mean;
  double baseline_variance;
  double halflife_mean;
  double halflife_variance;

  double mass_mean;
  double mass_variance;
  double width_mean;
  double width_variance;

  // Population only
  double mass_p2p_sd_param;
  double mass_s2s_sd_param;
  double width_p2p_sd_param;
  double width_s2s_sd_param;
  double halflife_sd_param; 
  double baseline_sd_param;  

  double error_alpha;
  double error_beta;

  double error_mean_pulse_count;

  double num_orderstat;
  double num_pulses;              // prior number of pulses, i.e. strauss_rate/beta
  double strauss_repulsion;       // strauss gamma for secondary/non-hc interaction
  double strauss_hardcore_range;  // range of hardcore interaction
  double strauss_repulsion_range; // range of secondary/non-hardcore interaction

  // Constructor
  PopulationPriors(double prior_mass_mean,
                   double prior_mass_var,                
                   double prior_mass_p2p_sd,     
                   double prior_mass_s2s_sd,     
                   double prior_width_mean,              
                   double prior_width_var,               
                   double prior_width_p2p_sd,    
                   double prior_width_s2s_sd,
                   double prior_baseline_mean,
                   double prior_baseline_var,
                   double prior_baseline_sd, 
                   double prior_halflife_mean,           
                   double prior_halflife_var,            
                   double prior_halflife_sd, 
                   double prior_error_alpha,             
                   double prior_error_beta,              
                   double prior_error_mean_pulse_count){

    // All models
    mass_mean         = prior_mass_mean;
    mass_variance     = prior_mass_var;
    width_mean        = prior_width_mean;
    width_variance    = prior_width_var;
    baseline_mean     = prior_baseline_mean;
    baseline_variance = prior_baseline_var;
    halflife_mean     = prior_halflife_mean;
    halflife_variance = prior_halflife_var;
    error_alpha       = prior_error_alpha;
    error_beta        = prior_error_beta;

    // Population uniform prior maximums (likely to be altered)
    mass_p2p_sd_param = prior_mass_s2s_sd;
    mass_s2s_sd_param = prior_mass_s2s_sd;
    width_p2p_sd_param = prior_width_p2p_sd;
    width_s2s_sd_param = prior_width_s2s_sd;
    baseline_sd_param = prior_baseline_sd;
    halflife_sd_param = prior_halflife_sd;

    // Other population-only variables
    error_mean_pulse_count = prior_error_mean_pulse_count;

    // Set single-subject only variables to 0
    num_orderstat           = 0;
    num_pulses              = 0;
    strauss_repulsion       = 0;
    strauss_repulsion_range = 0;

  }
};

struct PopulationEstimates {

  // Used in all models
  arma::vec baseline_halflife;
  double mass_mean;
  double width_mean;
  //int    pulse_count; // function of linked list instead?;
  // Always use these functions to get these values.  removed them as separate
  // member variables to ensure the result is always up-to-date
  double get_decay() { return log(2) / baseline_halflife(1); }
  // double get_logerrorsq() { return log(errorsq); }

  // Population model only:
  double mass_p2p_sd_var;
  double mass_s2s_sd_var;
  double width_p2p_sd_var;
  double width_s2s_sd_var;
  double baseline_mean;
  double baseline_s2s_sd_var;
  double halflife_mean;
  double halflife_s2s_sd_var;
  double error_var;

  PopulationEstimates(double sv_mass_mean,
                      double sv_mass_p2p_sd_var,
                      double sv_mass_s2s_sd_var,
                      double sv_width_mean,
                      double sv_width_p2p_sd_var,
                      double sv_width_s2s_sd_var,
                      double sv_baseline_mean,
                      double sv_baseline_s2s_sd_var,
                      double sv_halflife_mean,
                      double sv_halflife_s2s_sd_var,
                      double sv_error_var) {

     
    baseline_halflife = {sv_baseline_mean, sv_halflife_mean};

    mass_mean           = sv_mass_mean;
    mass_p2p_sd_var     = sv_mass_p2p_sd_var;
    mass_s2s_sd_var     = sv_mass_s2s_sd_var;
    width_mean          = sv_width_mean;
    width_p2p_sd_var    = sv_width_p2p_sd_var;
    width_s2s_sd_var    = sv_width_s2s_sd_var;
    baseline_s2s_sd_var = sv_baseline_s2s_sd_var;
    halflife_s2s_sd_var = sv_halflife_s2s_sd_var;
    error_var           = sv_error_var;

    //pulse_count = 1;

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

    patients     = in_patients;
    priors       = in_priors;
    estimates    = in_estimates;

  }

};



#endif

