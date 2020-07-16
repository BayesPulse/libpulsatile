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
  double width_s2s_sd_param;        //subj to subj SD of the subj mean pulse widht
    

  // Population only
  double mass_p2p_sd_param;
  double mass_s2s_sd_param;
  double width_p2p_sd_param;
  double width_s2s_sd_param;
  double halflife_sd_param; 
  double baseline_sd_param;  

  double error_alpha;
  double error_beta;

  double pulse_count;

  double num_orderstat;
  double num_pulses;              // prior number of pulses, i.e. strauss_rate/beta
  double strauss_repulsion;       // strauss gamma for secondary/non-hc interaction
  double strauss_hardcore_range;  // range of hardcore interaction
  double strauss_repulsion_range; // range of secondary/non-hardcore interaction

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
                   double prior_halflife_sd_param
                   ){

    mass_mean         = prior_mass_mean;
    mass_variance     = prior_mass_var;
    mass_p2p_sd_param = prior_mass_p2p_sd_param;
    mass_s2s_sd_param = prior_mass_s2s_sd_param;

    width_mean        = prior_width_mean;
    width_variance    = prior_width_var;
    width_p2p_sd_param = prior_width_p2p_sd_param;
    width_s2s_sd_param = prior_width_s2s_sd_param;
      
    baseline_mean     = prior_baseline_mean;
    baseline_variance = prior_baseline_var;
    baseline_sd_param = prior_baseline_sd_param;
      
    halflife_mean     = prior_halflife_mean;
    halflife_variance = prior_halflife_var;
    halflife_sd_param = prior_halflife_sd_param;

    // Population uniform prior maximums (likely to be altered)
      // Comment NEC 5/3/20: These are actually the parameters defining the half-Cauchy distn in the new coding
    mass_p2p_sd_param = prior_mass_s2s_sd;
    mass_s2s_sd_param = prior_mass_s2s_sd;
    width_p2p_sd_param = prior_width_p2p_sd;
    width_s2s_sd_param = prior_width_s2s_sd;
    baseline_sd_param = prior_baseline_sd;
    halflife_sd_param = prior_halflife_sd;

    // Other population-only variables
    //Comment NEC 5/3/20: This is the prior on the "number of pulses" (or rate parameter) in the Strauss process
    // Need to change the name to match population.cpp but I am not sure if this code is needed then.
    error_mean_pulse_count = prior_error_mean_pulse_count;

    // Set single-subject only variables to 0
    num_orderstat           = 0;
    num_pulses              = 0;
    strauss_repulsion       = 0;
    strauss_repulsion_range = 0;

  }
};

/*
//Max: Not sure what any of the below is.  Is it necessary or hold over from previous ideas of Matt's?
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
  double mass_p2p_sd;
  double mass_s2s_sd;
  double width_p2p_sd;
  double width_s2s_sd;
  double baseline_mean;
  double baseline_s2s_sd;
  double halflife_mean;
  double halflife_s2s_sd;
  double error_var;

  PopulationEstimates(double sv_mass_mean,
                      double sv_mass_p2p_sd,
                      double sv_mass_s2s_sd,
                      double sv_width_mean,
                      double sv_width_p2p_sd,
                      double sv_width_s2s_sd,
                      double sv_baseline_mean,
                      double sv_baseline_s2s_sd,
                      double sv_halflife_mean,
                      double sv_halflife_s2s_sd,
                      double sv_error_var) {

     
    baseline_halflife = {sv_baseline_mean, sv_halflife_mean};

    mass_mean       = sv_mass_mean;
    mass_p2p_sd     = sv_mass_p2p_sd;
    mass_s2s_sd     = sv_mass_s2s_sd;
    width_mean      = sv_width_mean;
    width_p2p_sd    = sv_width_p2p_sd;
    width_s2s_sd    = sv_width_s2s_sd;
    baseline_s2s_sd = sv_baseline_s2s_sd;
    halflife_s2s_sd = sv_halflife_s2s_sd;
    error_var       = sv_error_var;

    //pulse_count = 1;

  }
};*/



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
  PopulationPriors popPriors;
  PatientPriors patPriors;
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
             PatientPriors in_estimates) : 
         popPriors(in_priors),
         patPriors(in_estimates) {

    patients     = in_patients;
    popPriors    = in_priors;
    patPriors    = in_estimates;


  }

  int get_patientcount() { return patients.size(); };

};



#endif

