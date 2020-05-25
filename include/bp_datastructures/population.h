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
    mass_p2p_sd_param = prior_mass_s2s_sd_param;
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

  }
};


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

