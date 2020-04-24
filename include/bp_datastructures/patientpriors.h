#ifndef GUARD_patientpriors_h
#define GUARD_patientpriors_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/utils.h>

//
// patientpriors.h
//   defining the patientpriors object -- holds priors at patient and pulse
//   level
//
// Authors: Matt Mulvahill
//          Max McGrath
// Created: 12/30/17
// Notes:
//

using namespace Rcpp;


//
// PatientPriors structures
//
//  TODO: abandoned inheritance approach -- clean up comments
//   a. PatientPriors is the base structure containing the components common to
//      all data structures representing the hierarchical level for
//      patient-level priors/population estimates
//
//   1. PopulationEstimates is the structure representing this level for the
//      population model, containing two additional parameters. This object
//      holds parameters that are updated via mcmc.
//
//   2. PatientPriors_Single is the structure representing this level in the
//      single-subject model.  This object contains several additional
//      parameters that are otherwise set at the population level in the pop
//      model. These are all user-set parameters to the prior distributions.
//

// a. PatientPriors structure
// TODO: Will need constructor for Xsingle-subject, Xpopulation, and possibly
// separate ones for trigger and response hormones (so up to 4)

struct PatientPriors {

  //
  // Used in all models
  //
  double baseline_mean; //*
  double baseline_variance; //*
  double halflife_mean; //*
  double halflife_variance; //*
  double mass_mean; //*
  double mass_variance; //*
  double width_mean; //*
  double width_variance; //*
  double error_alpha;
  double error_beta;
  PulseUtils pu;

  //
  // For Population models:
  //

  // Population uniform prior maximums (likely to be altered)
  double mass_p2p_sd_var_max;
  double mass_s2s_sd_var_max;
  double width_p2p_sd_var_max;
  double width_s2s_sd_var_max;
  double baseline_s2s_sd_var_max; //*
  double halflife_s2s_sd_var_max; //*
  double error_mean_pulse_count;

  //
  // For Single-subject model:
  //

  // Member variables not used in Population model:
  double mass_sd_param;  //parameter in the Cauchy prior on the pulse-to-pulse sd mass
  double width_sd_param; //parameter in the Cauchy prior on the pulse-to-pulse sd width
  double num_orderstat;
  int    pulse_count;             // prior number of pulses, i.e. strauss_rate/beta
  double strauss_repulsion;       // strauss gamma for secondary/non-hc interaction
  double strauss_repulsion_range; // range of secondary/non-hardcore interaction

  //
  // Population model constructor:
  //
  PatientPriors(double prior_mass_mean,
                double prior_mass_var,                
                double prior_mass_p2p_sd_var_max,     
                double prior_mass_s2s_sd_var_max,     
                double prior_width_mean,              
                double prior_width_var,               
                double prior_width_p2p_sd_var_max,    
                double prior_width_s2s_sd_var_max,
                double prior_baseline_mean,
                double prior_baseline_var,
                double prior_baseline_s2s_sd_var_max, 
                double prior_halflife_mean,           
                double prior_halflife_var,            
                double prior_halflife_s2s_sd_var_max, 
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
    mass_p2p_sd_var_max = prior_mass_s2s_sd_var_max;
    mass_s2s_sd_var_max = prior_mass_s2s_sd_var_max;
    width_p2p_sd_var_max = prior_width_p2p_sd_var_max;
    width_s2s_sd_var_max = prior_width_s2s_sd_var_max;
    baseline_s2s_sd_var_max = prior_baseline_s2s_sd_var_max;
    halflife_s2s_sd_var_max = prior_halflife_s2s_sd_var_max;

    // Other population-only variables
    error_mean_pulse_count = prior_error_mean_pulse_count;

    // Set single-subject only variables to 0
    mass_sd_param           = 0;
    width_sd_param          = 0;
    num_orderstat           = 0;
    pulse_count             = 0;
    strauss_repulsion       = 0;
    strauss_repulsion_range = 0;

  };

  // Null constructor for population, patient object
  PatientPriors() { };


  // Single-subject constructor (Base + 8 parameters)
  PatientPriors(double prior_baseline_mean,
                double prior_baseline_variance,
                double prior_halflife_mean,
                double prior_halflife_variance,
                double prior_mass_mean,
                double prior_mass_variance,
                double prior_width_mean,
                double prior_width_variance,
                double prior_mass_sd_param,
                double prior_width_sd_param,
                double prior_error_alpha,
                double prior_error_beta,
                int    prior_pulse_count,  // does this prior really need to be an int? Or can we say 12.2 for prior pulse count?
                double prior_strauss_repulsion,
                //double prior_strauss_hardcore_range, // not in single-subj
                double prior_strauss_repulsion_range) {

    // All models
    baseline_mean     = prior_baseline_mean;
    baseline_variance = prior_baseline_variance;
    halflife_mean     = prior_halflife_mean;
    halflife_variance = prior_halflife_variance;
    mass_mean         = prior_mass_mean;
    mass_variance     = prior_mass_variance;
    width_mean        = prior_width_mean;
    width_variance    = prior_width_variance;

    // Set population model-only parms to 0
    mass_p2p_sd_var_max     = 0;
    mass_s2s_sd_var_max     = 0;
    width_p2p_sd_var_max    = 0;
    width_s2s_sd_var_max    = 0;
    baseline_s2s_sd_var_max = 0;
    halflife_s2s_sd_var_max = 0;
    error_mean_pulse_count = 0;

    // Single-subject only
    mass_sd_param           = prior_mass_sd_param;
    width_sd_param          = prior_width_sd_param;
    error_alpha             = prior_error_alpha;
    error_beta              = 1 / prior_error_beta;
    num_orderstat           = pu.orderstat_default();
    pulse_count             = prior_pulse_count;
    strauss_repulsion       = prior_strauss_repulsion; //gamma
    strauss_repulsion_range = prior_strauss_repulsion_range;

  }


// A more appropriate name for the population version
typedef struct PatientPriors PopulationEstimates;

};

#endif
