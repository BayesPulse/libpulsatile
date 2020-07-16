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
  double baseline_mean;         //(Pop) mean of the prior on ss baseline
  double halflife_mean;         //(Pop) mean of the priors on ss half-life
  double mass_mean;             //(Pop) mean of the priors on ss pulse mass
  double width_mean;            // (Pop) mean of the priors on ss pulse width
  double error_alpha;           //Shape Parameter in Gamma prior on model error variance
  double error_beta;            //Scale parameter in Gamma prior on model error variance
  double num_orderstat;          //Order statistic input.  Not implemented yet.
  int    pulse_count;            // prior number of pulses, i.e. strauss_rate/beta
  double strauss_repulsion;      // strauss gamma for secondary/non-hc interaction
  double strauss_repulsion_range;// range of secondary/non-hardcore interaction
  PulseUtils pu;

  //
  // For Population models:
  //
  
  // Population uniform prior maximums (likely to be altered)
  double mass_p2p_sd;           //Subj to subj SD in the mean pulse mass
  double mass_s2s_sd;           //Pulse to Pulse SD in the ind. pulse masses
  double width_p2p_sd;          //Subj to subj SD in the mean pulse width
  double width_s2s_sd;          //Pulse to Pulse SD in the ind. pulse widths
  double baseline_sd;           //Subj to subj SD in the baseline estimates
  double halflife_sd;           //Subj to subj SD in the half-life estimates
  

  //
  // For Single-subject model:
  //

  // Member variables not used in Population model:
  double mass_variance; //variance of the prior on the mean pulse mass
  double mass_sd_param;  //parameter in the Cauchy prior on the pulse-to-pulse sd mass
  double width_variance; //variance of the prior on the mean pulse width
  double width_sd_param; //parameter in the Cauchy prior on the pulse-to-pulse sd width
  double baseline_variance; //Variance of the priors on the baseline;
  double halflife_variance; //Variance of the priors on the halflife;

  //
  // Population model constructor:
  //
  PatientPriors(double sv_mass_mean,
                double sv_mass_p2p_sd,
                double sv_mass_s2s_sd,
                double sv_width_mean,
                double sv_width_p2p_sd,
                double sv_width_s2s_sd,
                double sv_baseline_mean,
                double sv_baseline_sd,
                double sv_halflife_mean,
                double sv_halflife_sd,
                double prior_error_alpha,             
                double prior_error_beta){

    // All models
    mass_mean         = sv_mass_mean;
    mass_p2p_sd       = sv_mass_p2p_sd;
    mass_s2s_sd       = sv_mass_s2s_sd;
    width_mean        = sv_width_mean;
    width_p2p_sd      = sv_width_p2p_sd;
    width_s2s_sd      = sv_width_s2s_sd;
    baseline_mean     = sv_baseline_mean;
    baseline_sd       = sv_baseline_sd;
    halflife_mean     = sv_halflife_mean;
    halflife_sd      = sv_halflife_sd;
    error_alpha       = prior_error_alpha;
    error_beta        = 1 / prior_error_beta;

    // Set single-subject only variables to 0
    mass_sd_param           = 0;
    width_sd_param          = 0;
    mass_variance           = 0;
    width_variance          = 0;
    baseline_variance       = 0;
    halflife_variance       = 0;

  };

  //
  // Second population model constructor (holds data for individual patients to be passed to MH fns)
  //
  PatientPriors(int sv_pulse_count,
                double sv_strauss_repulsion,
                double sv_strauss_repulsion_range) {

    pulse_count       = sv_pulse_count;
    strauss_repulsion = sv_strauss_repulsion;
    strauss_repulsion_range = sv_strauss_repulsion_range;
    num_orderstat     = pu.orderstat_default();
  }

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
    mass_p2p_sd     = 0;
    mass_s2s_sd     = 0;
    width_p2p_sd    = 0;
    width_s2s_sd    = 0;
    
    // Single-subject only
    mass_sd_param           = prior_mass_sd_param;
    width_sd_param          = prior_width_sd_param;
    error_alpha             = prior_error_alpha;
    error_beta              = 1 / prior_error_beta;
    num_orderstat           = pu.orderstat_default();
    pulse_count             = prior_pulse_count;
    strauss_repulsion       = prior_strauss_repulsion; //gamma
    strauss_repulsion_range = prior_strauss_repulsion_range;

  };

};
/*
// A more appropriate name for the population version
typedef struct PatientPriors PopulationEstimates;

};
*/

#endif
