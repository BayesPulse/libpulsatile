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
// Author: Matt Mulvahill
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
  double baseline_mean;
  double baseline_variance;
  double halflife_mean;
  double halflife_variance;
  double mass_mean;
  double mass_variance;
  double width_mean;
  double width_variance;
  PulseUtils pu;

  //
  // For Population models:
  //

  //   Haven't quite wrapped my head around these (is it sd or variance),
  //   believe they correspond to mass_sd and width_sd in patient estimates for
  //   single-subj model.
  double mass_mean_sd;
  double mass_mean_variance;
  double width_mean_sd;
  double width_mean_variance;

  // Population constructor (Base + 2 parameters)
  PatientPriors(double prior_baseline_mean,
                double prior_baseline_variance,
                double prior_halflife_mean,
                double prior_halflife_variance,
                double prior_mass_mean,
                double prior_mass_variance,
                double prior_width_mean,
                double prior_width_variance,
                double prior_mass_mean_sd,
                double prior_width_mean_sd) {

    // All models
    baseline_mean     = prior_baseline_mean;
    baseline_variance = prior_baseline_variance;
    halflife_mean     = prior_halflife_mean;
    halflife_variance = prior_halflife_variance;
    mass_mean         = prior_mass_mean;
    mass_variance     = prior_mass_variance;
    width_mean        = prior_width_mean;
    width_variance    = prior_width_variance;

    // Population-only variables
    mass_mean_sd  = prior_mass_mean_sd;
    width_mean_sd = prior_width_mean_sd;

    // Set single-subject only variables to 0
    mass_sd_max             = 0;
    width_sd_max            = 0;
    error_alpha             = 0;
    error_beta              = 0;
    num_orderstat           = 0;
    pulse_count             = 0;
    strauss_repulsion       = 0;
    strauss_repulsion_range = 0;

  }

  //
  // For Single-subject model:
  //

  // Member variables not used in Population model:
  double mass_sd_max;
  double width_sd_max;
  double error_alpha;
  double error_beta;
  double num_orderstat;
  int    pulse_count;             // prior number of pulses, i.e. strauss_rate/beta
  double strauss_repulsion;       // strauss gamma for secondary/non-hc interaction
  double strauss_repulsion_range; // range of secondary/non-hardcore interaction

  // Single-subject constructor (Base + 8 parameters)
  PatientPriors(double prior_baseline_mean,
                double prior_baseline_variance,
                double prior_halflife_mean,
                double prior_halflife_variance,
                double prior_mass_mean,
                double prior_mass_variance,
                double prior_width_mean,
                double prior_width_variance,
                double prior_mass_sd_max,
                double prior_width_sd_max,
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
    mass_mean_sd  = 0;
    width_mean_sd = 0;

    // Single-subject only
    mass_sd_max             = prior_mass_sd_max;
    width_sd_max            = prior_width_sd_max;
    error_alpha             = prior_error_alpha;
    error_beta              = prior_error_beta;
    num_orderstat           = pu.orderstat_default();
    pulse_count             = prior_pulse_count;
    strauss_repulsion       = prior_strauss_repulsion; //gamma
    strauss_repulsion_range = prior_strauss_repulsion_range;

  }

};

// A more appropriate name for the population version
typedef struct PatientPriors PopulationEstimates;



#endif
