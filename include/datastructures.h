#ifndef GUARD_datastructures_h
#define GUARD_datastructures_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include "utils.h"

//
// datastructures_patient.h
//   defining the datastructures used at the patient and pulse levels
//
// Author: Matt Mulvahill
// Created: 12/30/17
// Notes:
//

using namespace Rcpp;


//
// PatientPriors structures
//
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
struct PatientPriors {

  // Used in all models
  double baseline_mean;
  double baseline_variance;
  double halflife_mean;
  double halflife_variance;
  double mass_mean;
  double mass_variance;
  double width_mean;
  double width_variance;

  // Base constructor
  PatientPriors(double prior_baseline_mean,
                double prior_baseline_variance,
                double prior_halflife_mean,
                double prior_halflife_variance,
                double prior_mass_mean,
                double prior_mass_variance,
                double prior_width_mean,
                double prior_width_variance) {

    baseline_mean     = prior_baseline_mean;
    baseline_variance = prior_baseline_variance;
    halflife_mean     = prior_halflife_mean;
    halflife_variance = prior_halflife_variance;
    mass_mean         = prior_mass_mean;
    mass_variance     = prior_mass_variance;
    width_mean        = prior_width_mean;
    width_variance    = prior_width_variance;

  }

};


// 1. PopulationEstimates struct
struct PopulationEstimates : PatientPriors {
  // Only used in Population models:
  // Haven't quite wrapped my head around these (is it sd or variance), beleive
  // they correspond to mass_sd and width_sd in patient estimates for
  // single-subj model.
  double mass_mean_sd;
  double mass_mean_variance;
  double width_mean_sd;
  double width_mean_variance;

  // Population constructor (Base + 2 parameters)
  PopulationEstimates(double prior_baseline_mean,
                      double prior_baseline_variance,
                      double prior_halflife_mean,
                      double prior_halflife_variance,
                      double prior_mass_mean,
                      double prior_mass_variance,
                      double prior_width_mean,
                      double prior_width_variance,
                      double prior_mass_mean_sd,
                      double prior_width_mean_sd
                     ) :
    PatientPriors(prior_baseline_mean, prior_baseline_variance,
                  prior_halflife_mean, prior_halflife_variance, prior_mass_mean,
                  prior_mass_variance, prior_width_mean, prior_width_variance) {

    mass_mean_sd  = prior_mass_mean_sd;
    width_mean_sd = prior_width_mean_sd;

  }

};


// 2. PatientPriors_Single struct
struct PatientPriors_Single : PatientPriors {

  // Member variables not used in Population model:
  double mass_sd_max;
  double width_sd_max;
  double error_alpha;
  double error_beta;
  double num_orderstat;
  int    pulse_count;             // prior number of pulses, i.e. strauss_rate/beta
  double strauss_repulsion;       // strauss gamma for secondary/non-hc interaction
  double strauss_repulsion_range; // range of secondary/non-hardcore interaction
  //double strauss_hardcore_range;  // range of hardcore interaction (only need
  //                                //one for single subject model and doesn't
  //                                //exist in this struct in pop model)

  // TODO: Will need constructor for Xsingle-subject, Xpopulation, and possibly
  // separate ones for trigger and response hormones (so up to 4)

  // Single-subject constructor (Base + 8 parameters)
  PatientPriors_Single(double prior_baseline_mean,
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
                       int    prior_pulse_count,
                       double prior_strauss_repulsion,
                       //double prior_strauss_hardcore_range, // not in single-subj
                       double prior_strauss_repulsion_range
                      ) :
    PatientPriors(prior_baseline_mean, prior_baseline_variance,
                  prior_halflife_mean, prior_halflife_variance, prior_mass_mean,
                  prior_mass_variance, prior_width_mean, prior_width_variance) {

      mass_sd_max  = prior_mass_sd_max;
      width_sd_max = prior_width_sd_max;
      error_alpha  = prior_error_alpha;
      error_beta   = prior_error_beta;

      num_orderstat = pulseutils::orderstat_default();

      pulse_count             = prior_pulse_count;
      strauss_repulsion       = prior_strauss_repulsion;
      //strauss_hardcore_range  = prior_strauss_hardcore_range;
      strauss_repulsion_range = prior_strauss_repulsion_range;

    }

};





//
// PatientEstimates structure
//
//   1. PatientEstimates_Pop is the version of the 'common parms' mcmc chain
//   (patient-level estimates) that is used in the population model.
//
//   2. PatientEstimates_Single is the single-subject version, which adds two
//   parameters. These two parameters are estimated at the population level in
//   the pop model.
//

// 1. PatientEstimates_Pop
struct PatientEstimates_Pop {

  // Used in all models
  double baseline;
  double halflife;
  double errorsq;    // model error (variance)
  double mass_mean;
  double width_mean;
  int    pulse_count;

  // Population constructor
  PatientEstimates_Pop(double sv_baseline,
                       double sv_halflife,
                       double sv_errorsq,
                       double sv_mass_mean,
                       double sv_width_mean,
                       int    sv_pulse_count) {
    baseline    = sv_baseline,
    halflife    = sv_halflife;
    errorsq     = sv_errorsq;
    mass_mean   = sv_mass_mean;
    width_mean  = sv_width_mean;
    pulse_count = sv_pulse_count;

  }

  // Always use these functions to get these values.  removed them as separate
  // member variables to ensure the result is always up-to-date
  double get_decay() { return log(2) / halflife; }
  double get_logerrorsq() { return log(errorsq); }

};


// 2. PatientEstimates_Single
struct PatientEstimates_Single : PatientEstimates_Pop {
  // Not used in Population models:
  double mass_sd;
  double width_sd;

  // Single-subject constructor:
  PatientEstimates_Single(double sv_baseline,
                          double sv_halflife,
                          double sv_errorsq,
                          double sv_mass_mean,
                          double sv_width_mean,
                          int    sv_pulse_count,
                          double sv_mass_sd,
                          double sv_width_sd) :
    PatientEstimates_Pop(sv_baseline, sv_halflife, sv_errorsq, sv_mass_mean,
                         sv_width_mean, sv_pulse_count) {
      mass_sd  = sv_mass_sd;
      width_sd = sv_width_sd;
    }

};





//
// PulseEstimate structure
//   aka PulseParms -- pulse chain updated by mcmc algorithm
//   called
//   TODO: Check  on lambda, does fsh need separate definition?
//

struct PulseEstimate {

  // Variables used in all models
  double time;
  double mass;
  double width;
  double tvarscale_mass;   // variance scale for mass t-dist (eta)
  double tvarscale_width;  // variance scale for width t-dist (eta)
  arma::vec mean_contribution;

  // Used only in 2-hormone models
  //double lambda; // for fsh pulse only, denomsum - NOT SURE WHAT THIS TERM IS
                   //FOR (from Karen's code)

};



//
// PatientData structure
//
// TODO: Start here.  should it be updated to fit in with  inheritance approach?
//  (difference is probably too small to warrant it (just response_concentration
//  column))
// TODO: Then move on to patient class, then mmh, then mcmc for single subject,
// then add in gibbs.
//
struct PatientData {

  arma::vec time;
  arma::vec concentration;
  arma::vec response_concentration;
  int number_of_obs;
  double avg_period_of_obs; // in minutes
  double duration_of_obs;   // in minutes

  // Constructor for single hormone data
  PatientData(NumericVector in_time,
              NumericVector in_conc,
              int in_numobs,
              double in_period) {

    time              = as<arma::vec>(in_time);
    concentration     = as<arma::vec>(in_conc);
    number_of_obs     = time.size();
    duration_of_obs   = (time.end() - 1) - time.begin();
    avg_period_of_obs = duration_of_obs / number_of_obs;

  }

  // Constructor for two-hormone data
  //   using C++11 delegating constructors feature to avoid repeating code
  PatientData(NumericVector in_time,
              NumericVector in_conc,
              NumericVector in_responseconc,
              int in_numobs,
              double in_period) 
    : PatientData(in_time, in_conc, in_numobs, in_period) {

      response_concentration = as<arma::vec>(in_responseconc);

    }

};



#endif
