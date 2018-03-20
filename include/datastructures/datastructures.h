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





//
// PatientEstimates structure
//
//  TODO: abandoned inheritance approach -- clean up comments
//   1. PatientEstimates_Pop is the version of the 'common parms' mcmc chain
//   (patient-level estimates) that is used in the population model.
//
//   2. PatientEstimates_Single is the single-subject version, which adds two
//   parameters. These two parameters are estimated at the population level in
//   the pop model.
//

struct PatientEstimates {

  // Used in all models
  arma::vec baseline_halflife;
  //double baseline;
  //double halflife;
  double errorsq;    // model error (variance)
  double mass_mean;
  double width_mean;
  int    pulse_count; // function of linked list instead?;
  // Always use these functions to get these values.  removed them as separate
  // member variables to ensure the result is always up-to-date
  double get_decay() { return log(2) / baseline_halflife(1); }
  double get_logerrorsq() { return log(errorsq); }

  //
  // For population-only model (just constructor is unique)
  //
  PatientEstimates(double sv_baseline,
                   double sv_halflife,
                   double sv_errorsq,
                   double sv_mass_mean,
                   double sv_width_mean,
                   int    sv_pulse_count) {

    baseline_halflife = { sv_baseline, sv_halflife };
    errorsq     = sv_errorsq;
    mass_mean   = sv_mass_mean;
    width_mean  = sv_width_mean;
    pulse_count = 1;

  }


  //
  // For single-subject model only
  //
  double mass_sd;
  double width_sd;

  // Single-subject constructor:
  PatientEstimates(double sv_baseline,
                   double sv_halflife,
                   double sv_errorsq,
                   double sv_mass_mean,
                   double sv_width_mean,
                   int    sv_pulse_count,
                   double sv_mass_sd,
                   double sv_width_sd)
    : PatientEstimates(sv_baseline, sv_halflife, sv_errorsq, sv_mass_mean,
                       sv_width_mean, sv_pulse_count)
  {
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

class PulseEstimate {

  public:

    // Variables used in all models
    double time;
    double mass;
    double width;
    double tvarscale_mass;   // variance scale for mass t-dist (eta)
    double tvarscale_width;  // variance scale for width t-dist (eta)
    double lambda; // for fsh pulse only, denomsum - NOT SURE WHAT THIS TERM IS
                   //FOR (from Karen's code)

    // Public user-facing access function for mean_contrib. Only calculates
    // mean_contrib if input values have changed.
    // TODO: change this return type to a pointer
    arma::vec get_mean_contribution(const arma::vec & data_time, double decay_rate)
    {

      if ((prev_mass != mass) | (prev_width != width) | (prev_time != time) |
          (prev_decay_rate != decay_rate)) {
        calc_mean_contribution(data_time, decay_rate);
        prev_time = time;
        prev_mass = mass;
        prev_width = width;
        prev_decay_rate = decay_rate;
      }

      return mean_contribution;

    }

    // Constructor for new PulseEstimate objects
    PulseEstimate(double in_time,
                  double in_mass,
                  double in_width,
                  double in_tvarscale_mass,
                  double in_tvarscale_width,
                  //double fshlambda,
                  double patient_decay,
                  const arma::vec &data_time)
        : time            (in_time)
        , mass            (in_mass)
        , width           (in_width)
        , tvarscale_mass  (in_tvarscale_mass)
        , tvarscale_width (in_tvarscale_width)
        //, lambda          (fshlambda)
        , mean_contribution(data_time.n_elem)
      {
        mean_contribution.fill(0.);
        calc_mean_contribution(data_time, patient_decay);
        prev_time       = time;
        prev_mass       = mass;
        prev_width      = width;
        prev_decay_rate = patient_decay;
      }
    // Constructor for empty pulse object
    PulseEstimate()
      : time(0), mass(0), width(0), tvarscale_mass(0),
        tvarscale_width(0), mean_contribution(1) 
      {
        mean_contribution.fill(0);
      }

  private:

    arma::vec mean_contribution;
    double prev_time, prev_mass, prev_width, prev_decay_rate;

    // mean_contribution() of each pulse to the total mean_concentration
    void calc_mean_contribution(const arma::vec &data_time, double decay_rate)
    {

      double y, z, w;
      arma::vec x(data_time.n_elem);
      x.fill(0.);

      z  = width * decay_rate;
      y  = decay_rate * (0.5 * z  + time);
      z += time;
      w  = sqrt(2. * width);
      x = ((data_time - z) / w) * sqrt(2);

      double N = data_time.n_elem;
      // NOTE: potentially slow piece
      for (int i = 0; i < N; i++) {
        x(i) = Rf_pnorm5(x(i), 0.0, 1.0, 1, 0);
        //x = arma::normpdf(x); // doesn't give same results
      }

      // Finish calculating mean_contrib w/ vectorized ops
      //    mass * x is a vector, as is exp(), so use element-wise
      //    multiplication via %
      mean_contribution = (mass * x) % exp(y - data_time * decay_rate);
      // Truncate <0 = 0
      mean_contribution.for_each( [](arma::vec::elem_type& val) { val = std::max(val, 0.); } );

    }


};





//
// PatientData structure
//
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
  double fitstart;
  double fitend;

  // Constructor for single hormone data
  PatientData(NumericVector in_time,
              NumericVector in_conc) {

    time              = as<arma::vec>(in_time);
    concentration     = as<arma::vec>(in_conc);
    concentration     = log(concentration);  // store data on log scale
    number_of_obs     = time.size();
    duration_of_obs   = time(number_of_obs - 1) - time(0); // NOTE: 1430 for typical 24 hour dataset (not 1440)
    avg_period_of_obs = duration_of_obs / (number_of_obs - 1);

    fitstart = -(avg_period_of_obs * 4);
    fitend   =  time(number_of_obs - 1) + (avg_period_of_obs * 2);

  }

  // Constructor for two-hormone data
  //   using C++11 delegating constructors feature to avoid repeating code
  PatientData(NumericVector in_time,
              NumericVector in_conc,
              NumericVector in_responseconc
             ) : PatientData(in_time, in_conc) {

      response_concentration = as<arma::vec>(in_responseconc);
      response_concentration = log(response_concentration);

    }

};



#endif

