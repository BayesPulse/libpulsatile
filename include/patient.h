#ifndef GUARD_patient_h
#define GUARD_patient_h

#include <RcppArmadillo.h>
#include <RInside.h>                    // for the embedded R via RInside
//#include <list>

//
// patient.h
//   defining the patient class and subclasses
//
// Author: Matt Mulvahill
// Created: 11/14/17
// Notes:
//

using namespace Rcpp;
using namespace RcppArmadillo;


//------------------------------------------------------------------------------
// Data structures
//------------------------------------------------------------------------------


//
// Patient priors and const when used in single subject model
// PopulationEstimates when in the pop model (via typedef)
//
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

  // Only used in Population models:
  // Haven't quite wrapped my head around these (is it sd or variance), beleive
  // they correspond to mass_sd and width_sd in patient estimates for
  // single-subj model.
  double mass_mean_sd;
  double mass_mean_variance;
  double width_mean_sd;
  double width_mean_variance;

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


  // TODO: Will need constructor for single-subject, population, and possibly
  // separate ones for trigger and response hormones (so up to 4)

  // Base constructor
  PatientPriors(double prior_prior_baseline_mean,
                double prior_baseline_variance,
                double prior_halflife_mean,
                double prior_halflife_variance,
                double prior_mass_mean,
                double prior_mass_variance,
                double prior_width_mean,
                double prior_width_variance) {

    prior_baseline_mean = prior_prior_baseline_mean
    baseline_variance   = prior_baseline_variance
    halflife_mean       = prior_halflife_mean
    halflife_variance   = prior_halflife_variance
    mass_mean           = prior_mass_mean
    mass_variance       = prior_mass_variance
    width_mean          = prior_width_mean
    width_variance      = prior_width_variance

  }

  // Population constructor
  PatientPriors(double prior_prior_baseline_mean,
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
    PatientPriors(prior_prior_baseline_mean,
                  prior_baseline_variance,
                  prior_halflife_mean,
                  prior_halflife_variance,
                  prior_mass_mean,
                  prior_mass_variance,
                  prior_width_mean,
                  prior_width_variance) {

      mass_mean_sd  = prior_mass_mean_sd;
      width_mean_sd = prior_width_mean_sd;

    }

  // Single-subject constructor
  PatientPriors(double prior_prior_baseline_mean,
                double prior_baseline_variance,
                double prior_halflife_mean,
                double prior_halflife_variance,
                double prior_mass_mean,
                double prior_mass_variance,
                double prior_width_mean,
                double prior_width_variance,
                double prior_width_sd_max,
                double prior_error_alpha,
                double prior_error_beta,
                double prior_num_orderstat,
                int    prior_pulse_count,
                double prior_strauss_repulsion,
                double prior_strauss_hardcore_range,
                double prior_strauss_repulsion_range
               ) :
    PatientPriors(prior_prior_baseline_mean, prior_baseline_variance,
                  prior_halflife_mean, prior_halflife_variance, prior_mass_mean,
                  prior_mass_variance, prior_width_mean, prior_width_variance) {

      mass_sd_max  = prior_mass_sd_max;
      width_sd_max = prior_width_sd_max;
      error_alpha  = prior_error_alpha;
      error_beta   = prior_error_beta;

      num_orderstat = pulseutils::orderstat_default();

      pulse_count             = prior_pulse_count;
      strauss_repulsion       = prior_strauss_repulsion;
      strauss_hardcore_range  = prior_strauss_hardcore_range;
      strauss_repulsion_range = prior_strauss_repulsion_range;


    }

};

// aka PatientParms -- common chain updated by mcmc algorithm
// Called Common_parms in my single-subject
struct PatientEstimates {

  // Used in all models
  double baseline;
  double halflife;
  double decay = log(2) / halflife; // decay rate converted from above half-life
  double errorsq;                   // model error (variance)
  double logerrorsq;                // log of model error (may not be used)
  double mass_mean;
  double width_mean;
  int    pulse_count;

  // Not used in Population models:
  double mass_sd;
  double width_sd;

  // Population constructor
  PatientEstimates(double sv_baseline,
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

    decay      = log(2) / halflife;
    logerrorsq = log(errorsq);
  }

  // Single-subject constructor:
  PatientEstimates(double sv_baseline,
                   double sv_halflife,
                   double sv_errorsq,
                   double sv_mass_mean,
                   double sv_width_mean,
                   int    sv_pulse_count,
                   double mass_sd,
                   double width_sd) :
    PatientEstimates(sv_baseline, sv_halflife, sv_errorsq, sv_mass_mean,
                     sv_width_mean, sv_pulse_count) {
      mass_sd  = sv_mass_sd;
      width_sd = sv_width_sd;
    }

};


// aka PulseParms -- pulse chain updated by mcmc algorithm
// called
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



//--------------------------------------
// PatientData structure
//--------------------------------------
struct PatientData {

  vec time;
  vec concentration;
  vec response_concentration;
  int number_of_obs;
  double avg_period_of_obs; // in minutes
  double duration_of_obs;   // in minutes

  // Constructor for single hormone data
  PatientData(NumericVector in_time,
              NumericVector in_conc,
              int in_numobs,
              double in_period) {

    time              = as<vec>(in_time);
    concentration     = as<vec>(in_conc);
    number_of_obs     = as<int>(in_numobs);
    avg_period_of_obs = as<double>(in_period);
    duration_of_obs   = number_of_obs * avg_period_of_obs;

  }

  // Constructor for two-hormone data
  //   using C++11 delegating constructors feature to avoid repeating code
  PatientData(NumericVector in_time,
              NumericVector in_conc,
              NumericVector in_responseconc,
              int in_numobs,
              double in_period) 
    : PatientData(in_time, in_conc, in_numobs, in_period) {

      response_concentration = as<vec>(in_responseconc);

    }

};



//--------------------------------------
// Patient class
//--------------------------------------

class Patient {

  public:

    Patient(NumericMatrix data,
            NumericVector priors,
            NumericVector starting_values);
    ~Patient();
    void read_data(NumericMatrix data);
    int get_pulsecount() { return estimates.pulse_count };

  private:

    // need access to pulses, pulse_count, priors, data.. basically all of this
    //double likelihood; // calculate as-needed?
    PatientPriors priors;
    PatientData data;
    PatientEstimates estimates;
    AssocEstimates association; // Note: should move to Population class, check how its estimated.
    std::list<PulseEstimate> pulses;
    void read_data_in(arma::vec in_time,
                      arma::vec in_conc);
    void read_data_in(arma::vec in_time,
                      arma::vec in_conc,
                      arma::vec in_response_conc);

};




////////////////////////////////////////////////////////////
// Chains -- use R data types? YES
////////////////////////////////////////////////////////////
//struct PulseChain { NumericMatrix pulse };
//struct PatientChain { NumericMatrix patient };
//struct PopulationChain { NumericMatrix population };
//struct AssociationChain { NumericMatrix association };

#endif
