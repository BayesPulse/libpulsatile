#ifndef GUARD_patient_h
#define GUARD_patient_h

#include <armadillo>
#include <RInside.h>                    // for the embedded R via RInside
//#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <list>
//using namespace Rcpp;

//
// patient.h
//   defining the patient class and subclasses
//
// Author: Matt Mulvahill
// Created: 11/14/17
// Notes:
//

class Patient {

  public:

    Patient(NumericMatrix data,
            NumericVector priors,
            NumericVector starting_values);
    ~Patient();
    void read_data(NumericMatrix data);

  private:

    double likelihood;
    PatientPriors priors;
    PatientData data;
    PatientEstimates estimates;
    AssocEstimates association;
    int pulse_count;
    list<PulseEstimate> pulses;
    void read_data_in(arma::vec in_time,
                      arma::vec in_conc);
    void read_data_in(arma::vec in_time,
                      arma::vec in_conc,
                      arma::vec in_response_conc);


};

struct PatientData {

  arma::vec time;
  arma::vec concentration;
  arma::vec response_concentration;
  int number_of_obs;
  double avg_period_of_obs;
  double duration_of_obs = number_of_obs * period_of_obs;

};

// Patient priors and const when used in single subject model
// Estimated population chain 
struct PopulationEstimates {

  double num_pulses;   // prior number of pulses
  double baseline_mean;
  double baseline_variance;
  double halflife_mean;
  double halflife_variance;
  double mass_mean;
  double mass_variance;
  double width_mean;
  double width_variance;
  double mass_sdmax;
  double width_sdmax;
  double error_alpha;
  double error_beta;
  double num_orderstat;
  double strauss_gamma;
  double strauss_range;

};

// aka PatientParms -- values updated by mcmc algorithm
struct PatientEstimates {

  double baseline;
  double halflife;
  double decay;     // decay rate converted from above half-life
  double error2;    // model error (variance)
  double logerror2; // log of model error (may not be used)
  double mass;
  double width;
  double mass_sd;
  double width_sd;
  int pulse_count;

};

// linked list -- values updated by algorithm and objects created/destroyed by
//                birth-death
struct PulseEstimate {

  double time;
  double mass;
  double width;
  arma::vec mean_contribution(PatientData.number_of_obs);
  // the length needs to be moved -- likely need a 'Patient' class to do some of
  // the calculations required to enforce some of these rules.

};


struct AssocEstimates {

  // Parameters in the cox process for the response hormone intensity
  double cluster_size;      // cluster size (rho)
  double log_cluster_size;  // log scale cluster size
  double cluster_width;     // cluster width (nu)
  double log_cluster_width;

};


// Chains -- use R data types?
struct PulseChain {

};

struct PatientChain {

};

struct PopulationChain {

};

#endif
