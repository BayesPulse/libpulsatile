#ifndef GUARD_patient_h

#include <RcppArmadillo.h>
#include <Rcpp.h>

// patient.h
struct PatientData {

  arma::vec concentration;
  arma::vec time;
  int number_of_obs;
  double period_of_obs;
  double duration_of_obs = number_of_obs * period_of_obs;

};

// always const when used
struct PatientPriors {

    double npulses;   // prior number of pulses
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

// aka PatientParms
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

};

// linked list -- 
struct PulseEstimate {

  double time;
  double mass;
  double width;
  vector<double> mean_contribution(PatientData.number_of_obs);
  // the length needs to be moved -- likely need a 'Patient' class to do some of
  // the calculations required to enforce some of these rules.

};


#endif
