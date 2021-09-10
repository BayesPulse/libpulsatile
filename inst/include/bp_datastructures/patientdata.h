#ifndef GUARD_bp_datastructures_patientdata_h
#define GUARD_bp_datastructures_patientdata_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/utils.h>

//
// patientdata.h
//   defining the object for holding pulse level estimates
//
// Authors: Matt Mulvahill
//          Max McGrath
// Notes:
//

using namespace Rcpp;



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

    //fitstart = -(avg_period_of_obs * 4);
    fitstart = in_time(0) - (avg_period_of_obs * 5);
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

