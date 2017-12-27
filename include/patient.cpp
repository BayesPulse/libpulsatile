/*******************************************************************************
* PatientData functions and struct
*
* Member functions for creating PatientData structs. Not classes, so all member
* objects can be read/updated directly -- PatientData is constant, so is this
* what we want?
*******************************************************************************/
#include <RInside>
#include <RcppArmadillo>
//#include <Rcpp>
#include <testthat.h>
#include "patient.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Constructor
Patient::Patient(NumericMatrix data,
                 NumericVector priors,
                 NumericVector starting_values)
{

  read_data(data);
  PatientPriors priors();
  PatientData data();
  PatientEstimates estimates();


}


// Read in data
void Patient::read_data(NumericMatrix data)
{

  // col 0 = time; col 1 = (driver) conc.; col 2 = (response conc.)
  int ncol = data.ncol();
  if (ncol == 2) {
    read_data_in(data(_, 0), data(_, 1));
  } else if (ncol == 3) {
    read_data_in(data(_, 0),  data(_, 1),  data(_, 2));
  } else {
    // **error handling**
  }

}




// read_data() -- read data into PatientData object
void PatientData::read_data(std::vector<double> in_time,
                            std::vector<double> in_conc)
{

  time = in_time;
  conc = in_conc;

  number_of_obs     = time.size();
  duration_of_obs   = (time.end()-1) - time.begin(); // pseudo code.. use iterator?
  avg_period_of_obs = duration_of_obs / number_of_obs;

}


// read_data() -- read data into PatientData object, for driver-response data
void PatientData::read_data(std::vector<double> in_time,
                            std::vector<double> in_conc,
                            std::vector<double> in_response_conc)
{

  time = in_time;
  concentration = in_conc;
  response_concentration = in_response_conc;

  number_of_obs     = time.size();
  duration_of_obs   = (time.end()-1) - time.begin(); // pseudo code.. use iterator?
  avg_period_of_obs = duration_of_obs / number_of_obs;

}
