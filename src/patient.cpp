/*******************************************************************************
* PatientData functions and struct
*
* Member functions for creating PatientData structs. Not classes, so all member
* objects can be read/updated directly -- PatientData is constant, so is this
* what we want?
*******************************************************************************/
#include <RcppArmadillo>
#include <Rcpp>
#include "patient.h"

// [[Rcpp::depends(RcppArmadillo)]]

void PatientData::read_data(std::vector<double> in_time,
                            std::vector<double> in_conc)
{
  time = in_time;
  conc = in_conc;

  number_of_obs     = time.size();
  duration_of_obs   = (time.end()-1) - time.begin(); // pseudo code.. use iterator?
  avg_period_of_obs = duration_of_obs / number_of_obs;
}


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
