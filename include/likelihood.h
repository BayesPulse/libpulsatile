#ifndef GUARD_likelihood_h
#define GUARD_likelihood_h

//
// likelihood.h
//

#include <RcppArmadillo.h>
#include <RInside.h>
#include "patient.h"
#include "datastructures.h"

namespace pulselikelihood {

  //
  // likelihood() 
  //   computes the current likelihood using the observed log-concentrations and
  //   mean concentration
  //
  //  Args: Node_type *list; this is the current list of pulses that exist;
  //        double **ts; this is the matrix of observed data (a column of
  //                times and a column of log(concentration);
  //        Common_parms *parms; the current values of the common parameters;
  //        int N; the number of observations in **ts;
  //        Node_type *node_out; if we want, we can ignore a pulse;
  //
  double likelihood(Patient *patient,
                    Pulse *pulse_excluded,
                    bool response_hormone = FALSE)
  {

    int i;
    int N = patient.get_pulsecount();  // function not yet defined **
    double like = 0; // likelihood to be calculated
    arma::vec mean_conc;
    arma::vec *data;

    if (response_hormone) {
      data = patient->data.response_concentration;
    } else {
      data = patient->data.concentration;
    }

    // Sum across mean_contribs
    mean_conc = mean_concentration(patient, pulse_excluded);
    for (i = 0; i < N; i++) {
      like += pow(data(i) - mean_conc(i), 2);  // should be able to get rid of the loop by element diffing these vectors and squaring
    }

    like /= (-2.0 * patient->estimates.errorsq);
    like += -0.5 * N * (1.8378771 + patient->estimates.logerrorsq);

    return like;

  }


  //
  // mean_concentration()
  //   this takes each pulse's mean_contrib vector and sums across them
  //
  arma::vec mean_concentration(Patient *patient, Pulse *pulse_excluded)
  {

    int i;
    Pulse *pulse;
    int N = patient.get_pulsecount();  // function not yet defined **
    arma::vec mean_conc(N);

    // Add the contribution to the mean from each pulse
    pulse = Patient->pulses->next; // move to first pulse, not sure if works**
    while (pulse != NULL) {
      if (pulse != pulse_excluded) {
        for (i = 0; i < N; i++)
          mean_conc(i) += pulse->mean_contrib(i); // mean contrib not yet fully defined **
      }
      node = pulse->next;
    }

    // Add the baseline contribution and log
    mean_conc += baseline;
    mean_conc = log(mean_conc);

    return mean_conc;
  }


  //
  // mean_contribution()
  //   this updates a pulse's mean_contrib vector based on current values of
  //   parameters
  // This function updates a pulse's mean_contrib vector based on inputted parms
  //
  // Node_type *node; what pulse are we updating;
  // double **ts; this is the matrix of observed data (a column of times and a column of log(concentration);
  // Common_parms *parms; the current values of the common parameters;
  // int N; the number of observations in **ts;
  //
  void mean_contribution(PulseEstimate *pulse, arma::vec *concentration, PatientEstimates *patest)
  {

    int i;             // generic counter
    double x, y, z, w; // part of arithmetic used in calculating mean contrib
    double decay = patest->get_decay();

    z  = pulse->theta[1] * decay;
    y  = decay * (0.5 * z  + pulse->time);
    z += pulse->time;
    w  = sqrt(2. * pulse->theta[1]);

    for (i = 0; i < N; i++) {
      x = (ts[i][0] - z) / w;
      x = Rf_pnorm5(x * sqrt(2), 0.0, 1.0, 1, 0);

      if (x == 0) {
        pulse->mean_contrib[i] = 0; 
      } else {
        pulse->mean_contrib[i] = pulse->theta[0] * x * exp(y - ts[i][0] * decay);
      }
    }

  }

}


#endif
