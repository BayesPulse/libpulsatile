#ifndef GUARD_utils_h
#define GUARD_utils_h

//
// calculations.c
//

#include <RcppArmadillo.h>
#include <RInside.h>

namespace pulseutils {


  // 
  // orderstat_default()
  //   Defining in a single location the order-statistic on pulse location used
  //   in all versions of the algirithm.
  //
  int orderstat_default() { return 3 };

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
                    //Node_type *list,
                    //double **ts,
                    //Common_parms *parms,
                    //~~int N,
                    //Node_type *node_out,
                    //double baseline)
  {

    int i;
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

    like /= (-2.0 * patient->estimates.error2);
    like += -0.5 * N * (1.8378771 + patient->estimates.logerror2);

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
  // rmvnorm
  //
  arma::vec rmvnorm(arma::vec mean,    // mean of distr to sample from (curr values)
                    arma::mat cholvar) // cholesky decomposed varcovar matrix (psd)
  {

    int i, j;
    arma::vec runiv = { ::Rf_rnorm(0, 1), ::Rf_rnorm(0, 1) };
    arma::vec result = { 0, 0 };

    for (i = 0; i < mean.n_elem; i++) {
      for (j = 0; j <= i; j++)  result[i] += cholvar(i, j) * runiv(j);
      result(i) += mean(i);
    }

    return result;

  }


  //
  // Single random multinomial
  //  TODO: UPDATE for C++
  //
  //int one_rmultinom(double *cumprobs, int n_probs) {

  //  int i;
  //  int rtn = 0;
  //  int *ans;
  //  ans = (int *)calloc(n_probs, sizeof(int));
  //  double *probs;
  //  probs = (double *)calloc(n_probs, sizeof(double));

  //  for (i = 0; i < n_probs; i++) {
  //    if (i == 0) probs[i] = cumprobs[i];
  //    else probs[i] = cumprobs[i] - cumprobs[i-1];
  //    ans[i] = 0;
  //  }

  //  Rf_rmultinom(1, probs, n_probs, ans);

  //  for (i = 0; i < n_probs; i++) {
  //    if (ans[i] == 1) rtn = i;
  //  }

  //  return(rtn);
  //}


  //
  // set_seed -- set R seed from C++ for unit testing functions using R's RNGs
  //

  void set_seed(unsigned int seed) {

    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);

  }

}

//----------------------------------------------------------------------------//
// END OF FILE
//----------------------------------------------------------------------------//

#endif

