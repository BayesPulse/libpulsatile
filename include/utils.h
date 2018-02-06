#ifndef GUARD_utils_h
#define GUARD_utils_h

#include <RcppArmadillo.h>
#include <RInside.h>


//
// utils.h
//    Defines the PulseUtils class, containing miscellaneous functions and
//    calculations needed for estimation.
//

class PulseUtils {

  public:

    // 
    // orderstat_default()
    //   Defining in a single location the order-statistic on pulse location used
    //   in all versions of the algirithm. (always 3, but if need to be changed,
    //   here is where it should be done.)
    //
    int orderstat_default() { return 3; };


    //
    // rmvnorm()
    //   Random multivariate normal.  Accepts a mean vector and a
    //   cholesky-decomposed variance-covariance matrix (psd from PV class). And
    //   returns a vector of results. Currently only works with bivariate
    //   normal, but could be generatlized if needed.
    //
    arma::vec rmvnorm(arma::vec mean,    // mean of distr to sample from (curr values)
                      arma::mat cholvar) // cholesky decomposed varcovar matrix (psd)
    {

      unsigned int i, j;
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

    //
    // calc_sr_strauss()
    //   Calculates sum(S(R)), the exponent on the gamma parameter in the Strauss
    //   process/prior for pulse location. Used for Strauss prior in birth_death
    //   and mh_time_strauss.
    //
    int calc_sr_strauss(double location,     // location to test pulses against
                        Patient *patient,
                        PulseEstimate *pulse_excluded) {

      int s_r = 0;       // Sum of indicators where diff < 20
      double difference; // Time difference
      std::list<PulseEstimate>::const_iterator this_piter = patient.pulses.begin();

      while (this_piter != NULL) {
        if (this_piter != excl_pulse) {
          // skip if node is same that location is from;
          difference = fabs(location - node->time);
          // increment by 1 if diff<R
          s_r = (difference < priors->range) ? s_r + 1 : s_r; 
        }
        this_piter++;
      }

      // sum(S(R)) - scalar value for sum of # pulses too close to each other
      return(s_r); 

    }


    //
    // set_seed() 
    //   set R seed from C++ for unit testing functions using R's RNGs
    //
    void set_seed(unsigned int seed) {

      Rcpp::Environment base_env("package:base");
      Rcpp::Function set_seed_r = base_env["set.seed"];
      set_seed_r(seed);

    }

};

//----------------------------------------------------------------------------//
// END OF FILE
//----------------------------------------------------------------------------//

#endif
