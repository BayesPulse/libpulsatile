#ifndef GUARD_utils_h
#define GUARD_utils_h

#include <RcppArmadillo.h>
#include <RInside.h>

//
// utils.h
//

class PulseUtils {

  public:

    // 
    // orderstat_default()
    //   Defining in a single location the order-statistic on pulse location used
    //   in all versions of the algirithm.
    //
    int orderstat_default() { return 3; };

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

};

//----------------------------------------------------------------------------//
// END OF FILE
//----------------------------------------------------------------------------//

#endif
