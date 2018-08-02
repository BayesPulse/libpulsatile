#ifndef GUARD_bp_mcmc_utils_h
#define GUARD_bp_mcmc_utils_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif


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
    //
    int one_rmultinom(arma::vec cumprobs) {

      int i;
      int rtn = 0;
      int n_probs = cumprobs.size();
      //std::cout << "size of cumprobs = " << n_probs << std::endl;
      arma::Col<int> ans(n_probs);
      arma::Col<double> probs(n_probs);

      for (i = 0; i < n_probs; i++) {

        if (i == 0) probs(i) = cumprobs(i);
        else probs(i) = cumprobs(i) - cumprobs(i-1);

        //std::cout << "i = " << i << " and cumprob = " << cumprobs(i) << "and prob = " << probs(i) << std::endl;
      }

      ::Rf_rmultinom(1, probs.begin(), n_probs, ans.begin());
      //std::cout << "ans = " << ans << std::endl;

      for (int j = 0; j < n_probs; j++) {
        //std::cout << "j = " << j << " and ans = " << ans(j) << std::endl;
        if (ans(j) == 1) rtn = j;
      }

      //std::cout << "return = " << rtn << std::endl;
      return rtn;
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
