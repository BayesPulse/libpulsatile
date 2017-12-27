
#include <RcppArmadillo.h>
#include <RInside.h>
//#include "catch.h"
#include "counter.h"
#include "proposalvariance.h"
#include "mh.h"
#include "patient.h"
#include "utils.h"

int main(int argc, char **argv) {

  std::cout << "Hello world\n";

  RInside R;

  //int i;
  //int len = 1;
  //arma::vec myvec(len);

  //for (i = 0; i < len; i++) {
  //  myvec(i) = ::Rf_rnorm(0, 1);
  //  std::cout << "Random normal variate " << i << " = " << myvec(i)  << "\n";
  //}

  pulseutils::set_seed(171227);

  arma::vec initial_means = { 2, 3 };
  arma::vec initial_pvs = { 0.7, 0.1 };
  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);

  std::cout << "rmvnorm = " << pulseutils::rmvnorm(initial_means, pv.getpsd());

  return 0;

}


