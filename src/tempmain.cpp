
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

  // Create R instance (for RNGs) and set seed for reproducing test results
  RInside R;
  pulseutils::set_seed(171227);

  // Create objects for checking/testing RNG implementation
  arma::vec initial_means = { 2, 3 };
  arma::vec initial_pvs   = { 0.7, 0.1 };

  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);
  arma::mat cholmat = pv.getpsd();
  double cholcov    = cholmat(0, 1);

  // Test/Print RNG results
  std::cout << "PV = " << pv.getpv() << "\n";
  std::cout << "Choldecomp of PV (i.e. PSD, upper triangle form)  = " << pv.getpsd() << "\n";
  std::cout << "PSD test that actually upper triangular form = " << cholcov << "\n";
  std::cout << "rmvnorm = " << pulseutils::rmvnorm(initial_means, pv.getpsd()) << "\n";

  return 0;

}


