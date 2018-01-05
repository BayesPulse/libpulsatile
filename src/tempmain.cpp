
#include <RcppArmadillo.h>
#include <RInside.h>
//#include "catch.h"
#include "counter.h"
#include "proposalvariance.h"
#include "mh.h"
#include "patient.h"
#include "population.h"
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

  //double a = 11;
  //arma::vec testvec = initial_means;
  //testvec += a;
  //std::cout << "vector += scalar: " << initial_means << " += " << a <<  " = " << testvec << "\n";
  //// TRUE

  // Testing using R objects in c++
  NumericVector time(100); //, conc(100);
  NumericVector conc     = rnorm(100, 3, 0.1);
  NumericVector response = conc + rnorm(100, 1, 0.001);
  for (int i = 0; i < time.size(); i++) time(i) = i + 1;


  std::cout << "Time meta: size=" << time.size() << ", begin=" << time(0) << ", end=" << time(time.size() - 1) << "\n";

  //std::cout << "NumericVector w/ values from 1 to 100 and 100 random normals and 100 more\n";
  //for (int i = 0; i < time.size(); i++) {
  //  std::cout << time(i) << ", " << conc(i) << ", " << response(i) << "\n";
  //}


  return 0;

}


