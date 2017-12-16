#include <testthat.h>
#include <RcppArmadillo.h>
#include "proposalvariance.h"

//
// Test the ProposalVariance class
//

context( "ProposalVariance class") {

  test_that("Test 2-parameter version") {

  //arma::dvec initial_pvs = {0.7, 0.1};
  double x = 0.3;

  ProposalVariance <double, double> pv(x, 200, 15000, 0.25);
  expect_true(pv.getpv() == 0.3);

  //pv.counter.addreject();
  //expect_true(pv.counter.getiter() == 1);

  //cnt.addaccept();
  //expect_true( cnt.getaccept() == 1   );
  //expect_true( cnt.getiter()   == 2   );
  //expect_true( cnt.getratio()  == 0.5 );

  //cnt.resetratio();
  //expect_true( cnt.getaccept() == 0 );
  //expect_true( cnt.getiter()   == 0 );
  //expect_true( cnt.getratio()  != cnt.getratio() );

  }
}

