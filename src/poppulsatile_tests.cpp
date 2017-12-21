#include <RcppArmadillo.h>
#include <testthat.h>
#include "counter.h"
#include "proposalvariance.h"

//
// Test the ProposalVariance class
//

context( "ProposalVariance classes") {

  test_that("1-parameter version, counter implementation.") {

    //arma::dvec initial_pvs = {0.7, 0.1};
    double x = 0.3;

    ProposalVariance pv(x, 500, 25000, 0.35);
    expect_true(pv.getpv() == 0.3);

    pv.addreject();
    expect_true(pv.getiter() == 1);

    pv.addaccept();
    expect_true( pv.getaccept() == 1   );
    expect_true( pv.getiter()   == 2   );
    expect_true( pv.getratio()  == 0.5 );

    pv.resetratio();
    expect_true( pv.getaccept() == 0 );
    expect_true( pv.getiter()   == 0 );
    expect_true( pv.getratio()  != pv.getratio() ); // nan, in c++ is never equal
                                                    // to itself, so testing that
                                                    // this is nan...

  }

  // test_that("Test 2-parameter version, counter implementation.") {
  // }

  test_that("1-parameter version, core functions.") {

    ProposalVariance pv(0.7, 500, 25000, 0.35);
    expect_true(pv.getpv() == 0.7);

    int i;
    for (int i = 0; i < 100; i++) {
      if (i % 4 == 0) pv.addaccept();
      pv.addreject();
    }

    // double check math on this..
    expect_true(pv.getratio() == 0.2);

  }

}


//
// Test the Counter class
//

context( "Counter class") {

  test_that("Counter can count") {

  Counter cnt;

  cnt.addreject();
  expect_true(cnt.getiter() == 1);

  cnt.addaccept();
  expect_true( cnt.getaccept() == 1   );
  expect_true( cnt.getiter()   == 2   );
  expect_true( cnt.getratio()  == 0.5 );

  cnt.resetratio();
  expect_true( cnt.getaccept() == 0 );
  expect_true( cnt.getiter()   == 0 );
  expect_true( cnt.getratio()  != cnt.getratio() );

  }
}

