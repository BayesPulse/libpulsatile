#include <RcppArmadillo.h>
#include <testthat.h>
#include "counter.h"
#include "proposalvariance.h"

//
// Test the ProposalVariance class
//

context( "ProposalVariance classes") {

  test_that("2-parameter version, counter implementation.") {

    // Initialize object for testing
    arma::vec initial_pvs = { 0.7, 0.1 };
    ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);

    // Build matrix for testing that full matrix is initialized properly from
    // vector
    double x, y, xy;
    x = 0.7; y = 0.1; xy = -0.9 * sqrt(x * y);
    arma::mat checkpv = { { x, xy }, { xy, y } };
    // Test it
    expect_true(arma::approx_equal(pv.getpv(), checkpv, "absdiff", 0.0000001));

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

  test_that("2-parameter version, core functions.") {

    // Initialize object for testing
    arma::vec initial_pvs = { 0.7, 0.1 };
    ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);
    arma::mat initial_mat = pv.getpv();

    int i;
    for (int i = 0; i < 100; i++) {
      if (i % 4 == 0) pv.addaccept();
      pv.addreject();
    }

    // double check math on this..
    expect_true(pv.getratio() == 0.2);

    // Test adjustpv();
    pv.adjustpv(-0.90);
    expect_true(arma::approx_equal(pv.getpv(), initial_mat, "absdiff", 0.0000001));

    //double y = 1.0 + 1000.0 * pow(0.2 - 0.25, 3);
    //if (y < 0.9) y = 0.9;
    //if (y > 1.1) y = 1.1;
    //if (y == 1.1 | y == 0.9) {
    //  pv = pv % mydiag;
    //}



  }


}

