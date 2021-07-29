#include <testthat.h>
#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/counter.h>
#include <bp_mcmc/proposalvariance.h>
#include <bp_mcmc/mh.h>
// #include <testing/catch.h>

//
// Test the Counter class
//

context( "Counter" ) {
  
  Counter cnt;
  
  cnt.addreject();
  test_that( "can iterate/add reject" ) {
    
    expect_true(cnt.getiter_since_reset() == 1);
    expect_true(cnt.getaccept() == 0);
    expect_true(cnt.getratio() == 0);
  }
  
  cnt.addaccept();
  test_that( "can add accept" ) {
    
    expect_true( cnt.getaccept() == 1   );
    expect_true( cnt.getiter_since_reset() == 2 );
    expect_true( cnt.getratio()  == 0.5 );
  }
  
  cnt.resetratio();
  test_that( "can reset" ) {
    
    expect_true( cnt.getaccept() == 0 );
    expect_true( cnt.getiter_since_reset() == 0 );
    expect_true( cnt.getratio()  == 0 ); // way to test 0
  }
  
  test_that( "1-parameter version, core functions." ) {
    
    expect_true( cnt.getaccept() == 0 );
    expect_true( cnt.getiter_since_reset() == 0 );
    
    for (int i = 0; i < 100; ++i) {
      if (i % 4 == 0) cnt.addaccept();
      else cnt.addreject();
    }
    
    // double check math on this..
    expect_true( cnt.getaccept() == 25 );
    expect_true( cnt.getiter_since_reset() == 100 );
    expect_true( cnt.getratio() == 0.25 );
  }
  
  
}



//
// Test the ProposalVariance 1p class
//

context( "ProposalVariance - One Param" ) {
  
  int iteration = 0;
  double x = 10;
  ProposalVariance pv(x, 500, 25000, 0.35);
  
  test_that( "can initilize and return pv and psd values" ) {
    expect_true( pv.getpv() == Approx(10) );
    expect_true( pv.getpsd() == sqrt(10) );
  }
  
  pv.addreject(iteration); iteration++;
  test_that( "can iterate/add reject" ) {
    
    expect_true( pv.getaccept() == 0 );
    expect_true( pv.getiter_since_reset()   == 1 );
    expect_true( pv.getratio()  == 0 );
  }
  
  pv.addaccept(iteration); iteration++;
  test_that( "can add accept" ) {
    
    expect_true( pv.getaccept() == 1   );
    expect_true( pv.getiter_since_reset()   == 2 );
    expect_true( pv.getratio()  == 0.5 );
  }
  
  pv.resetratio();
  test_that( "can reset" ) {
    
    expect_true( pv.getaccept() == 0 );
    expect_true( pv.getiter_since_reset()   == 0 );
    expect_true( pv.getratio()  == 0 );
  }
  
  test_that( "1-parameter version, core functions." ) {
    ProposalVariance pv2(12, 500, 25000, 0.35);
    
    expect_true( pv2.getpv() == Approx(12) );
    expect_true( pv2.getpsd() == sqrt(12) );
    
    iteration = 0;
    while (iteration < 100) {
      if (iteration % 4 == 0) pv2.addaccept(iteration);
      else pv2.addreject(iteration);
      ++iteration;
    }
    
    expect_true( pv2.getaccept() == 25 );
    expect_true( pv2.getiter_since_reset()   == 100 );
    expect_true( pv2.getratio() == 0.25 );
  }
  
  test_that( "can adjust on adjust_iter" ) {
    
    pv.resetratio();
    iteration = 0;
    while (iteration < 501) {
      // 500 = adjust_iter
      if ((iteration % 2) == 0) pv.addreject(iteration);
      else pv.addaccept(iteration);
      ++iteration;
    }
    
    expect_true( pv.getpv() == Approx(11) );
    
  }
}



//
// Test the ProposalVariance2p class
//

context( "ProposalVariance - 2 param") {
  
  int iteration = 0;
  arma::vec initial_pvs = { 0.7, 0.1 };
  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);
  
  test_that( "can initialize pv and psd matrices" ) {
    double x, y, xy;
    x = 0.7; y = 0.1; xy = -0.9 * sqrt(x * y);
    arma::mat checkpv = { { x, xy }, { xy, y } };
    arma::mat checkchol = arma::chol(checkpv);
    
    expect_true( arma::approx_equal(pv.getpv(), checkpv, "absdiff", 0.0000001) );
    expect_true( arma::approx_equal(pv.getpsd(), checkchol, "absdiff", 0.0000001) );
  }
  
  test_that( "Cholesky decomposed matrix is upper triangle form" ) {
    
    arma::mat cholmat = pv.getpsd();
    double upperright = cholmat(0, 1);
    double lowerleft  = cholmat(1, 0);
    
    expect_true( lowerleft == 0 );
    expect_true( upperright != 0 );
    
  }
  
  pv.addreject(iteration); ++iteration;
  test_that( "can iterate/add reject" ) {
    
    expect_true( pv.getiter_since_reset()   == 1 );
    expect_true(pv.getaccept() == 0);
    expect_true(pv.getratio() == 0);
    
  }
  
  pv.addaccept(iteration); ++iteration;
  test_that( "can add accept" ) {
    
    expect_true( pv.getaccept() == 1   );
    expect_true( pv.getiter_since_reset()   == 2 );
    expect_true( pv.getratio()  == 0.5 );
  }
  
  pv.resetratio();
  test_that( "can reset" ) {
    
    expect_true( pv.getaccept() == 0 );
    expect_true( pv.getiter_since_reset()   == 0 );
    expect_true( pv.getratio()  == 0 );
  }
  
  
  test_that( "2-parameter version, core functions." ) {
    
    // Initialize object for testing
    arma::vec initial_pvs = { 0.7, 0.1 };
    ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);
    arma::mat initial_mat = pv.getpv();
    
    iteration = 0;
    while (iteration < 100) {
      if (iteration % 4 == 0) pv.addaccept(iteration);
      else (pv.addreject(iteration));
      ++iteration;
    }
    
    expect_true( pv.getaccept() == 25 );
    expect_true( pv.getiter_since_reset()   == 100 );
    expect_true( pv.getratio() == 0.25 );
    
    // Test adjustpv();
    pv.adjustpv(-0.90);
    expect_true( arma::approx_equal(pv.getpv(), initial_mat, "absdiff", 0.0000001) );
    
    // game out case scenarios where > 0.9 and < 1.1.
    //double y = 1.0 + 1000.0 * pow(0.2 - 0.25, 3);
    //if (y < 0.9) y = 0.9;
    //if (y > 1.1) y = 1.1;
    //if (y == 1.1 | y == 0.9) {
    //  pv = pv % mydiag;
    //}
    
  }
  
}

