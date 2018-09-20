#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/counter.h>
#include <bp_mcmc/proposalvariance.h>
#include <bp_mcmc/mh.h>
#include <testing/catch.h>


//
// Test the Counter class
//

TEST_CASE( "Counter class", "[counter]" ) {

  Counter cnt;

  cnt.addreject();
  SECTION( "can iterate/add reject" ) {

    REQUIRE(cnt.getiter_since_reset() == 1);
    REQUIRE(cnt.getaccept() == 0);
    REQUIRE(cnt.getratio() == 0);
  }

  cnt.addaccept();
  SECTION( "can add accept" ) {

    REQUIRE( cnt.getaccept() == 1   );
    REQUIRE( cnt.getiter_since_reset() == 2 );
    REQUIRE( cnt.getratio()  == 0.5 );
  }

  cnt.resetratio();
  SECTION( "can reset" ) {

    REQUIRE( cnt.getaccept() == 0 );
    REQUIRE( cnt.getiter_since_reset() == 0 );
    REQUIRE( cnt.getratio()  == 0 ); // way to test 0
  }

  SECTION( "1-parameter version, core functions." ) {

    REQUIRE( cnt.getaccept() == 0 );
    REQUIRE( cnt.getiter_since_reset() == 0 );

    for (int i = 0; i < 100; ++i) {
      if (i % 4 == 0) cnt.addaccept();
      else cnt.addreject();
    }

    // double check math on this..
    REQUIRE( cnt.getaccept() == 25 );
    REQUIRE( cnt.getiter_since_reset() == 100 );
    REQUIRE( cnt.getratio() == 0.25 );
  }


}



//
// Test the ProposalVariance 1p class
//

TEST_CASE( "ProposalVariance class (1 param)", "[ProposalVariance]" ) {

  int iteration = 0;
  double x = 10;
  ProposalVariance pv(x, 500, 25000, 0.35);

  SECTION( "can initilize and return pv and psd values" ) {
    REQUIRE( pv.getpv() == Approx(10) );
    REQUIRE( pv.getpsd() == sqrt(10) );
  }

  pv.addreject(iteration); iteration++;
  SECTION( "can iterate/add reject" ) {

    REQUIRE( pv.getaccept() == 0 );
    REQUIRE( pv.getiter_since_reset()   == 1 );
    REQUIRE( pv.getratio()  == 0 );
  }

  pv.addaccept(iteration); iteration++;
  SECTION( "can add accept" ) {

    REQUIRE( pv.getaccept() == 1   );
    REQUIRE( pv.getiter_since_reset()   == 2 );
    REQUIRE( pv.getratio()  == 0.5 );
  }

  pv.resetratio();
  SECTION( "can reset" ) {

    REQUIRE( pv.getaccept() == 0 );
    REQUIRE( pv.getiter_since_reset()   == 0 );
    REQUIRE( pv.getratio()  == 0 );
  }

  SECTION( "1-parameter version, core functions." ) {
    ProposalVariance pv2(12, 500, 25000, 0.35);

    REQUIRE( pv2.getpv() == Approx(12) );
    REQUIRE( pv2.getpsd() == sqrt(12) );

    iteration = 0;
    while (iteration < 100) {
      if (iteration % 4 == 0) pv2.addaccept(iteration);
      else pv2.addreject(iteration);
      ++iteration;
    }

    REQUIRE( pv2.getaccept() == 25 );
    REQUIRE( pv2.getiter_since_reset()   == 100 );
    REQUIRE( pv2.getratio() == 0.25 );
  }

  SECTION( "can adjust on adjust_iter" ) {

    pv.resetratio();
    iteration = 0;
    while (iteration < 501) {
      // 500 = adjust_iter
      if ((iteration % 2) == 0) pv.addreject(iteration);
      else pv.addaccept(iteration);
      ++iteration;
    }

    REQUIRE( pv.getpv() == Approx(11) );

  }
}



//
// Test the ProposalVariance2p class
//

TEST_CASE( "ProposalVariance class (2 param)", "[ProposalVariance2p]" ) {

  int iteration = 0;
  arma::vec initial_pvs = { 0.7, 0.1 };
  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);

  SECTION( "can initialize pv and psd matrices" ) {
    double x, y, xy;
    x = 0.7; y = 0.1; xy = -0.9 * sqrt(x * y);
    arma::mat checkpv = { { x, xy }, { xy, y } };
    arma::mat checkchol = arma::chol(checkpv);

    REQUIRE( arma::approx_equal(pv.getpv(), checkpv, "absdiff", 0.0000001) );
    REQUIRE( arma::approx_equal(pv.getpsd(), checkchol, "absdiff", 0.0000001) );
  }

  SECTION( "Cholesky decomposed matrix is upper triangle form" ) {

    arma::mat cholmat = pv.getpsd();
    double upperright = cholmat(0, 1);
    double lowerleft  = cholmat(1, 0);

    REQUIRE( lowerleft == 0 );
    REQUIRE( upperright != 0 );

  }

  pv.addreject(iteration); ++iteration;
  SECTION( "can iterate/add reject" ) {

    REQUIRE( pv.getiter_since_reset()   == 1 );
    REQUIRE(pv.getaccept() == 0);
    REQUIRE(pv.getratio() == 0);

  }

  pv.addaccept(iteration); ++iteration;
  SECTION( "can add accept" ) {

    REQUIRE( pv.getaccept() == 1   );
    REQUIRE( pv.getiter_since_reset()   == 2 );
    REQUIRE( pv.getratio()  == 0.5 );
  }

  pv.resetratio();
  SECTION( "can reset" ) {

    REQUIRE( pv.getaccept() == 0 );
    REQUIRE( pv.getiter_since_reset()   == 0 );
    REQUIRE( pv.getratio()  == 0 );
  }


  SECTION( "2-parameter version, core functions." ) {

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

    REQUIRE( pv.getaccept() == 25 );
    REQUIRE( pv.getiter_since_reset()   == 100 );
    REQUIRE( pv.getratio() == 0.25 );

    // Test adjustpv();
    pv.adjustpv(-0.90);
    REQUIRE( arma::approx_equal(pv.getpv(), initial_mat, "absdiff", 0.0000001) );

    // game out case scenarios where > 0.9 and < 1.1.
    //double y = 1.0 + 1000.0 * pow(0.2 - 0.25, 3);
    //if (y < 0.9) y = 0.9;
    //if (y > 1.1) y = 1.1;
    //if (y == 1.1 | y == 0.9) {
    //  pv = pv % mydiag;
    //}

   }

}



