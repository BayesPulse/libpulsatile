#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.h"
#include <armadillo>
#include "counter.h"
#include "proposalvariance.h"


//
// Test the Counter class
//

TEST_CASE( "Counter class", "[counter]" ) {

  Counter cnt;

  cnt.addreject();
  SECTION( "can iterate/add reject" ) {

    REQUIRE(cnt.getiter() == 1);
    REQUIRE(cnt.getaccept() == 0);
    REQUIRE(cnt.getratio() == 0);
  }

  cnt.addaccept();
  SECTION( "can add accept" ) {

    REQUIRE( cnt.getaccept() == 1   );
    REQUIRE( cnt.getiter()   == 2   );
    REQUIRE( cnt.getratio()  == 0.5 );
  }

  cnt.resetratio();
  SECTION( "can reset" ) {

    REQUIRE( cnt.getaccept() == 0 );
    REQUIRE( cnt.getiter()   == 0 );
    REQUIRE( cnt.getratio()  != cnt.getratio() );
  }

  SECTION( "1-parameter version, core functions." ) {

    REQUIRE( cnt.getaccept() == 0 );
    REQUIRE( cnt.getiter() == 0 );

    for (int i = 0; i < 100; ++i) {
      if (i % 4 == 0) cnt.addaccept();
      cnt.addreject();
    }

    // double check math on this..
    REQUIRE( cnt.getaccept() == 25 );
    REQUIRE( cnt.getiter() == 100 );
    REQUIRE( cnt.getratio() == 0.25 );
  }


}




//
// Test the ProposalVariance 1p class
//

TEST_CASE( "ProposalVariance class (1 param)", "[ProposalVariance]" ) {

  double x = 0.3;
  ProposalVariance pv(x, 500, 25000, 0.35);

  SECTION( "can initilize and return pv and psd values" ) {
    REQUIRE( pv.getpv() == 0.3 );
    REQUIRE( pv.getpsd() == sqrt(0.3) );
  }

  pv.addreject();
  SECTION( "can iterate/add reject" ) {

    REQUIRE(pv.getiter() == 1);
    REQUIRE(pv.getaccept() == 0);
    REQUIRE(pv.getratio() == 0);
  }

  pv.addaccept();
  SECTION( "can add accept" ) {

    REQUIRE( pv.getaccept() == 1   );
    REQUIRE( pv.getiter()   == 2   );
    REQUIRE( pv.getratio()  == 0.5 );
  }

  pv.resetratio();
  SECTION( "can reset" ) {

    REQUIRE( pv.getaccept() == 0 );
    REQUIRE( pv.getiter()   == 0 );
    REQUIRE( pv.getratio()  != pv.getratio() );
  }

  SECTION( "1-parameter version, core functions." ) {
    ProposalVariance pv2(0.7, 500, 25000, 0.35);

    REQUIRE( pv2.getpv() == 0.7 );
    REQUIRE( pv2.getpsd() == sqrt(0.7) );

    for (int i = 0; i < 100; ++i) {
      if (i % 4 == 0) pv2.addaccept();
      pv2.addreject();
    }

    // double check math on this..
    REQUIRE( pv2.getratio() == 0.25 );
  }

}




//
// Test the ProposalVariance2p class
//

//TEST_CASE( "ProposalVariance class (2 param)", "[ProposalVariance2p]" ) {
//
//  arma::vec initial_pvs = { 0.7, 0.1 };
//  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);
//
//  SECTION( "can initialize pv and psd matrices" ) {
//    double x, y, xy;
//    x = 0.7; y = 0.1; xy = -0.9 * sqrt(x * y);
//    arma::mat checkpv = { { x, xy }, { xy, y } };
//    arma::mat checkchol = chol(checkpv);
//
//    REQUIRE( arma::approx_equal(pv.getpv(), checkpv, "absdiff", 0.0000001) );
//    REQUIRE( arma::approx_equal(pv.getpsd(), checkchol, "absdiff", 0.0000001) );
//  }
//
//  SECTION( "can iterate/add reject" ) {
//    pv.addreject();
//
//    REQUIRE(pv.getiter() == 1);
//    REQUIRE(pv.getaccept() == 0);
//    REQUIRE(pv.getratio() == 0);
//  }
//
//  SECTION( "can add accept" ) {
//    pv.addaccept();
//
//    REQUIRE( pv.getaccept() == 1   );
//    REQUIRE( pv.getiter()   == 2   );
//    REQUIRE( pv.getratio()  == 0.5 );
//  }
//
//  SECTION( "can reset" ) {
//    pv.resetratio();
//
//    REQUIRE( pv.getaccept() == 0 );
//    REQUIRE( pv.getiter()   == 0 );
//    REQUIRE( pv.getratio()  != pv.getratio() );
//  }
//
//
//  SECTION( "2-parameter version, core functions." ) {
//
//    // Initialize object for testing
//    arma::vec initial_pvs = { 0.7, 0.1 };
//    ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);
//    arma::mat initial_mat = pv.getpv();
//
//    int i;
//    for (int i = 0; i < 100; i++) {
//      pv.addreject();
//      if (i % 4 == 0) pv.addaccept();
//    }
//
//    // // double check math on this..
//    // REQUIRE( pv.getaccept() == 20 );
//    // REQUIRE( pv.getiter() == 100 );
//    // REQUIRE( pv.getratio() == 0.2 );
//
//    // // Test adjustpv();
//    // pv.adjustpv(-0.90);
//    // REQUIRE( arma::approx_equal(pv.getpv(), initial_mat, "absdiff", 0.0000001) );
//
//    // //double y = 1.0 + 1000.0 * pow(0.2 - 0.25, 3);
//    // //if (y < 0.9) y = 0.9;
//    // //if (y > 1.1) y = 1.1;
//    // //if (y == 1.1 | y == 0.9) {
//    // //  pv = pv % mydiag;
//    // //}
//
//  }
//
//}
