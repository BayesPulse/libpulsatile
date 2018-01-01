#include <RcppArmadillo.h>
#include <RInside.h>
#include "datastructures.h"
#include "utils.h"
#include "proposalvariance.h"
#include "catch.h"

//
// datastructures_tests.cpp
//  For now, also contains tests of utils.h
//


//
// Test utility functions
//

TEST_CASE( "orderstat_default", "[utils]" ) {

  SECTION( "equal to 3" ) {

    REQUIRE(pulseutils::orderstat_default() == 3);

  }

}

TEST_CASE( "rmvnorm Function", "[utils]" ) {

  RInside R;

  arma::vec initial_means = { 2, 3 };
  arma::vec initial_pvs = { 0.7, 0.1 };
  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);

  SECTION( "can generate random value" ) {

    arma::vec answer = { 1.3264, 3.1485 };

    pulseutils::set_seed(171227);
    REQUIRE(
            arma::approx_equal(pulseutils::rmvnorm(initial_means, pv.getpsd()),
                               answer, "absdiff", 0.0001) 
            );
  }

}


//
// Test datastructures structs
//

TEST_CASE( "PatientPriors base constructor works (never otherwise used directly)", 
           "[datastructures]" ) {

  PatientPriors ppriors(1.5, 100, 45, 100, 3.5, 100, 30, 100);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    // Test use pointers just for fun.
    //PatientPriors * pppointer = &ppriors;

    REQUIRE(ppriors.baseline_mean     == 1.5);
    REQUIRE(ppriors.baseline_variance == 100.0);
    REQUIRE(ppriors.halflife_mean     == 45.0);
    REQUIRE(ppriors.halflife_variance == 100.0);
    REQUIRE(ppriors.mass_mean         == 3.5);
    REQUIRE(ppriors.mass_variance     == 100.0);
    REQUIRE(ppriors.width_mean        == 30.0);
    REQUIRE(ppriors.width_variance    == 100.0);

  }

  SECTION( "Other variables are initialized to 0 by default." ) {

    REQUIRE(ppriors.mass_mean_sd      == 0.0);
    REQUIRE(ppriors.width_mean_sd     == 0.0);
    REQUIRE(ppriors.mass_sd_max       == 0.0);
    REQUIRE(ppriors.width_sd_max      == 0.0);
    REQUIRE(ppriors.error_alpha       == 0.0);
    REQUIRE(ppriors.error_beta        == 0.0);
    REQUIRE(ppriors.pulse_count       == 0);
    REQUIRE(ppriors.strauss_repulsion == 0.0);
    REQUIRE(ppriors.strauss_repulsion_range == 0.0);

  }

}


//TEST_CASE( "PatientPriors Population constructor works", 
//           "[datastructures]" ) {
//
//  PatientPriors pppop(1.5, 100, 45, 100, 3.5, 100, 30, 100, 5, 10);
//
//  SECTION( "Variables included in constructor are initialized as expected." ) {
//
//    REQUIRE(pppop.baseline_mean     == 1.5);
//    REQUIRE(pppop.baseline_variance == 100.0);
//    REQUIRE(pppop.halflife_mean     == 45.0);
//    REQUIRE(pppop.halflife_variance == 100.0);
//    REQUIRE(pppop.mass_mean         == 3.5);
//    REQUIRE(pppop.mass_variance     == 100.0);
//    REQUIRE(pppop.width_mean        == 30.0);
//    REQUIRE(pppop.width_variance    == 100.0);
//    REQUIRE(pppop.mass_mean_sd      == 5.0);
//    REQUIRE(pppop.width_mean_sd     == 10.0);
//
//  }
//
//  SECTION( "Other variables are initialized to 0 by default." ) {
//
//    REQUIRE(pppop.mass_sd_max       == 0);
//    REQUIRE(pppop.width_sd_max      == 0);
//    REQUIRE(pppop.error_alpha       == 0);
//    REQUIRE(pppop.error_beta        == 0);
//    REQUIRE(pppop.pulse_count       == 0);
//    REQUIRE(pppop.strauss_repulsion == 0);
//    REQUIRE(pppop.strauss_repulsion_range == 0);
//
//  }
//
//}
//
//
//TEST_CASE( "PatientPriors single-subject constructor works",
//           "[datastructures]" ) {
//
//  PatientPriors ppsingle(1.5, 100, 45, 100, 3.5, 100, 30, 100,
//                         10, 100, 1000, 1000, 12, 0, 40);
//
//  SECTION( "Variables included in constructor are initialized as expected." ) {
//
//    REQUIRE(ppsingle.baseline_mean     == 1.5);
//    REQUIRE(ppsingle.baseline_variance == 100.0);
//    REQUIRE(ppsingle.halflife_mean     == 45.0);
//    REQUIRE(ppsingle.halflife_variance == 100.0);
//    REQUIRE(ppsingle.mass_mean         == 3.5);
//    REQUIRE(ppsingle.mass_variance     == 100.0);
//    REQUIRE(ppsingle.width_mean        == 30.0);
//    REQUIRE(ppsingle.width_variance    == 100.0);
//    REQUIRE(ppsingle.mass_sd_max       == 10);
//    REQUIRE(ppsingle.width_sd_max      == 100);
//    REQUIRE(ppsingle.error_alpha       == 1000);
//    REQUIRE(ppsingle.error_beta        == 1000);
//    REQUIRE(ppsingle.pulse_count       == 12);
//    REQUIRE(ppsingle.strauss_repulsion == 0);
//    REQUIRE(ppsingle.strauss_repulsion_range == 40);
//
//  }
//
//  SECTION( "Other variables are initialized to 0 by default." ) {
//
//    REQUIRE(ppsingle.mass_mean_sd      == 0);
//    REQUIRE(ppsingle.width_mean_sd     == 0);
//
//  }
//
//}


