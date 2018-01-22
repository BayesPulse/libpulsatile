#include <RcppArmadillo.h>
#include <RInside.h>
#include "utils.h"
#include "proposalvariance.h"
#include "catch.h"


//
// utils_tests.cpp
//     Test utility functions
//

PulseUtils pu;


TEST_CASE( "orderstat_default", "[utils]" ) {

  SECTION( "equal to 3" ) {

    REQUIRE(pu.orderstat_default() == 3);

  }

}

TEST_CASE( "rmvnorm Function", "[utils]" ) {

  RInside R;

  arma::vec initial_means = { 2, 3 };
  arma::vec initial_pvs = { 0.7, 0.1 };
  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);

  SECTION( "can generate random value" ) {

    arma::vec answer = { 1.3264, 3.1485 };

    pu.set_seed(171227);
    REQUIRE(
            arma::approx_equal(pu.rmvnorm(initial_means, pv.getpsd()),
                               answer, "absdiff", 0.0001) 
            );
  }

}


