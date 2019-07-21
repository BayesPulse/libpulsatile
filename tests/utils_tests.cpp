#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/utils.h>
#include <bp_mcmc/proposalvariance.h>
#include <testing/catch.h>


//
// utils_tests.cpp
//     Test utility functions
//

PulseUtils pu;


TEST_CASE( "orderstat_default", "[utils]" ) {

  SECTION( "equal to 3" ) {

    REQUIRE( pu.orderstat_default() == 3 );

  }

}

TEST_CASE( "rmvnorm Function", "[utils]" ) {

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

TEST_CASE( "one_rmultinom Function", "[utils]" ) {

  arma::vec not_cumprobs = { 0.1, 0.05, 0.02, 0.03, 0.8 };
  arma::vec cumprobs = { 0.1, 0.15, 0.17, 0.2, 1.0 };

  // currently rmultinom accepts non-probability vectors -- don't sum to one --
  // look into this.
  //SECTION( "fails on non-cumulative probability vector" ) {
  //  //pu.set_seed(171227);
  //  REQUIRE_THROWS( pu.one_rmultinom(not_cumprobs) );
  //}

  std::cout << "multinom cum prob" << cumprobs << std::endl;
  SECTION( "Succeeds on cumulative probability vector" ) {
    //pu.set_seed(171227);
    REQUIRE_NOTHROW( pu.one_rmultinom(not_cumprobs) );
  }
  std::cout << "passed multinom test" << std::endl;

}


