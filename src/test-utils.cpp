#include <testthat.h>
#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/utils.h>
#include <bp_mcmc/proposalvariance.h>
// #include <testing/catch.h>

//
// utils_tests.cpp
//     Test utility functions
//

PulseUtils pu;


context( "orderstat_default" ) {
  
  test_that( "equal to 3" ) {
    
    expect_true( pu.orderstat_default() == 3 );
    
  }
  
}

context( "rmvnorm Function" ) {
  
  arma::vec initial_means = { 2, 3 };
  arma::vec initial_pvs = { 0.7, 0.1 };
  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);
  
  test_that( "can generate random value" ) {
    
    arma::vec answer = { 1.3264, 3.1485 };
    
    pu.set_seed(171227);
    expect_true(
      arma::approx_equal(pu.rmvnorm(initial_means, pv.getpsd()),
                         answer, "absdiff", 0.0001) 
    );
  }
  
}

// context( "one_rmultinom Function" ) {
//   
//   arma::vec not_cumprobs = { 0.1, 0.05, 0.02, 0.03, 0.8 };
//   arma::vec cumprobs = { 0.1, 0.15, 0.17, 0.2, 1.0 };
//   
//   // currently rmultinom accepts non-probability vectors -- don't sum to one --
//   // look into this.
//   //test_that( "fails on non-cumulative probability vector" ) {
//   //  //pu.set_seed(171227);
//   //  expect_true_THROWS( pu.one_rmultinom(not_cumprobs) );
//   //}
//   
//   std::cout << "multinom cum prob" << cumprobs << std::endl;
//   test_that( "Succeeds on cumulative probability vector" ) {
//     //pu.set_seed(171227);
//     REQUIRE_NOTHROW( pu.one_rmultinom(not_cumprobs) );
//   }
//   std::cout << "passed multinom test" << std::endl;
//   
// }
