//// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>
//#include <cmath>
//#include "proposalvariance.h"
//#include <cstdio>

//
// proposalvariance.cpp
//   Methods for the ProposalVariance class
//
// Author: Matt Mulvahill
// Created: 10/13/17
// Notes:
//   Outstanding questions:
//    - Where does this class get implemented? -- implementation handled here,
//    instantiation handled in either mcmc function or MMH classes?
//    - draw_proposal() takes SD's - straighten this out.
//    - are defined destructors needed?




