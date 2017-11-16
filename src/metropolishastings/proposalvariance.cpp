#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <cmath>
#include "proposalvariance.h"

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
//


// destructor -- need to define here? - nothing so far - not yet used
//ProposalVariance::~ProposalVariance() { }


// Constructor
// TODO: - need separate constructors for Marginal and Joint
//       - Joint needs cholesky decomposition and covariance initialization for
//       proposal varcov matrix
MarginalProposalVariance::MarginalProposalVariance(double initial_pv,
                                                   int adjust_at_iter = 500,
                                                   int max_iters = 25000); 
{

  pv             = initial_pv;
  adjust_at_iter = adjust_at_iter;
  max_iters      = max_iters;
  accept_ct      = 0;
  iter_ct        = 0;
  ratio          = 0;


}



JointProposalVariance::JointProposalVariance(double initial_pv,
                                             int adjust_at_iter = 500,
                                             int max_iters = 25000); 
{

  adjust_at_iter = adjust_at_iter;
  max_iters      = max_iters;
  accept_ct      = 0;
  iter_ct        = 0;
  ratio          = 0;

  // Code for cholesky decomp and covariance initialization
  // but seems to be no cholesky decomp (pop_mcmc.c:350)
  pv(0, 0) = initial_pv(0);
  pv(1, 1) = initial_pv(1);
  pv(0, 1) = pv(1, 0) = -0.90 * sqrt(pv(0, 0) * pv(1, 1));

}




// Acceptance adjustment routine for one-parameter modified MH
// args: x   = current acceptance rate
//       *X  = current proposal variance
// returns: none -- updates internally
void MarginalProposalVariance::adjustpv(double x, double *X, double target_pv = 0.35)
{

  double y = 1.0 + 1000.0 * pow(get_ratio() - target_pv, 3);

  if (y < 0.9) {
    pv *= 0.9;
  } else if (y > 1.1) {
    pv *= 1.1;
  }

}


// Acceptance adjustment routine for two-parameter modified MH
//   used previously for baseline and halflife
// args: [desired?] corr correlation between the two variances
// returns: none -- updates internally
void JointProposalVariance::adjustpv(double corr, double target_pv = 0.25)
{

  // y - new diagonal elements of proposal variance-covariance matrix based on
  // inputs
  double y = 1.0 + 1000.0 * pow(get_ratio() - target_pv, 3);

  if (y < .90) {

    y = .90;
    pv(0, 0) *= y;
    pv(1, 1) *= y;
    pv(0, 1)  = pv(1, 0) = corr * sqrt(pv(0, 0) * pv(1, 1));

  } else if (y > 1.1) {

    y = 1.1;
    pv(0, 0) *= y;
    pv(1, 1) *= y;
    pv(0, 1)  = pv(1, 0) = corr * sqrt(pv(0, 0) * pv(1, 1));

  }

}

