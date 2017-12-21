// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "proposalvariance.h"


//--------------------------------------------------------//
// ProposalVariance (1-parm)
//--------------------------------------------------------//

// -------------------------------------
// initialize_proposals()
// -------------------------------------
// One for one-variable and one for two-variable MMH/pv.
// create/set pv and psd objects (proposal var-covar and correlation matrices)
void ProposalVariance::initialize_proposals(double initial_pv)
{

  set_proposals(initial_pv);

}


// -------------------------------------
// set_proposals()
// -------------------------------------
void ProposalVariance::set_proposals(double this_pv)
{

  pv = this_pv;
  psd = sqrt(pv);

}


// -------------------------------------
// adjustpv()
// -------------------------------------
// Acceptance adjustment routine for one-parameter modified MH
void ProposalVariance::adjustpv()
{

  double y = 1.0 + 1000.0 * pow(getratio() - target_ratio, 3);
  if (y < 0.9)      set_proposals(pv * 0.9);
  else if (y > 1.1) set_proposals(pv * 1.1);

}




//--------------------------------------------------------//
// ProposalVariance2p (2-parm)
//--------------------------------------------------------//

//--------------------------------------
// initialize_proposals() for each type
//--------------------------------------

ProposalVariance2p::ProposalVariance2p()
  : count(), adjust_iter(500), max_iter(25000), target_ratio(0.35)
{

  arma::vec in_pv(2);
  in_pv.fill(0);
  initialize_proposals(in_pv);

}


ProposalVariance2p::ProposalVariance2p(arma::vec in_pv,
                                       int in_adjust_iter, // adjust pv on multiples of adjust_iter
                                       int in_max_iter,    // maximum iteration to adjust pv
                                       double in_target_ratio) 
{

  adjust_iter  = in_adjust_iter;
  max_iter     = in_max_iter;
  target_ratio = in_target_ratio;

  initialize_proposals(in_pv);

}


//--------------------------------------
// initialize_proposals() for each type
//--------------------------------------

void ProposalVariance2p::initialize_proposals(arma::vec initial_pv) {

  arma::mat this_pv = arma::diagmat(initial_pv);
  set_proposals(this_pv, -0.90);

}


//--------------------------------------
// set_proposals()
//   Consider renaming this. It always takes in a diagonal matrix and fills in
//   the off-diagonal. (And calculates/creates the chol decomp matrix.)
//--------------------------------------
void ProposalVariance2p::set_proposals(arma::mat this_pv, double corr) {

  pv = this_pv;
  pv(0, 1) = pv(1, 0) = calc_covariance(pv, corr);
  psd = chol(pv);

}

double ProposalVariance2p::calc_covariance(arma::mat pv, double corr) {

  return corr * sqrt(pv(0, 0) * pv(1, 1));

}


//--------------------------------------
// adjustpv()
//--------------------------------------

void ProposalVariance2p::adjustpv(double corr = -0.90)
{

  // identity matrix
  arma::mat mydiag(2, 2, arma::fill::eye);

  // y - new diagonal elements of proposal variance-covariance matrix based
  // on inputs
  double y = 1.0 + 1000.0 * pow(getratio() - target_ratio, 3);

  if (y < .90) {

    y = .90;
    pv  = pv % mydiag; // set off-diag to 0 with Schur/Hadamard multiplication
    set_proposals(pv * y, corr);

  } else if (y > 1.1) {

    y = 1.1;
    pv  = pv % mydiag; // set off-diag to 0 with Schur/Hadamard multiplication
    set_proposals(pv * y, corr);

  }

}


