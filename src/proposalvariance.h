#ifndef GUARD_proposalvariance_h
#define GUARD_proposalvariance_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "counter.h"

//
// proposalvariance.h
//   Definitions for the ProposalVariance class
//
// Author: Matt Mulvahill
// Created: 10/13/17
//
// NOTE: Sec9.3.1 in Acc C++ define here to tell compiler to avoid function-call
// overhead, so all simple fn's are defined here -- not sure if this works for
// non-const functions or virtual functions.
//
// NOTE: ProposalVariance has a constructor only to document the requirements
// for the constructor.  It can't actually be used since its an abstract class.
//

template <typename input_type, typename internal_type>
class ProposalVariance {

  public:
    ProposalVariance(input_type initial_pv,  // proposal variance
                     int adjust_iter,   // adjust pv on multiples of adjust_iter
                     int max_iters,     // maximum iteration to adjust pv
                     double target_ratio); // target accept ratio
    ~ProposalVariance();
    void adjustpv();
    void adjustpv(double corr = -0.90);
    internal_type getpv() const { return pv; };
    internal_type getpsd() const { return psd; };
    Counter counter;

  private:
    int adjust_iter;  // iteration to adjust on
    int max_iter;     // iteration to stop adjusting
    double target_ratio; // target proposal variance
    internal_type pv;           // proposal variance
    internal_type psd;          // proposal standard deviation

    void initialize_proposals(input_type initial_pv);
    void set_proposals(double this_pv);
    void set_proposals(arma::mat this_pv, double corr);
    double calc_covariance(arma::mat pv, double corr);
//    arma::mat::fixed<2, 2> pv; // proposal var-covar matrix 
//    arma::mat::fixed<2, 2> psd;  // proposal correlation matrix

};


// -------------------------------------
// Constructor
// -------------------------------------
template <typename input_type, typename internal_type>
ProposalVariance<input_type, internal_type>::ProposalVariance(input_type initial_pv,
                                                              int adjust_at_iter,
                                                              int max_iters,
                                                              double target_ratio) // default 0.35 for 1, 0.25 for 2-parameter version
  : adjust_iter = adjust_at_iter
  , max_iters   = max_iters
  , target_ratio = target_ratio
  , counter()
{

  initialize_proposals(initial_pv);

}


// -------------------------------------
// initialize_proposals() for each type
// -------------------------------------
// One for one-variable and one for two-variable MMH/pv.
// create/set pv and psd objects (proposal var-covar and correlation matrices)
void ProposalVariance::initialize_proposals(double this_pv)
{

  set_proposals(this_pv);

}

void ProposalVariance::initialize_proposals(arma::vec initial_pv)
{

  arma::mat this_pv = arma::diagmat(initial_pv);
  set_proposals(this_pv, -0.90);

}

// -------------------------------------
// set_proposals() for each type
// -------------------------------------
void ProposalVariance::set_proposals(double this_pv)
{

  pv = this_pv;
  psd = sqrt(pv);

}

void ProposalVariance::set_proposals(arma::mat this_pv, double corr)
{

  pv = this_pv;
  pv(0, 1) = pv(1, 0) = calc_covariance(pv, corr);
  psd = chol(pv);

}

double ProposalVariance::calc_covariance(arma::mat pv, double corr)
{

  return corr * sqrt(pv(0, 0) * pv(1, 1));

}


// -------------------------------------
// adjustpv() for each type
// -------------------------------------
// Acceptance adjustment routine for one-parameter modified MH
void ProposalVariance::adjustpv()
{

  double y = 1.0 + 1000.0 * pow(getratio() - target_ratio, 3);
  if (y < 0.9)      set_proposals(pv * 0.9);
  else if (y > 1.1) set_proposals(pv * 1.1);

}

// Acceptance adjustment routine for two-parameter modified MH used previously
// for baseline and halflife
//   - args: corr correlation between the two variances
void ProposalVariance::adjustpv(double corr = -0.90)
{

  arma::mat mydiag(2, 2, arma::fill::eye);

  // y - new diagonal elements of proposal variance-covariance matrix based on
  // inputs
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





#endif


//////////////////////////////////////////////////////////////////////////

//class MarginalProposalVariance : public ProposalVariance
//{
//
//  public:
//    MarginalProposalVariance(double initial_pv, int adjust_iter, int max_iters, double target_ratio);
//    ~ProposalVariance(); // Destructor
//
//    void adjustpv();
//    double getpv  const { return pv; };
//    double getpsd const { return psd; };
//
//  private:
//    void initialize_proposals(double this_pv);
//    void set_proposals(double this_pv);
//    double pv; // proposal var-covar matrix
//    double psd; // proposal var-covar matrix
//    double target_ratio = 0.35;
//
//};
//
//
//class JointProposalVariance : public ProposalVariance
//{
//
//  public:
//    JointProposalVariance(arma::vec initial_pv, int adjust_iter, int max_iters, double target_ratio);
//    ~ProposalVariance(); // Destructor
//
//    void adjustpv(double corr = -0.90);
//    arma::mat getpv()  const { return pv; };
//    arma::mat getpsd() const { return psd; };
//
//  private:
//    void initialize_pv();
//    void update_pv();
//    double set_covariance(double var00, double var11, double corr);
//    double set_covariance(double var00, double var11);
//    arma::mat::fixed<2, 2> pv; // proposal var-covar matrix 
//    arma::mat::fixed<2, 2> psd;  // proposal correlation matrix
//    double target_ratio = 0.25;
//
//};
