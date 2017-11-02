#ifndef GUARD_proposalvariance_h
#define GUARD_proposalvariance_h

#include <Rcpp.h>
#include <RcppArmadillo.h>

//
// proposalvariance.h
//   Definitions for the ProposalVariance class
//
// Author: Matt Mulvahill
// Created: 10/13/17
//
// NOTE: Sec9.3.1 in Acc C++ define here to tell compiler to avoid function-call
// overhead, so all simple fn's are defined here -- not sure if this works for
// non-const functions.
//

class ProposalVariance {

  public:
    // Constructor
    ProposalVariance(double initial_pv, // proposal variance
                     int adjust_iter,   // adjust pv on multiples of adjust_iter
                     int max_iters);    // maximum iteration to adjust pv
    ~ProposalVariance(); // Destructor
    virtual adjustpv();
    void addaccept() { ++accept_ct; ++iter_ct; };            // Add to acceptance count
    void addreject() { ++iter_ct; };                         // Add to iters but not accept count
    double getratio() const { return accept_ct / iter_ct; };
    void resetratio() { accept_ct = 0; iter_ct = 0; };
    virtual getpv() const { return pv; };

  private:
    int accept_ct; // acceptance count
    int iter_ct;   // iteration count
    virtual pv;     // proposal variance
    int adjust_iter; // iteration to adjust on
    int max_iter;

};

class MarginalProposalVariance : public ProposalVariance
{
  public:
    void adjustpv(double target_pv = 0.25);
    double getpv const { return pv; };
  private:
    double pv; // proposal var-covar matrix
}

class JointProposalVariance : public ProposalVariance
{
  public:
    void adjustpv(double corr = -0.90, double target_pv = 0.25);
    arma::mat getpv() const { return pv; };
  private:
    arma::mat pv(2, 2); // proposal var-covar matrix
}

#endif

