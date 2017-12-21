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

class ProposalVariance {

  public:
    ProposalVariance()
      : pv(0)
      , count()
      , adjust_iter(500)
      , max_iter(25000)
      , target_ratio(0.35) { };
    ProposalVariance(double in_pv,
                     int in_adjust_iter,   // adjust pv on multiples of adjust_iter
                     int in_max_iter,     // maximum iteration to adjust pv
                     double in_target_ratio) {
      adjust_iter  = in_adjust_iter;
      max_iter    = in_max_iter;
      target_ratio = in_target_ratio;
      initialize_proposals(in_pv);
    }

    // ProposalVariance functions
    double getpv() const { return pv; };
    void adjustpv();

    // Counter implementation -- works, but keep an eye out for a better option
    void addreject() { count.addreject(); };
    void addaccept() { count.addaccept(); };
    double getratio() { return count.getratio(); };
    void resetratio() { count.resetratio(); };
    int getiter() { return count.getiter(); };
    int getaccept() { return count.getaccept(); };

  private:
    double pv;
    double psd;          // proposal standard deviation
    Counter count;
    int adjust_iter;  // iteration to adjust on
    int max_iter;     // iteration to stop adjusting
    double target_ratio; // target proposal variance

    // ProposalVariance internal functions
    void initialize_proposals(double initial_pv);
    void set_proposals(double this_pv);

};


class ProposalVariance2p {

  public:
    ProposalVariance2p();
    ProposalVariance2p(arma::vec in_pv,
                           int in_adjust_iter,   // adjust pv on multiples of adjust_iter
                           int in_max_iter,      // maximum iteration to adjust pv
                           double in_target_ratio);

    // ProposalVariance functions
    arma::mat getpv() const  { return pv; };
    arma::mat getpsd() const { return psd; };
    void adjustpv(double corr);

    // Counter implementation -- works, but keep an eye out for a better option
    void addreject() { count.addreject(); };
    void addaccept() { count.addaccept(); };
    double getratio() { return count.getratio(); };
    void resetratio() { count.resetratio(); };
    int getiter() { return count.getiter(); };
    int getaccept() { return count.getaccept(); };

  private:
    arma::mat::fixed<2, 2> pv;
    arma::mat::fixed<2, 2> psd; // proposal standard deviation
    Counter count;
    int adjust_iter;  // iteration to adjust on
    int max_iter;     // iteration to stop adjusting
    double target_ratio; // target proposal variance

    // ProposalVariance internal functions
    void initialize_proposals(arma::vec initial_pv);
    void set_proposals(arma::mat this_pv, double corr);
    double calc_covariance(arma::mat pv, double corr);

};



#endif
