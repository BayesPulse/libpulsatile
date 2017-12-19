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
    ProposalVariance();
    ProposalVariance(double inpv,
                     int adjust_iter,   // adjust pv on multiples of adjust_iter
                     int max_iters,     // maximum iteration to adjust pv
                     double target_ratio) {
      pv           = inpv;
      adjust_iter  = adjust_iter;
      max_iters    = max_iters;
      target_ratio = target_ratio;
    }
    double getpv() const { return pv; };

    // Counter implementation -- works, but keep an eye out for a better option
    void addreject() { count.addreject(); };
    void addaccept() { count.addaccept(); };
    double getratio() { return count.getratio(); };
    void resetratio() { count.resetratio(); };
    int getiter() { return count.getiter(); };
    int getaccept() { return count.getaccept(); };

    // 

  private:
    double pv;
    Counter count;
    int adjust_iter;  // iteration to adjust on
    int max_iter;     // iteration to stop adjusting
    double target_ratio; // target proposal variance
};

ProposalVariance::ProposalVariance()
  : pv(0)
  , count()
  , adjust_iter(500)
  , max_iter(25000)
  , target_ratio(0.35)

{
}


#endif
