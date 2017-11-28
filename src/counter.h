#ifndef GUARD_proposalvariance_h
#define GUARD_proposalvariance_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

//
// counter.h
//   Counter class for use in proposal variance classes
//
// Author: Matt Mulvahill
// Created: 10/13/17
// Notes:
//

class Counter {
  public:
    Counter(double target_ratio, int adjust_at_iter, int max_iters);
    void addaccept() { ++accept_ct; ++iter_ct; }; // Add to acceptance count
    void addreject() { ++iter_ct; };              // Add to iters but not accept count
    double getratio() const { return (double) accept_ct / iter_ct; };
    void resetratio() { accept_ct = 0; iter_ct = 0; };
    int getaccept() const { return accept_ct; };
    int getiter() const { return iter_ct; };

  private:
    int accept_ct;       // acceptance count
    int iter_ct;         // iteration count
    int adjust_at_iter;  // iteration to adjust on
    int max_iters;       // iteration to stop adjusting
    double target_ratio; // target proposal variance
};

Counter::Counter(double target_ratio, int adjust_at_iter, int max_iters)
{

  adjust_at_iter = adjust_at_iter;
  max_iters      = max_iters;
  accept_ct      = 0;
  iter_ct        = 0;
  target_ratio   = target_ratio;

}


#endif

