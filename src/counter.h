#ifndef GUARD_counter_h
#define GUARD_counter_h

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
    Counter() : accept_ct(0), iter_ct(0) { };
    void addaccept() { ++accept_ct; ++iter_ct; }; // Add to acceptance count
    void addreject() { ++iter_ct; };              // Add to iters but not accept count
    double getratio() const { return (double) accept_ct / iter_ct; };
    void resetratio() { accept_ct = 0; iter_ct = 0; };
    int getaccept() const { return accept_ct; };
    int getiter() const { return iter_ct; };

  private:
    int accept_ct; // acceptance count
    int iter_ct;   // iteration count
};

#endif

