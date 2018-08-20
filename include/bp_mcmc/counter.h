#ifndef GUARD_bp_mcmc_counter_h
#define GUARD_bp_mcmc_counter_h


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
    Counter() : accept_ct(0), iter_ct_since_reset(0), iter_ct(0) { };

    // Add to acceptance count
    void addaccept() { ++accept_ct; ++iter_ct_since_reset; ++iter_ct; };
    // Add to iters but not accept count
    void addreject() { ++iter_ct; ++iter_ct_since_reset; };

    double getratio() const {
      if (iter_ct_since_reset > 0) {
        return (double)accept_ct / (double)iter_ct_since_reset;
      } else {
        return 0.0;
      }
    };

    void resetratio() { accept_ct = 0; iter_ct_since_reset = 0; };
    int getaccept() const { return accept_ct; };
    int getiter() const { return iter_ct; };
    int getiter_since_reset() const { return iter_ct_since_reset; };

  private:
    int accept_ct;           //  acceptance count
    int iter_ct_since_reset; //  iteration count
    int iter_ct;             //  iteration count
};

#endif

