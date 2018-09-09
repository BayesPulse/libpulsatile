#ifndef GUARD_bp_birthdeath_deathrate_h
#define GUARD_bp_birthdeath_deathrate_h


//
// deathrate.h
//   Deathrate class for use in proposal variance classes -- only counts iters
//   since last reset, requires argument for master iteration count (this allows
//   for treating one loop over random effect parameters (pulse time, mass,
//   width) as one iteration, rather than num_pulses iterations)
//
// Author: Matt Mulvahill
// Created: 10/13/17
// Notes:
//

class Counter {

  public:
    Counter() : accept_ct(0), iter_ct_since_reset(0) { };

    // Add to acceptance count
    void addaccept() { ++accept_ct; ++iter_ct_since_reset; };
    // Add to iters but not accept count
    void addreject() { ++iter_ct_since_reset; };

    double getratio() const {
      if (iter_ct_since_reset > 0) {
        return (double)accept_ct / (double)iter_ct_since_reset;
      } else {
        return 0.0;
      }
    };

    void resetratio() { accept_ct = 0; iter_ct_since_reset = 0; };
    int getaccept() const { return accept_ct; };
    int getiter_since_reset() const { return iter_ct_since_reset; };

  private:
    int accept_ct;           //  acceptance count
    int iter_ct_since_reset; //  iteration count since last adjustment

};

#endif

