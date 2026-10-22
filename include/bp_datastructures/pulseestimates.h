#ifndef GUARD_bp_datastructures_pulseestimates_h
#define GUARD_bp_datastructures_pulseestimates_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
//#include <bp_mcmc/utils.h>

//
// pulseestimates.h
//   defining the object for holding pulse level estimates
//
// Author: Matt Mulvahill
// Notes:
//

using namespace Rcpp;


//
// PulseEstimates structure
//   aka PulseParms -- pulse chain updated by mcmc algorithm called
//   TODO: Check  on lambda, does fsh need separate definition?
//

class PulseEstimates {

  public:

    // Variables used in all models
    double time;
    double mass;
    double width;
    double tvarscale_mass;   // variance scale for mass t-dist (eta)
    double tvarscale_width;  // variance scale for width t-dist (eta)

    // For joint hormone model (response pulses only)
    double lambda;  // Coupling intensity from driver hormone pulses

    // Public user-facing access function for mean_contrib. Only calculates
    // mean_contrib if input values have changed. Returns a copy; hot paths use
    // mean_contribution_ref() to avoid it.
    arma::vec get_mean_contribution(const arma::vec & data_time, double decay_rate)
    {
      return mean_contribution_ref(data_time, decay_rate);
    }

    // mean_contribution_ref()
    //   Two-slot cache of this pulse's contribution to the mean concentration,
    //   keyed on (time, mass, width, decay_rate). MH samplers set a proposed
    //   value, evaluate the likelihood, then restore the current value; with
    //   two slots the restored value is a cache hit, not a recomputation. The
    //   slot committed to a Patient's running sum (Patient::sync_pulse_sum())
    //   is never overwritten, so the Patient can always subtract it.
    const arma::vec & mean_contribution_ref(const arma::vec & data_time,
                                            double decay_rate)
    {
      if (slot_matches(cur_slot_, decay_rate)) return contrib_[cur_slot_];

      int other = 1 - cur_slot_;
      if (slot_matches(other, decay_rate)) {
        cur_slot_ = other;
        return contrib_[cur_slot_];
      }

      // Miss: compute into a slot that is not committed to a running sum
      int target = (cur_slot_ == committed_slot_) ? other : cur_slot_;
      first_nonzero_[target] =
        calc_mean_contribution(data_time, decay_rate, contrib_[target]);
      slot_generation_[target]++;
      key_[target][0] = time;
      key_[target][1] = mass;
      key_[target][2] = width;
      key_[target][3] = decay_rate;
      slot_valid_[target] = true;
      cur_slot_ = target;
      return contrib_[cur_slot_];
    }

    // Running-sum bookkeeping used by Patient (not for samplers)
    const void * sum_owner() const { return sum_owner_; }
    const arma::vec & committed_contribution() const {
      return contrib_[committed_slot_];
    }
    // Every entry of the committed contribution before this index is 0
    arma::uword committed_first_nonzero() const {
      return first_nonzero_[committed_slot_];
    }
    int committed_slot() const { return committed_slot_; }
    int current_slot() const { return cur_slot_; }
    // Changes every time a slot's contents are recomputed
    unsigned long slot_generation(int slot) const {
      return slot_generation_[slot];
    }
    void commit_to(const void * owner) {
      sum_owner_      = owner;
      committed_slot_ = cur_slot_;
    }
    // True when the current parameters (and decay) match the committed slot
    bool committed_matches(double decay_rate) const {
      return committed_slot_ >= 0 && slot_matches(committed_slot_, decay_rate);
    }

    // For returning to chain (single hormone model)
    arma::rowvec get_vector_of_values() const
    {
      arma::rowvec out { time, mass, width, tvarscale_mass, tvarscale_width };
      return out;
    }

    // For returning to chain (joint hormone model - includes lambda)
    arma::rowvec get_vector_of_values_with_lambda() const
    {
      arma::rowvec out { time, mass, width, tvarscale_mass, tvarscale_width, lambda };
      return out;
    }

    // Constructor for new PulseEstimates objects
    PulseEstimates(double in_time,
                   double in_mass,
                   double in_width,
                   double in_tvarscale_mass,
                   double in_tvarscale_width,
                   double patient_decay,
                   const arma::vec &data_time,
                   double in_lambda = 0.0)  // lambda defaults to 0 for single hormone
        : time            (in_time)
        , mass            (in_mass)
        , width           (in_width)
        , tvarscale_mass  (in_tvarscale_mass)
        , tvarscale_width (in_tvarscale_width)
        , lambda          (in_lambda)
      {
        mean_contribution_ref(data_time, patient_decay);
      }

    // Constructor for empty pulse object
    PulseEstimates()
      : time(0), mass(0), width(0), tvarscale_mass(0),
        tvarscale_width(0), lambda(0.0)
      {
        contrib_[0].zeros(1);
      }

  private:

    arma::vec contrib_[2];
    arma::uword first_nonzero_[2] = {0, 0};      // entries before are 0
    unsigned long slot_generation_[2] = {0, 0};  // bumped on recompute
    double key_[2][4];
    bool slot_valid_[2]     = {false, false};
    int cur_slot_           = 0;
    int committed_slot_     = -1;      // slot included in owner's running sum
    const void * sum_owner_ = nullptr; // Patient whose running sum holds it

    bool slot_matches(int slot, double decay_rate) const {
      return slot_valid_[slot] && key_[slot][0] == time &&
        key_[slot][1] == mass && key_[slot][2] == width &&
        key_[slot][3] == decay_rate;
    }

    // mean_contribution() of each pulse to the total mean_concentration
    // Returns the index of the first nonzero entry (n_elem if all are zero).
    arma::uword calc_mean_contribution(const arma::vec &data_time,
                                       double decay_rate,
                                       arma::vec &mean_contribution)
    {

      double y, z, w;

      z  = width * decay_rate;
      y  = decay_rate * (0.5 * z  + time);
      z += time;
      w  = sqrt(2. * width);
      const double sqrt2 = sqrt(2.);

      // Single pass, no temporaries: contribution(t) =
      //   mass * pnorm(((t - z) / w) * sqrt(2)) * exp(y - t * decay_rate),
      // with the same operation order as the former vectorized expression, so
      // results are bit-for-bit unchanged. Where pnorm() underflows to 0 the
      // entry is 0 after the clamp below regardless of exp(), so exp() is
      // skipped there.
      const arma::uword N = data_time.n_elem;
      mean_contribution.set_size(N);
      const double * t = data_time.memptr();
      double * out     = mean_contribution.memptr();
      arma::uword first_nonzero = N;
      for (arma::uword i = 0; i < N; i++) {
        // R's pnorm() returns exactly 0 for x <= -37.5193 (lower tail), so it
        // is skipped below -37.52; x == NaN fails the test and is passed on.
        double x = ((t[i] - z) / w) * sqrt2;
        double p = (x < -37.52) ? 0. : Rf_pnorm5(x, 0.0, 1.0, 1, 0);
        double val = (p == 0.) ? 0. : (mass * p) * exp(y - t[i] * decay_rate);
        // Truncate <0 = 0, and sanitize non-finite entries to 0.
        //
        // A pulse's contribution is a bounded, non-negative secretion value.
        // For an extreme (near-degenerate) pulse -- e.g. a width so large the
        // random-walk sampler has let it diffuse far past the observation
        // window -- the factors above can overflow: pnorm(x) underflows to 0
        // while exp(y - t*decay) grows to +Inf, giving 0 * Inf = NaN (or a
        // spurious +Inf). Physically the true contribution of such a pulse at
        // any finite time is ~0, so any non-finite entry is a floating-point
        // artifact and is clamped to 0. Without this a single pathological
        // width poisons mean_concentration -> the error-variance Gibbs draw's
        // sum-of-squares -> errorsq = NaN, which then breaks downstream coda
        // diagnostics (na.fail). Finite, in-range pulses are unaffected.
        out[i] = (std::isfinite(val) && val > 0.) ? val : 0.;
        if (out[i] != 0. && first_nonzero == N) first_nonzero = i;
      }
      return first_nonzero;

    }


};



#endif
