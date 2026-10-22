#ifndef GUARD_patient_h
#define GUARD_patient_h

#include <RcppArmadillo.h>
#include <utility>
#include <vector>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include "bp_datastructures/patientdata.h"
#include "bp_datastructures/patientestimates.h"
#include "bp_datastructures/patientpriors.h"
#include "bp_datastructures/pulseestimates.h"



//
// patient.h
//   defining the patient class and subclasses
//
//   NOTE: getting too big -- clean up/reorganize
//


//
// Parent class of Patients and Populations, object that MMH and Gibbs
// recognize
//
//struct MCMCSamplingUnit { };


using namespace Rcpp;

typedef std::list<PulseEstimates> PulseList;
typedef PulseList::iterator PulseIter;
typedef PulseList::const_iterator PulseConstIter;


//
// Patient struct
//
struct Patient {

  // Member objects for all models
  PatientData data;
  PatientEstimates estimates;
  PulseList pulses;
  PulseIter piter = pulses.begin();
  PulseList responses;
  PulseIter riter = responses.begin();

  // Random-effects distribution for pulse mass/width. When false (default),
  // pulse random effects are Student-t via the per-pulse t-scale (tvarscale,
  // kappa) scale-mixture. When true, kappa is held at 1 for every pulse and the
  // t-scale is never drawn, giving Gaussian pulse random effects. The
  // birth-death process reads this to initialize new pulses, and the MCMC loop
  // reads it to skip the t-scale draw.
  bool gaussian_random_effects = false;

  // Scale of the pulse mass/width random effects. When true (the published
  // Mulvahill-thesis / Horton parameterization), the random effect is modeled on
  // the LOG scale: log(mass), log(width) ~ N(mean, sd^2 / kappa), so mass_mean /
  // width_mean are means of the log values, mass_sd / width_sd are SDs of the log
  // values, positivity is automatic (no truncation at 0), and no truncation
  // normalizing constant appears. When false (the current natural-scale research
  // parameterization), the random effect is a truncated normal on the natural
  // scale. The pulse `mass`/`width` members always store the natural-scale value
  // used by the secretion function; only the prior/draw geometry changes.
  // C++ default is false to keep existing standalone behavior/tests stable; the R
  // model specs default this to true (log-normal) to match the papers.
  bool lognormal_pulses = false;

  // Prior on the pulse-to-pulse SD of the random effects (mass_sd / width_sd).
  // When false (default), the SD has the current half-Cauchy prior with scale
  // mass_sd_param / width_sd_param, contributing a non-trivial prior ratio to the
  // SD-of-random-effects MH acceptance ratio. When true (the published Mulvahill-
  // thesis / Horton parameterization), the SD has a Uniform(0, mass_sd_max) /
  // Uniform(0, width_sd_max) prior: the density is constant on its support, so the
  // prior ratio cancels to 0 and support is enforced as 0 < sigma < max. This axis
  // is orthogonal to lognormal_pulses (it only affects the prior ratio and the
  // support bound, not the log-vs-natural residual geometry). C++ default is false
  // (half-Cauchy) to keep existing standalone behavior/tests stable; the R model
  // specs default this to true (uniform) to match the papers.
  bool uniform_sd_prior = false;

  //
  // For single-subject model
  //

  // Member objects
  PatientPriors priors;

  // Constructor
  Patient(PatientData in_data,
          PatientPriors in_priors,
          PatientEstimates in_parms) :
      data(in_data)
    , estimates(in_parms)
    , priors(in_priors) {

    PulseEstimates firstpulse(in_data.fitstart,
                              1, 1, 1, 1,
                              in_parms.get_decay(), 
                              in_data.time);
    pulses.push_back(firstpulse);

  }


  //
  // For population models
  //

  // Constructor
  Patient(PatientData in_data,
          PatientEstimates in_parms) :
      data(in_data)
    , estimates(in_parms) 
    , priors() {
      // priors member obj not used
  }


  //
  // Methods for both types/all-models
  //

  // get_pulsecount()
  //   Get current number of pulses
  int get_pulsecount() const { return pulses.size(); };

  // get_sumerrorsquared()
  //   Sums of squared error for error gibbs
  double get_sumerrorsquared(bool response_hormone) {

    arma::vec mean = mean_concentration(response_hormone);
    arma::vec conc;

    if (response_hormone) {
      conc = data.response_concentration;
    } else {
      conc = data.concentration;
    }

    arma::vec diff = conc - mean;
    arma::vec squareddiff =  diff % diff; // % Schur product (elementwise multiplication)
    double ssq = arma::accu(squareddiff);

    //Rcpp::Rcout << "Conc vector:\n" << conc << std::endl;
    //Rcpp::Rcout << "Mean vector:\n" << mean << std::endl;
    //Rcpp::Rcout << "Diff vector:\n" << diff << std::endl;
    //Rcpp::Rcout << "Square diff vector:\n" << squareddiff << std::endl;
    //Rcpp::Rcout << "Sums of squares:" << ssq << std::endl;

    return ssq;

  }




  // likelihood()
  //   computes the current log-likelihood using the observed log-concentrations
  //   and mean concentration. Not stored. For the driver/single hormone
  //   (response_hormone == false) this uses the running sum of pulse
  //   contributions (O(n) per call instead of O(K * n)); see sync_pulse_sum().
  //   The pulse-excluding overload and likelihood_reference() keep the original
  //   full recomputation and are used to verify the fast path.
  double likelihood(bool response_hormone) {
    if (response_hormone) return likelihood_reference(response_hormone);
    sync_pulse_sum();
    double like;
    const double baseline = estimates.baseline_halflife(0);
    if (memo_.valid && memo_.baseline == baseline &&
        memo_.errorsq == estimates.errorsq) {
      like = memo_.like;
    } else {
      like = loglik_from_sum();
      memo_.valid    = true;
      memo_.like     = like;
      memo_.baseline = baseline;
      memo_.errorsq  = estimates.errorsq;
    }
#ifdef BP_VERIFY_LIKELIHOOD
    verify_against_reference(like, likelihood_reference(false), "likelihood");
#endif
    return like;
  }

  double likelihood_reference(bool response_hormone) {
    PulseIter emptyiter;
    return likelihood(response_hormone, emptyiter);
  }

  double likelihood(bool response_hormone, PulseIter pulse_excluded) {

    double like = 0;
    arma::vec conc;
    arma::vec diffs;

    if (response_hormone) {
      conc = data.response_concentration;
    } else {
      conc = data.concentration;
    }

    // Calculate likelihood
    diffs = conc - mean_concentration(response_hormone, pulse_excluded);
    diffs = arma::square(diffs);
    like  = arma::accu(diffs);
    like /= (-2.0 * estimates.errorsq);
    like += -0.5 * conc.n_elem * (1.8378771 + estimates.get_logerrorsq());

    return like;

  }

  // Calculate all partial likelihoods (excluding each pulse(i))
  //   Each partial is computed from the running sum minus that pulse's
  //   contribution: O(K * n) in total rather than O(K^2 * n).
  arma::vec get_partial_likelihood(bool response_hormone) {

    arma::vec partials(get_pulsecount());

    if (response_hormone) {
      PulseIter exclude_pulse = responses.begin();
      int i = 0;
      while (exclude_pulse != responses.end()) {
        partials(i) = likelihood(response_hormone, exclude_pulse);
        exclude_pulse++;
        i++;
      }
      return partials;
    }

    sync_pulse_sum();

    // A pulse contributes exactly 0 before its onset, so there the partial
    // residuals equal the full ones: precompute prefix sums of the full
    // squared residuals and only evaluate each partial from the pulse's first
    // nonzero contribution on.
    const arma::uword n = data.concentration.n_elem;
    const double * conc = data.concentration.memptr();
    const double * sum  = pulse_sum_.memptr();
    const double baseline = estimates.baseline_halflife(0);
    full_ssq_prefix_.set_size(n + 1);
    double * prefix = full_ssq_prefix_.memptr();
    prefix[0] = 0.0;
    for (arma::uword j = 0; j < n; j++) {
      double r = conc[j] - std::log(sum[j] + baseline);
      prefix[j + 1] = prefix[j] + r * r;
    }

    int i = 0;
    for (auto &pulse : pulses) {
      const double * c = pulse.committed_contribution().memptr();
      arma::uword j = pulse.committed_first_nonzero();
      double ssq = prefix[j];
      for (; j < n; j++) {
        double r = conc[j] - std::log(sum[j] - c[j] + baseline);
        ssq += r * r;
      }
      partials(i) = ssq_to_loglik(ssq);
      i++;
    }
#ifdef BP_VERIFY_LIKELIHOOD
    i = 0;
    for (PulseIter it = pulses.begin(); it != pulses.end(); ++it, ++i) {
      verify_against_reference(partials(i), likelihood(false, it),
                               "partial likelihood");
    }
#endif

    return partials;

  }

  // mean_concentration()
  //   log of baseline plus the sum of each pulse's mean_contrib vector.
  //   There is a version for a) using all pulses (running sum for the driver
  //   hormone) and b) one for excluding one pulse (full recomputation).
  arma::vec mean_concentration(bool response_hormone) {

    if (!response_hormone) {
      sync_pulse_sum();
      return arma::log(pulse_sum_ + estimates.baseline_halflife(0));
    }
    PulseIter emptyiter;
    return mean_concentration(response_hormone, emptyiter);

  }

  arma::vec mean_concentration(bool response_hormone, PulseIter pulse_excluded)
  {

    arma::vec mean_conc(data.concentration.n_elem);
    mean_conc.fill(0);
    PulseIter pulse_iter;
    PulseConstIter pulselist_end;

    if (response_hormone) {
      pulse_iter    = responses.begin();
      pulselist_end = responses.end();
    } else {
      pulse_iter    = pulses.begin();
      pulselist_end = pulses.end();
    }

    // Add the contribution to the mean from each pulse
    while (pulse_iter != pulselist_end) {
      if (pulse_iter != pulse_excluded) {
        mean_conc += pulse_iter->mean_contribution_ref(data.time,
                                                       estimates.get_decay());
      }
      ++pulse_iter;
    }

    // Add the baseline contribution and log
    mean_conc += estimates.baseline_halflife(0);
    mean_conc = log(mean_conc);

    return mean_conc;

  }

  // rebuild_pulse_sum()
  //   Recompute the running sum of pulse contributions from scratch, in list
  //   order (the same order as the full recomputation, so the likelihood right
  //   after a rebuild matches the original bit for bit). Called whenever the
  //   pulse set changes in a way sync_pulse_sum() cannot apply incrementally,
  //   and once per MCMC iteration by the birth-death sampler to bound
  //   floating-point drift from incremental updates.
  void rebuild_pulse_sum() {

    const double decay = estimates.get_decay();
    pulse_sum_.zeros(data.concentration.n_elem);
    for (auto &pulse : pulses) {
      pulse_sum_ += pulse.mean_contribution_ref(data.time, decay);
      pulse.commit_to(this);
    }
    pulse_sum_count_ = pulses.size();
    pulse_sum_valid_ = true;
    undo_valid_      = false;
    memo_.valid      = false;
    memo_prev_.valid = false;

  }

  // sync_pulse_sum()
  //   Bring the running sum up to date with the current pulse parameters and
  //   decay. Samplers modify pulses in place (propose, evaluate, restore), so
  //   each pulse is checked against the parameters it was committed with and
  //   its change is applied as (new - committed) -- O(n) per changed pulse.
  //   The sum before the last such update is kept, so the restore after a
  //   rejected proposal swaps it back (bit-identical, with its likelihood
  //   still memoized) instead of recomputing. Pulses not yet committed to this
  //   Patient (births) are added. Removed pulses (or a copied Patient, whose
  //   pulses were committed to the original) are detected by count and
  //   trigger a full rebuild.
  void sync_pulse_sum() {

    if (!pulse_sum_valid_ ||
        pulse_sum_.n_elem != data.concentration.n_elem) {
      rebuild_pulse_sum();
      return;
    }

    const double decay = estimates.get_decay();
    int owned = 0;
    bool has_new = false;
    changed_.clear();
    for (auto &pulse : pulses) {
      if (pulse.sum_owner() == this) {
        owned++;
        if (!pulse.committed_matches(decay)) changed_.push_back(&pulse);
      } else {
        has_new = true;
      }
    }

    if (owned != pulse_sum_count_) {
      rebuild_pulse_sum();
      return;
    }

    if (!changed_.empty()) {

      for (auto pulse : changed_) {
        pulse->mean_contribution_ref(data.time, decay);
      }

      if (is_undo_of_last_update()) {

        // Exactly reverses the last update: restore the previous sum
        pulse_sum_.swap(pulse_sum_prev_);
        std::swap(memo_, memo_prev_);
        for (auto &entry : undo_) {
          // Record the slot being left, so the undo can itself be undone
          entry.slot       = entry.pulse->committed_slot();
          entry.generation = entry.pulse->slot_generation(entry.slot);
          entry.pulse->commit_to(this);
        }

      } else {

        pulse_sum_prev_ = pulse_sum_;
        memo_prev_      = memo_;
        memo_.valid     = false;
        undo_.clear();
        for (auto pulse : changed_) {
          int slot = pulse->committed_slot();
          undo_.push_back(UndoEntry{pulse, slot, pulse->slot_generation(slot)});
          pulse_sum_ -= pulse->committed_contribution();
          pulse->commit_to(this);
          pulse_sum_ += pulse->committed_contribution();
        }
        undo_valid_ = true;

      }
    }

    if (has_new) {
      for (auto &pulse : pulses) {
        if (pulse.sum_owner() != this) {
          pulse_sum_ += pulse.mean_contribution_ref(data.time, decay);
          pulse.commit_to(this);
          pulse_sum_count_++;
        }
      }
      undo_valid_      = false;
      memo_.valid      = false;
      memo_prev_.valid = false;
    }

  }

  // calc_sr_strauss()
  //   Calculates sum(S(R)), the exponent on the gamma parameter in the Strauss
  //   process/prior for pulse location. Used for Strauss prior in birth_death
  //   and mh_time_strauss.
  //   location is time/loc to test against other pulses
  int calc_sr_strauss(double location) {
    PulseIter emptyiter;
    return calc_sr_strauss(location, &(*emptyiter));
  }

  int calc_sr_strauss(double location, PulseEstimates * pulse_excluded) {

    int s_r = 0;       // Sum of indicators where diff < 20
    double difference = 0.; // Time difference
    PulseIter pulse = pulses.begin();
    PulseConstIter pulse_end = pulses.end();

    while (pulse != pulse_end) {
      if (&(*pulse) != pulse_excluded) { // TODO: Test that pulse is actually excluded!
        // skip if node is same that location is from;
        difference = abs(location - pulse->time);
        // increment by 1 if diff<R
        s_r = (difference < priors.strauss_repulsion_range) ? s_r + 1 : s_r; 
      }
      pulse++;
    }

    // sum(S(R)) - scalar value for sum of # pulses too close to each other
    return(s_r); 

  }


#ifdef BP_VERIFY_LIKELIHOOD
  // Debug builds (-DBP_VERIFY_LIKELIHOOD) check every running-sum likelihood
  // and partial likelihood against the original full recomputation and stop
  // on a mismatch beyond floating-point rounding.
  static void verify_against_reference(double fast, double reference,
                                       const char * what) {
    double tol = 1e-9 * (1.0 + std::abs(reference));
    if (!(std::abs(fast - reference) <= tol) &&
        !(std::isinf(fast) && fast == reference)) {
      Rcpp::stop("BP_VERIFY_LIKELIHOOD: %s mismatch: "
                 "fast=%.17g reference=%.17g", what, fast, reference);
    }
  }
#endif

  private:

  // Running sum of the contributions of all pulses (driver/single hormone)
  arma::vec pulse_sum_;
  bool pulse_sum_valid_ = false;
  int pulse_sum_count_  = 0;
  arma::vec full_ssq_prefix_;  // scratch for get_partial_likelihood()

  // One-step undo of the last incremental update (see sync_pulse_sum())
  struct UndoEntry {
    PulseEstimates * pulse;
    int slot;                  // committed slot before the update
    unsigned long generation;  // that slot's generation at the time
  };
  arma::vec pulse_sum_prev_;   // running sum before the last update
  std::vector<UndoEntry> undo_;
  std::vector<PulseEstimates *> changed_;  // scratch for sync_pulse_sum()
  bool undo_valid_ = false;

  // Memoized likelihood of the running sum (and of the previous sum)
  struct LikelihoodMemo {
    bool valid = false;
    double like = 0.0, baseline = 0.0, errorsq = 0.0;
  };
  LikelihoodMemo memo_, memo_prev_;

  // True when the pulses changed since the last commit are exactly those of
  // the last update, each now back at the slot it was committed to before it.
  bool is_undo_of_last_update() const {
    if (!undo_valid_ || changed_.size() != undo_.size()) return false;
    for (size_t k = 0; k < undo_.size(); k++) {
      const PulseEstimates * pulse = changed_[k];
      if (pulse != undo_[k].pulse ||
          pulse->current_slot() != undo_[k].slot ||
          pulse->slot_generation(undo_[k].slot) != undo_[k].generation) {
        return false;
      }
    }
    return true;
  }

  // loglik_from_sum()
  //   Gaussian log-likelihood of the log-concentrations given the running sum
  //   of pulse contributions. Accumulates the squared residuals the same way
  //   arma::accu() does (two interleaved accumulators), so right after a
  //   rebuild it matches likelihood_reference() bit for bit.
  double loglik_from_sum() {

    const double * conc = data.concentration.memptr();
    const double * sum  = pulse_sum_.memptr();
    const arma::uword n = data.concentration.n_elem;
    const double baseline = estimates.baseline_halflife(0);
    double acc1 = 0.0, acc2 = 0.0;
    arma::uword j;

    for (j = 1; j < n; j += 2) {
      double r1 = conc[j - 1] - std::log(sum[j - 1] + baseline);
      double r2 = conc[j] - std::log(sum[j] + baseline);
      acc1 += r1 * r1;
      acc2 += r2 * r2;
    }
    if ((j - 1) < n) {
      double r1 = conc[j - 1] - std::log(sum[j - 1] + baseline);
      acc1 += r1 * r1;
    }

    return ssq_to_loglik(acc1 + acc2);

  }

  // Gaussian log-likelihood from the sum of squared residuals
  double ssq_to_loglik(double ssq) {
    double like = ssq;
    like /= (-2.0 * estimates.errorsq);
    like += -0.5 * data.concentration.n_elem *
      (1.8378771 + estimates.get_logerrorsq());
    return like;
  }

};




#endif
