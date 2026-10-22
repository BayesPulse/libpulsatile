#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <random>
#include <bp_datastructures/patient.h>
#include <bp_datastructures/utils.h>
#include <bp_mcmc/rngscope.h>
#include <bp_mcmc/utils.h>
#include <catch.h>


//
// likelihood_cache_tests.cpp
//     Tests of the running-sum likelihood (Patient::sync_pulse_sum() and the
//     two-slot per-pulse contribution cache) against the original full
//     recomputation, and of BPRNGScope.
//


// Compare every fast-path quantity with its full-recomputation reference
static void require_matches_reference(Patient &pat) {

  double reference = pat.likelihood_reference(false);
  REQUIRE(pat.likelihood(false) ==
          Approx(reference).epsilon(1e-12).margin(1e-9));

  arma::vec partials = pat.get_partial_likelihood(false);
  REQUIRE(partials.n_elem == (arma::uword)pat.get_pulsecount());
  int i = 0;
  for (PulseIter it = pat.pulses.begin(); it != pat.pulses.end(); ++it, ++i) {
    REQUIRE(partials(i) ==
            Approx(pat.likelihood(false, it)).epsilon(1e-12).margin(1e-9));
  }

  PulseIter none;
  arma::vec mean_fast = pat.mean_concentration(false);
  arma::vec mean_full = pat.mean_concentration(false, none);
  REQUIRE(arma::approx_equal(mean_fast, mean_full, "reldiff", 1e-12));

  arma::vec diff = pat.data.concentration - mean_full;
  REQUIRE(pat.get_sumerrorsquared(false) ==
          Approx(arma::accu(diff % diff)).epsilon(1e-12));

}


TEST_CASE( "Running-sum likelihood matches the full recomputation",
           "[patient][likelihood_cache]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  utils.add_default_pulses(&pat);
  require_matches_reference(pat);

  SECTION( "through a random sequence of sampler-like updates" ) {

    std::mt19937 gen(20251022);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    auto pick = [&](int n) { return (int)(unif(gen) * n) % n; };

    for (int step = 0; step < 400; step++) {

      int op = pick(7);
      PulseIter pulse = pat.pulses.begin();
      std::advance(pulse, pick(pat.get_pulsecount()));

      if (op == 0) {
        // Propose a new location, evaluate, restore (rejected proposal)
        double current = pulse->time;
        pulse->time = current + 30 * (unif(gen) - 0.5);
        require_matches_reference(pat);
        pulse->time = current;
      } else if (op == 1) {
        // Propose a new mass and keep it (accepted proposal)
        double current = pulse->mass;
        pulse->mass = current * (0.5 + unif(gen));
        pat.likelihood(false);
        pulse->mass = current;
        pulse->mass = current * (0.5 + unif(gen));
      } else if (op == 2) {
        // Propose a new width, evaluate, restore
        double current = pulse->width;
        pulse->width = current * (0.5 + unif(gen));
        require_matches_reference(pat);
        pulse->width = current;
      } else if (op == 3) {
        // Baseline/half-life proposal (changes every pulse's contribution)
        arma::vec current = pat.estimates.baseline_halflife;
        double scale_b = 0.9 + 0.2 * unif(gen);
        double scale_h = 0.9 + 0.2 * unif(gen);
        pat.estimates.baseline_halflife = { current(0) * scale_b,
                                            current(1) * scale_h };
        require_matches_reference(pat);
        if (unif(gen) < 0.5) pat.estimates.baseline_halflife = current;
      } else if (op == 4 && pat.get_pulsecount() < 30) {
        // Birth
        PulseEstimates born(1440 * unif(gen), 1 + 4 * unif(gen),
                            10 + 60 * unif(gen), 1, 1,
                            pat.estimates.get_decay(), pat.data.time);
        pat.pulses.push_back(born);
      } else if (op == 5 && pat.get_pulsecount() > 2) {
        // Death
        pat.pulses.erase(pulse);
      } else {
        // Error variance draw
        pat.estimates.errorsq = 0.001 + 0.05 * unif(gen);
      }

      require_matches_reference(pat);

    }

  }

  SECTION( "restoring a rejected proposal restores the likelihood exactly" ) {

    double like0 = pat.likelihood(false);
    arma::vec mean0 = pat.mean_concentration(false);

    PulseIter pulse = std::next(pat.pulses.begin(), 3);
    double current = pulse->time;
    pulse->time = current + 7.5;
    double proposed = pat.likelihood(false);
    pulse->time = current;

    REQUIRE(proposed != like0);
    REQUIRE(pat.likelihood(false) == like0);
    REQUIRE(arma::all(pat.mean_concentration(false) == mean0));

  }

  SECTION( "moving a pulse to a third value is not mistaken for an undo" ) {

    PulseIter pulse = std::next(pat.pulses.begin(), 5);
    double v0 = pulse->width;
    pulse->width = v0 * 1.5;
    pat.likelihood(false);
    pulse->width = v0 * 0.7;
    require_matches_reference(pat);
    pulse->width = v0 * 1.2;
    require_matches_reference(pat);
    pulse->width = v0;
    require_matches_reference(pat);

  }

  SECTION( "a copied Patient keeps an independent, correct running sum" ) {

    double like0 = pat.likelihood(false);
    Patient copy = pat;
    copy.pulses.front().mass *= 2;
    copy.pulses.pop_back();
    require_matches_reference(copy);
    REQUIRE(pat.likelihood(false) == like0);
    require_matches_reference(pat);

  }

}


TEST_CASE( "Pulse contribution cache keeps the committed and proposed slots",
           "[pulseestimates][likelihood_cache]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  utils.add_default_pulses(&pat);
  pat.likelihood(false);  // commit every pulse to the running sum

  PulseEstimates &pulse = *std::next(pat.pulses.begin(), 2);
  double decay = pat.estimates.get_decay();
  int slot0 = pulse.current_slot();
  unsigned long generation0 = pulse.slot_generation(slot0);
  arma::vec contrib0 = pulse.mean_contribution_ref(pat.data.time, decay);

  SECTION( "restoring a proposal reuses the cached slot without recomputing" ) {
    double current = pulse.time;
    pulse.time = current + 12;
    pulse.mean_contribution_ref(pat.data.time, decay);
    REQUIRE(pulse.current_slot() != slot0);
    pulse.time = current;
    arma::vec restored = pulse.mean_contribution_ref(pat.data.time, decay);
    REQUIRE(pulse.current_slot() == slot0);
    REQUIRE(pulse.slot_generation(slot0) == generation0);
    REQUIRE(arma::all(restored == contrib0));
  }

  SECTION( "entries before the first nonzero index are exactly zero" ) {
    PulseEstimates late(1300, 3, 30, 1, 1, decay, pat.data.time);
    pat.pulses.push_back(late);
    pat.likelihood(false);
    const PulseEstimates &added = pat.pulses.back();
    arma::uword first = added.committed_first_nonzero();
    const arma::vec &contrib = added.committed_contribution();
    REQUIRE(first > 0);
    REQUIRE(first < contrib.n_elem);
    REQUIRE(arma::all(contrib.head(first) == 0));
    REQUIRE(contrib(first) > 0);
  }

}


TEST_CASE( "Initial pulse contribution is evaluated at the observation times",
           "[patient]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  double decay = pat.estimates.get_decay();

  PulseEstimates expected(pat.data.fitstart, 1, 1, 1, 1, decay, pat.data.time);
  arma::vec first = pat.pulses.front().get_mean_contribution(pat.data.time,
                                                             decay);
  REQUIRE(arma::approx_equal(first,
                             expected.get_mean_contribution(pat.data.time,
                                                            decay),
                             "absdiff", 0.0));

}


TEST_CASE( "BPRNGScope draws one continuous stream", "[rngscope]" ) {

  PulseUtils pu;
  double a, b, c;

  // Draws inside nested scopes and between them...
  pu.set_seed(424242);
  {
    BPRNGScope outer;
    { BPRNGScope inner; a = Rf_runif(0, 1); }
    b = Rf_runif(0, 1);
    { BPRNGScope inner; c = Rf_runif(0, 1); }
  }

  // ...match three consecutive draws: nothing is rewound or repeated
  pu.set_seed(424242);
  double a2, b2, c2;
  {
    BPRNGScope scope;
    a2 = Rf_runif(0, 1);
    b2 = Rf_runif(0, 1);
    c2 = Rf_runif(0, 1);
  }

  REQUIRE(a == a2);
  REQUIRE(b == b2);
  REQUIRE(c == c2);
  REQUIRE(b != c);
  REQUIRE(bp_rng_scope_depth() == 0);

}
