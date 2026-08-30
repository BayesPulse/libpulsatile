#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/population.h>
#include <bp_datastructures/patient.h>
#include <bp_datastructures/patientdata.h>
#include <bp_datastructures/patientpriors.h>
#include <bp_datastructures/patientestimates.h>
#include <bpmod_population/pop_draw_means.h>
#include <bpmod_population/pop_draw_baseline_halflife_means.h>
#include <bpmod_population/pop_draw_mean_sds.h>
#include <bpmod_population/pop_draw_baseline_halflife_sds.h>
#include <bpmod_population/pop_draw_pulse_sds.h>
#include <bpmod_population/pop_draw_error.h>
#include <testing/catch.h>

//
// population_tests.cpp
//     Tests for population model data structures and samplers
//

//
// PopulationEstimates tests
//

TEST_CASE( "PopulationEstimates constructor works", "[population]" ) {

  PopulationEstimates pop_est(
    3.5,    // mass_mean
    30.0,   // width_mean
    2.5,    // baseline_mean
    45.0,   // halflife_mean
    1.0,    // mass_mean_sd
    5.0,    // width_mean_sd
    0.5,    // baseline_sd
    10.0,   // halflife_sd
    2.0,    // mass_sd
    8.0,    // width_sd
    0.05    // errorsq
  );

  SECTION( "All parameters initialized correctly" ) {
    REQUIRE( pop_est.mass_mean      == 3.5 );
    REQUIRE( pop_est.width_mean     == 30.0 );
    REQUIRE( pop_est.baseline_mean  == 2.5 );
    REQUIRE( pop_est.halflife_mean  == 45.0 );
    REQUIRE( pop_est.mass_mean_sd   == 1.0 );
    REQUIRE( pop_est.width_mean_sd  == 5.0 );
    REQUIRE( pop_est.baseline_sd    == 0.5 );
    REQUIRE( pop_est.halflife_sd    == 10.0 );
    REQUIRE( pop_est.mass_sd        == 2.0 );
    REQUIRE( pop_est.width_sd       == 8.0 );
    REQUIRE( pop_est.errorsq        == 0.05 );
  }

  SECTION( "Utility functions work correctly" ) {
    REQUIRE( pop_est.get_decay() == Approx(log(2) / 45.0) );
    REQUIRE( pop_est.get_logerrorsq() == Approx(log(0.05)) );
    REQUIRE( pop_est.get_mass_mean_variance() == Approx(1.0) );
    REQUIRE( pop_est.get_width_mean_variance() == Approx(25.0) );
    REQUIRE( pop_est.get_baseline_variance() == Approx(0.25) );
    REQUIRE( pop_est.get_halflife_variance() == Approx(100.0) );
    REQUIRE( pop_est.get_mass_variance() == Approx(4.0) );
    REQUIRE( pop_est.get_width_variance() == Approx(64.0) );
  }
}


//
// PopulationPriors tests
//

TEST_CASE( "PopulationPriors constructor works", "[population]" ) {

  PopulationPriors pop_priors(
    3.5, 100.0,   // mass_mean prior (mean, var)
    30.0, 100.0,  // width_mean prior
    2.5, 100.0,   // baseline_mean prior
    45.0, 100.0,  // halflife_mean prior
    5.0,          // mass_mean_sd_max
    10.0,         // width_mean_sd_max
    2.0,          // baseline_sd_max
    20.0,         // halflife_sd_max
    5.0,          // mass_sd_max
    15.0,         // width_sd_max
    1000.0,       // error_alpha
    1000.0,       // error_beta
    12.0,         // pulse_count_prior
    0.0,          // strauss_repulsion
    40.0          // strauss_repulsion_range
  );

  SECTION( "All prior parameters initialized correctly" ) {
    REQUIRE( pop_priors.mass_mean_mean         == 3.5 );
    REQUIRE( pop_priors.mass_mean_var          == 100.0 );
    REQUIRE( pop_priors.width_mean_mean        == 30.0 );
    REQUIRE( pop_priors.width_mean_var         == 100.0 );
    REQUIRE( pop_priors.baseline_mean_mean     == 2.5 );
    REQUIRE( pop_priors.baseline_mean_var      == 100.0 );
    REQUIRE( pop_priors.halflife_mean_mean     == 45.0 );
    REQUIRE( pop_priors.halflife_mean_var      == 100.0 );
    REQUIRE( pop_priors.mass_mean_sd_max       == 5.0 );
    REQUIRE( pop_priors.width_mean_sd_max      == 10.0 );
    REQUIRE( pop_priors.baseline_sd_max        == 2.0 );
    REQUIRE( pop_priors.halflife_sd_max        == 20.0 );
    REQUIRE( pop_priors.mass_sd_max            == 5.0 );
    REQUIRE( pop_priors.width_sd_max           == 15.0 );
    REQUIRE( pop_priors.error_alpha            == 1000.0 );
    REQUIRE( pop_priors.error_beta             == 0.001 );  // Note: inverted
    REQUIRE( pop_priors.pulse_count_prior      == 12.0 );
    REQUIRE( pop_priors.strauss_repulsion      == 0.0 );
    REQUIRE( pop_priors.strauss_repulsion_range == 40.0 );
  }
}


//
// Population struct tests
//

TEST_CASE( "Population struct creates and manages subjects", "[population]" ) {

  // Create some simple patient data
  NumericVector time = NumericVector::create(0, 10, 20, 30, 40);
  NumericVector conc = NumericVector::create(2.5, 3.0, 2.8, 2.6, 2.4);

  // Create patient estimates
  PatientEstimates est1(2.5, 45, 0.05, 3.5, 30, 2.0, 8.0);
  PatientEstimates est2(2.3, 50, 0.05, 3.2, 28, 2.0, 8.0);

  // Create patient priors (single-subject style for now)
  PatientPriors priors(2.5, 100, 45, 100, 3.5, 100, 30, 100,
                      2, 8, 1000, 1000, 12, 0, 40);

  // Create patient data
  PatientData data1(time, conc);
  PatientData data2(time, conc);

  // Create patients (order: data, priors, estimates)
  Patient patient1(data1, priors, est1);
  Patient patient2(data2, priors, est2);

  // Create population estimates and priors
  PopulationEstimates pop_est(3.5, 30.0, 2.5, 45.0, 1.0, 5.0, 0.5, 10.0,
                              2.0, 8.0, 0.05);
  PopulationPriors pop_priors(3.5, 100, 30, 100, 2.5, 100, 45, 100,
                              5, 10, 2, 20, 5, 15, 1000, 1000, 12, 0, 40);

  // Create population
  std::vector<Patient> subjects = {patient1, patient2};
  Population pop(subjects, pop_est, pop_priors);

  SECTION( "Population initializes correctly" ) {
    REQUIRE( pop.num_subjects == 2 );
    REQUIRE( pop.subjects.size() == 2 );
  }

  SECTION( "Population utility methods work" ) {
    REQUIRE( pop.get_total_observations() == 10 );  // 5 obs * 2 subjects
    REQUIRE( pop.get_total_pulses() == 2 );  // Each patient gets 1 initial pulse

    arma::vec baselines = pop.get_all_subject_baselines();
    REQUIRE( baselines.n_elem == 2 );
    REQUIRE( baselines(0) == 2.5 );
    REQUIRE( baselines(1) == 2.3 );

    arma::vec halflives = pop.get_all_subject_halflives();
    REQUIRE( halflives.n_elem == 2 );
    REQUIRE( halflives(0) == 45.0 );
    REQUIRE( halflives(1) == 50.0 );

    arma::vec mass_means = pop.get_all_subject_mass_means();
    REQUIRE( mass_means.n_elem == 2 );
    REQUIRE( mass_means(0) == 3.5 );
    REQUIRE( mass_means(1) == 3.2 );

    arma::vec width_means = pop.get_all_subject_width_means();
    REQUIRE( width_means.n_elem == 2 );
    REQUIRE( width_means(0) == 30.0 );
    REQUIRE( width_means(1) == 28.0 );
  }
}


//
// Sampler instantiation tests
//

TEST_CASE( "Population samplers can be instantiated", "[population]" ) {

  SECTION( "Gibbs samplers instantiate correctly" ) {
    Pop_DrawMeans means_sampler;
    Pop_DrawBaselineHalflifeMeans bh_means_sampler;
    Pop_DrawError error_sampler;

    // If we get here, samplers were created successfully
    REQUIRE( true );
  }

  SECTION( "MH samplers for subject-to-subject SDs instantiate correctly" ) {
    // These require proposal variance parameters
    Pop_DrawMeanSDs mass_mean_sd_sampler(0.1, 100, 10000, 0.4, false, false, 1000);
    Pop_DrawMeanSDs width_mean_sd_sampler(0.1, 100, 10000, 0.4, true, false, 1000);

    Pop_DrawBaselineHalflifeSD baseline_sd_sampler(0.1, 100, 10000, 0.4, false, false, 1000);
    Pop_DrawBaselineHalflifeSD halflife_sd_sampler(0.1, 100, 10000, 0.4, true, false, 1000);

    REQUIRE( true );
  }

  SECTION( "MH samplers for pulse-to-pulse SDs instantiate correctly" ) {
    Pop_DrawPulseSDs mass_sd_sampler(0.1, 100, 10000, 0.4, false, false, 1000);
    Pop_DrawPulseSDs width_sd_sampler(0.1, 100, 10000, 0.4, true, false, 1000);

    REQUIRE( true );
  }
}


//
// Pop_DrawPulseSDs::posterior_function -- log-normal parameterization
//
// The pulse-to-pulse SD sigma is a pooled population-level parameter drawn from
// all pulses across all subjects. Under the papers' (Horton) log-normal model,
// log(theta_k) ~ N(subject_log_mean, sigma^2/kappa_k): the per-pulse residuals
// entering sigma's full conditional are formed on the LOG scale,
// (log(theta_k) - subject_log_mean), weighted by the per-pulse t-scale kappa_k,
// where subject_log_mean (mass_mean/width_mean) is exactly the quantity the
// per-subject SS_DrawFixedEffects targets under log-normal -- so the three
// population full conditionals (means, mean-SD, pulse-SD) are mutually
// consistent. There is no truncation at 0, so the truncated-normal normalizing
// constants (old_norm_const/new_norm_const, the pnorm terms) are dropped. The
// Gaussian-density n*log(sigma) term and the precision term are unchanged; the
// Uniform SD prior ratio cancels in both branches.
//
TEST_CASE( "Pop_DrawPulseSDs posterior_function: log-normal pulse-SD full conditional",
           "[population][lognormal]" ) {

  NumericVector time = NumericVector::create(0, 10, 20, 30, 40);
  NumericVector conc = NumericVector::create(2.5, 3.0, 2.8, 2.6, 2.4);

  // Values chosen so the log-scale and natural-scale residuals differ sharply.
  const double mu    = 1.1;   // subject mean of LOG mass
  const double sigma = 0.8;   // current pulse-to-pulse SD (population mass_sd)
  const double prop  = 1.3;   // proposed pulse-to-pulse SD
  const double theta = 5.0;   // pulse mass (natural scale, > 0)
  const double kappa = 1.4;   // per-pulse t-scale

  // Subject with mass_mean = mu; PopulationEstimates carries mass_sd = sigma.
  PatientEstimates est(2.5, 45, 0.05, mu, 3.0, sigma, 0.9);
  PatientPriors    priors(2.5, 100, 45, 100, 3.5, 100, 30, 100,
                          2, 8, 1000, 1000, 12, 0, 40);
  PatientData      data(time, conc);
  Patient          patient(data, priors, est);

  PopulationEstimates pop_est(mu, 3.0, 2.5, 45.0, 1.0, 5.0, 0.5, 10.0,
                              sigma, 0.9, 0.05);
  PopulationPriors    pop_priors(mu, 100, 3.0, 100, 2.5, 100, 45, 100,
                                 5, 10, 2, 20, 5, 15, 1000, 1000, 12, 0, 40);

  std::vector<Patient> subjects = { patient };
  Population pop(subjects, pop_est, pop_priors);
  pop.priors.lognormal_pulses = true;   // papers' parameterization

  REQUIRE( pop.get_total_pulses() == 1 );  // single pulse -> single-term sum

  // Configure the single pulse (subject mean already set to mu above).
  PulseEstimates & pulse = pop.subjects[0].pulses.front();
  pulse.mass           = theta;
  pulse.tvarscale_mass = kappa;

  // for_width = false -> sample the pulse-to-pulse SD of mass
  Pop_DrawPulseSDs mass_sd_sampler(0.1, 100, 10000, 0.4, false, false, 1000);
  double got = mass_sd_sampler.posterior_function(&pop, prop, &pop);

  // Closed form (log-normal): log residual, NO truncation normalizing constant.
  const double logtheta    = std::log(theta);
  const double ssq         = kappa * (logtheta - mu) * (logtheta - mu);
  const double first_part  = 1 * (std::log(sigma) - std::log(prop));
  const double second_part = 0.5 * ssq *
                             (1.0 / (sigma * sigma) - 1.0 / (prop * prop));
  const double expected    = first_part + second_part;

  SECTION( "matches the analytic log-normal pulse-SD log-ratio" ) {
    REQUIRE( got == Approx(expected) );
  }

  SECTION( "differs from the natural-scale parameterization" ) {
    // Natural-scale uses (theta - mu) residuals and keeps the pnorm truncation
    // terms, so the ratio must differ -- confirming the lognormal branch is taken.
    pop.priors.lognormal_pulses = false;
    double got_natural = mass_sd_sampler.posterior_function(&pop, prop, &pop);
    REQUIRE( got != Approx(got_natural) );
  }
}
