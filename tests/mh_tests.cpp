#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif

// datastructures headers
#include <bp_datastructures/patient.h>
#include <bp_datastructures/patientdata.h>
#include <bp_datastructures/patientpriors.h>
#include <bp_datastructures/patientestimates.h>
#include <bp_datastructures/pulseestimates.h>
#include <bp_datastructures/utils.h>

// single subject model headers
#include <bpmod_singlesubject/ss_draw_baselinehalflife.h>
#include <bpmod_singlesubject/ss_draw_fixedeffects.h>
#include <bpmod_singlesubject/ss_draw_sdrandomeffects.h>
#include <bpmod_singlesubject/ss_draw_locations.h>
#include <bpmod_singlesubject/ss_draw_randomeffects.h>
#include <bpmod_singlesubject/ss_draw_tvarscale.h>
#include <bpmod_singlesubject/ss_draw_error.h>

// mcmc routine headers
#include <bp_mcmc/utils.h>
#include <bp_mcmc/proposalvariance.h>

// testing header
#include <testing/catch.h>




//----------------------------------------------------------------------
// mh_tests.cpp
//     Test MH and child classes
//
// Note: Mor tests needed.
//----------------------------------------------------------------------



//
// Test first implementation of ModifiedMetropolisHastings in
// SS_DrawFixedEffects child class
//
TEST_CASE( "first mmh test -- SS_DrawFixedEffects", "[mmh-implementations]" ) {

  //
  // Testing setup
  //
  // Create patient object -- using the default test dataset/specs
  DataStructuresUtils utils;
  PulseUtils pu;
  pu.set_seed(9999999);
  Patient pat = utils.create_new_test_patient_obj();
  Patient * patient = &pat;
  patient = utils.add_default_pulses(patient);

  // Create sampler object  (small mass pv for testing purposes)
  SS_DrawFixedEffects draw_fixed_effects_mass(0.5, 500, 25000, 0.35, false);
  SS_DrawFixedEffects draw_fixed_effects_width(30, 500, 25000, 0.35, true);


  //
  // Now, time for tests
  //

  SECTION( "Check sub-functions" ) {

    REQUIRE( draw_fixed_effects_mass.pv.getpsd() == sqrt(0.5)  );
    REQUIRE( draw_fixed_effects_mass.pv.getpv() == Approx(0.5) );
    REQUIRE( draw_fixed_effects_width.pv.getpsd() == sqrt(30)  );
    REQUIRE( draw_fixed_effects_width.pv.getpv() == Approx(30) );

  }

  SECTION( "Check tracking iterations and adjusting pv/psd" ) {

    double initial_pv, adjusted_pv, final_pv;
    int iter = 0;
    initial_pv = draw_fixed_effects_mass.pv.getpv();

    while (iter < 501) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean, iter);
      iter++;
    }

    adjusted_pv = draw_fixed_effects_mass.pv.getpv();
    REQUIRE( adjusted_pv == (initial_pv * 1.1) );

    while (iter < 25000) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean, iter);
      iter++;
    }

    adjusted_pv = draw_fixed_effects_mass.pv.getpv();

    // Test before and after the final change
    adjusted_pv = draw_fixed_effects_mass.pv.getpv();
    REQUIRE( draw_fixed_effects_mass.pv.getpv() == adjusted_pv );
    draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean, iter);
    iter++;
    REQUIRE( draw_fixed_effects_mass.pv.getpv() == adjusted_pv );
    REQUIRE( iter == 25001 );

    final_pv = draw_fixed_effects_mass.pv.getpv();
    while (iter < 50000) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean, iter);
      iter++;
    }
    REQUIRE( draw_fixed_effects_mass.pv.getpv() == final_pv );
    REQUIRE( iter == 50000 );


  }

}


TEST_CASE( "second mmh test -- SS_DrawLocationsStrauss", "[mmh-implementations]" ) {

  //
  // Testing setup
  //
  // Create patient object -- using the default test dataset/specs
  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  Patient * patient = &pat;
  patient = utils.add_default_pulses(patient);

  // Create sampler object 
  SS_DrawLocationsStrauss draw_pulse_locations_strauss(10, 500, 25000, 0.35);

  //
  // Now, time for tests
  //

  SECTION( "Check sub-functions" ) {

    REQUIRE( draw_pulse_locations_strauss.pv.getpsd() == sqrt(10)    );
    REQUIRE( draw_pulse_locations_strauss.pv.getpv() == Approx(10.0) );

  }

  SECTION( "Check tracking iterations and adjusting pv/psd" ) {

    int iter = 0;
    double initial_pv, adjusted_pv, final_pv, adjusted_psd;
    initial_pv = draw_pulse_locations_strauss.pv.getpv();

    draw_pulse_locations_strauss.sample_pulses(patient, iter);
    ++iter;

    while (iter < 501) {
      draw_pulse_locations_strauss.sample_pulses(patient, iter);
      ++iter;
    }
    adjusted_pv = draw_pulse_locations_strauss.pv.getpv();
    adjusted_psd = draw_pulse_locations_strauss.pv.getpsd();
    REQUIRE( adjusted_pv == (initial_pv * 1.1) );
    REQUIRE( adjusted_psd == Approx(sqrt(initial_pv * 1.1)) );

    while (iter < 25000) {
      draw_pulse_locations_strauss.sample_pulses(patient, iter);
      ++iter;
    }

    // Test before and after the final change
    adjusted_pv = draw_pulse_locations_strauss.pv.getpv();
    REQUIRE( draw_pulse_locations_strauss.pv.getpv() == adjusted_pv );
    REQUIRE( draw_pulse_locations_strauss.pv.getpsd() == Approx(sqrt(adjusted_pv)) );
    draw_pulse_locations_strauss.sample_pulses(patient, iter);
    ++iter;
    REQUIRE( draw_pulse_locations_strauss.pv.getpv() != adjusted_pv );

    // Test final psd change
    final_pv = draw_pulse_locations_strauss.pv.getpv();
    while (iter < 50000) {
      draw_pulse_locations_strauss.sample_pulses(patient, iter);
      ++iter;
    }

    REQUIRE( draw_pulse_locations_strauss.pv.getpv() == final_pv );
    REQUIRE( draw_pulse_locations_strauss.pv.getpsd() == Approx(sqrt(final_pv)) );

  }


}


TEST_CASE( "Temporary/partial test of all mmh objects", "[mmh-implementations]" ) {

  //
  // Testing setup
  //
  // Create patient object -- using the default test dataset/specs
  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  Patient * patient = &pat;
  patient = utils.add_default_pulses(patient);

  // Create sampler object 
  arma::vec bhl_pv = { 0.5, 45 };
  SS_DrawFixedEffects draw_fixed_effects(1.1, 500, 25000, 0.35, false);
  SS_DrawSDRandomEffects draw_sd_pulse_masses(2, 500, 25000, 0.35, false);
  SS_DrawBaselineHalflife draw_baselinehalflife(bhl_pv, 500, 25000, 0.25);
  SS_DrawLocationsStrauss draw_pulse_locations_strauss(10, 500, 25000, 0.35);
  SS_DrawRandomEffects draw_pulse_masses(1.1, 500, 25000, 0.35, false);
  SS_DrawTVarScale draw_pulse_tvarscale(1.01, 500, 25000, 0.35, false);

  SS_DrawError draw_error;


  //
  // Now, time for testing
  //

  SECTION( "Check sub-functions" ) {

    REQUIRE( draw_fixed_effects.pv.getpsd() == sqrt(1.1)  );
    REQUIRE( draw_fixed_effects.pv.getpv() == Approx(1.1) );

    REQUIRE( draw_sd_pulse_masses.pv.getpsd() == sqrt(2)  );
    REQUIRE( draw_sd_pulse_masses.pv.getpv() == Approx(2) );

    double x, y, xy;
    x = 0.5; y = 45; xy = -0.9 * sqrt(x * y);
    arma::mat checkpv = { { x, xy }, { xy, y } };
    arma::mat checkchol = arma::chol(checkpv);
    REQUIRE( arma::approx_equal(draw_baselinehalflife.pv.getpsd(), checkchol,
                                "absdiff", 0.0000001) );
    REQUIRE( arma::approx_equal(draw_baselinehalflife.pv.getpv(), checkpv,
                                "absdiff", 0.0000001) );

    REQUIRE( draw_pulse_locations_strauss.pv.getpsd() == sqrt(10) );
    REQUIRE( draw_pulse_locations_strauss.pv.getpv() == Approx(10.0) );

    REQUIRE( draw_pulse_masses.pv.getpsd() == sqrt(1.1) );
    REQUIRE( draw_pulse_masses.pv.getpv() == Approx(1.1) );

    REQUIRE( draw_pulse_tvarscale.pv.getpsd() == sqrt(1.01) );
    REQUIRE( draw_pulse_tvarscale.pv.getpv() == Approx(1.01) );

  }

  SECTION( "Temporary test section - Run sampler for other MMH objects" ) {

    int iter = 0;
    double x, y, xy;
    x = 0.5; y = 45; xy = -0.9 * sqrt(x * y);
    arma::mat checkpv = { { x, xy }, { xy, y } };
    double pvfe     = draw_fixed_effects.pv.getpv();
    double pvsd     = draw_sd_pulse_masses.pv.getpv();
    double pvloc    = draw_pulse_locations_strauss.pv.getpv();
    double pvpmass  = draw_pulse_masses.pv.getpv();
    double pvpscale = draw_pulse_tvarscale.pv.getpv();

    while (iter < 1500) {

      draw_fixed_effects.sample(patient, &patient->estimates.mass_mean, iter);
      draw_sd_pulse_masses.sample(patient, &patient->estimates.mass_sd, patient,
                                  iter);
      draw_baselinehalflife.sample(patient,
                                   &patient->estimates.baseline_halflife, iter);
      draw_pulse_locations_strauss.sample_pulses(patient, iter);
      draw_pulse_masses.sample_pulses(patient, iter);
      draw_pulse_tvarscale.sample_pulses(patient, iter);
      draw_error.sample(patient);

      ++iter;

    }

    REQUIRE( draw_fixed_effects.pv.getpv()           != pvfe );
    REQUIRE( draw_sd_pulse_masses.pv.getpv()         != pvsd );
    REQUIRE( draw_pulse_locations_strauss.pv.getpv() != pvloc );
    REQUIRE( draw_pulse_masses.pv.getpv()            != pvpmass );
    REQUIRE( draw_pulse_tvarscale.pv.getpv()         != pvpscale );
    REQUIRE( !arma::approx_equal(draw_baselinehalflife.pv.getpv(), checkpv,
                                 "absdiff", 0.0000001) );

  }

}


