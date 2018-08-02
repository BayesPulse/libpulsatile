#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif

#include <bp_datastructures/patient.h>
#include <bp_datastructures/patientdata.h>
#include <bp_datastructures/patientpriors.h>
#include <bp_datastructures/patientestimates.h>
#include <bp_datastructures/pulseestimates.h>
#include <bp_datastructures/utils.h>

#include <bpmod_singlesubject/ss_draw_baselinehalflife.h>
#include <bpmod_singlesubject/ss_draw_fixedeffects.h>
#include <bpmod_singlesubject/ss_draw_sdrandomeffects.h>
#include <bpmod_singlesubject/ss_draw_locations.h>
#include <bpmod_singlesubject/ss_draw_randomeffects.h>
#include <bpmod_singlesubject/ss_draw_tvarscale.h>
#include <bpmod_singlesubject/ss_draw_error.h>

#include <bp_mcmc/utils.h>
#include <bp_mcmc/proposalvariance.h>

#include <testing/catch.h>

//
// mh_tests.cpp
//     Test MH and child classes
//
// Note: Mor tests needed.
//



//
// Test first implementation of ModifiedMetropolisHastings in
// SS_DrawFixedEffects child class
//

TEST_CASE( "first mmh test -- SS_DrawFixedEffects", "[mmh-implementations]" ) {

  // Create patient object -- using the default test dataset/specs
  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  Patient * patient = &pat;
  patient = utils.add_default_pulses(patient);

  // Create sampler object 
  SS_DrawFixedEffects draw_fixed_effects_mass(1.1, 500, 25000, 0.35, false);
  SS_DrawFixedEffects draw_fixed_effects_width(30, 500, 25000, 0.35, true);


  //
  // Now, time for testing
  //

  SECTION( "Check sub-functions" ) {

    REQUIRE(draw_fixed_effects_mass.pv.getpsd() == sqrt(1.1));
    REQUIRE(draw_fixed_effects_mass.pv.getpv() == Approx(1.1));
    REQUIRE(draw_fixed_effects_width.pv.getpsd() == sqrt(30));
    REQUIRE(draw_fixed_effects_width.pv.getpv() == Approx(30));

  }

  SECTION( "Check tracking iterations and adjusting pv/psd" ) {

    double initial_psd, adjusted_psd, final_psd;
    initial_psd = draw_fixed_effects_mass.pv.getpsd();
    for (int i = 0; i < 501; i++) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean);
    }
    adjusted_psd = draw_fixed_effects_mass.pv.getpsd();
    REQUIRE(adjusted_psd == sqrt((pow(initial_psd, 2)*1.1)));

    for (int i = 500; i < 24499; i++) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean);
    }

    // Test before and after the final change
    adjusted_psd = draw_fixed_effects_mass.pv.getpsd();
    REQUIRE(draw_fixed_effects_mass.pv.getpsd() == adjusted_psd);
    draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean);
    REQUIRE(draw_fixed_effects_mass.pv.getpsd() != adjusted_psd);

    final_psd = draw_fixed_effects_mass.pv.getpsd();
    for (int i = draw_fixed_effects_mass.pv.getiter(); i < 50000; i++) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean);
    }
    REQUIRE(draw_fixed_effects_mass.pv.getpsd() == final_psd);
    REQUIRE(draw_fixed_effects_mass.pv.getiter() == 50000);


  }

}


TEST_CASE( "second mmh test -- SS_DrawLocationsStrauss", "[mmh-implementations]" ) {

  // Create patient object -- using the default test dataset/specs
  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  Patient * patient = &pat;
  patient = utils.add_default_pulses(patient);

  // Create sampler object 
  SS_DrawLocationsStrauss draw_pulse_locations_strauss(10, 500*11, 25000*11, 0.35);

  //
  // Now, time for testing
  //

  SECTION( "Check sub-functions" ) {

    REQUIRE(draw_pulse_locations_strauss.pv.getpsd() == sqrt(10));
    REQUIRE(draw_pulse_locations_strauss.pv.getpv() == Approx(10.0));

  }

  // TODO: rewrite this to loop over pulses once, then a bunch
  SECTION( "Check tracking iterations and adjusting pv/psd" ) {

    double initial_psd, adjusted_psd, final_psd;
    initial_psd = draw_pulse_locations_strauss.pv.getpsd();

    draw_pulse_locations_strauss.sample_pulses(patient);

    for (int i = 0; i < 500; i++) {
      draw_pulse_locations_strauss.sample_pulses(patient);
    }
    adjusted_psd = draw_pulse_locations_strauss.pv.getpsd();
    REQUIRE(adjusted_psd == sqrt((pow(initial_psd, 2)*1.1)));

    for (int i = 500; i < 24499; i++) {
      draw_pulse_locations_strauss.sample_pulses(patient);
    }

    // Test before and after the final change
    adjusted_psd = draw_pulse_locations_strauss.pv.getpsd();
    REQUIRE(draw_pulse_locations_strauss.pv.getpsd() == adjusted_psd);
    draw_pulse_locations_strauss.sample_pulses(patient);
    REQUIRE(draw_pulse_locations_strauss.pv.getpsd() != adjusted_psd);

    // Test final psd change
    final_psd = draw_pulse_locations_strauss.pv.getpsd();
    for (int i = 24500; i < 50000; i++) {
      draw_pulse_locations_strauss.sample_pulses(patient);
    }
    REQUIRE(draw_pulse_locations_strauss.pv.getpsd() == final_psd);
    REQUIRE(draw_pulse_locations_strauss.pv.getiter() == 550011);

  }


}


TEST_CASE( "Temporary/partial test of all mmh objects", "[mmh-implementations]" ) {

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

  SS_DrawLocationsStrauss draw_pulse_locations_strauss(10, 500*11, 25000*11, 0.35);
  SS_DrawRandomEffects draw_pulse_masses(1.1, 500*11, 25000*11, 0.35, false);
  SS_DrawTVarScale draw_pulse_tvarscale(1.01, 500*11, 25000*11, 0.35);

  SS_DrawError draw_error;


  //
  // Now, time for testing
  //

  SECTION( "Check sub-functions" ) {

    REQUIRE( draw_fixed_effects.pv.getpsd() == sqrt(1.1) );
    REQUIRE( draw_fixed_effects.pv.getpv() == Approx(1.1) );

    REQUIRE( draw_sd_pulse_masses.pv.getpsd() == sqrt(2) );
    REQUIRE( draw_sd_pulse_masses.pv.getpv() == Approx(2) );

    double x, y, xy;
    x = 0.5; y = 45; xy = -0.9 * sqrt(x * y);
    arma::mat checkpv = { { x, xy }, { xy, y } };
    arma::mat checkchol = arma::chol(checkpv);
    REQUIRE( arma::approx_equal(draw_baselinehalflife.pv.getpsd(), checkchol, "absdiff", 0.0000001) );
    REQUIRE( arma::approx_equal(draw_baselinehalflife.pv.getpv(), checkpv, "absdiff", 0.0000001) );

    REQUIRE( draw_pulse_locations_strauss.pv.getpsd() == sqrt(10) );
    REQUIRE( draw_pulse_locations_strauss.pv.getpv() == Approx(10.0) );

    REQUIRE( draw_pulse_masses.pv.getpsd() == sqrt(1.1) );
    REQUIRE( draw_pulse_masses.pv.getpv() == Approx(1.1) );

    REQUIRE( draw_pulse_tvarscale.pv.getpsd() == sqrt(1.01) );
    REQUIRE( draw_pulse_tvarscale.pv.getpv() == Approx(1.01) );

  }

  SECTION( "Temporary test section - Run sampler for other MMH objects" ) {

    double x, y, xy;
    x = 0.5; y = 45; xy = -0.9 * sqrt(x * y);
    arma::mat checkpv = { { x, xy }, { xy, y } };
    double pvfe     = draw_fixed_effects.pv.getpv();
    double pvsd     = draw_sd_pulse_masses.pv.getpv();
    double pvloc    = draw_pulse_locations_strauss.pv.getpv();
    double pvpmass  = draw_pulse_masses.pv.getpv();
    double pvpscale = draw_pulse_tvarscale.pv.getpv();

    for (int i = 0; i < 1500; i++) {
      draw_fixed_effects.sample(patient, &patient->estimates.mass_mean);
      draw_sd_pulse_masses.sample(patient, &patient->estimates.mass_sd, patient);
      draw_baselinehalflife.sample(patient, &patient->estimates.baseline_halflife);
      draw_pulse_locations_strauss.sample_pulses(patient);
      draw_pulse_masses.sample_pulses(patient);
      draw_pulse_tvarscale.sample_pulses(patient);
      draw_error.sample(patient);

      //std::cout << "Draw " << i << " errorsq = " << patient->estimates.errorsq << std::endl;;
      //std::cout << "Draw " << i << "; Baseline = " << patient->estimates.baseline_halflife(0) << " ; Halflife = " << patient->estimates.baseline_halflife(1) << std::endl; 

    }

    REQUIRE( draw_fixed_effects.pv.getpv()           != pvfe );
    REQUIRE( draw_sd_pulse_masses.pv.getpv()         != pvsd );
    REQUIRE( draw_pulse_locations_strauss.pv.getpv() != pvloc );
    REQUIRE( draw_pulse_masses.pv.getpv()            != pvpmass );
    REQUIRE( draw_pulse_tvarscale.pv.getpv()         != pvpscale );
    REQUIRE( !arma::approx_equal(draw_baselinehalflife.pv.getpv(), checkpv, "absdiff", 0.0000001) );

  }

}


