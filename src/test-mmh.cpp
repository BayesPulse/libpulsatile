#include <testthat.h>
#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif

// single subject model headers
#include <bpmod_singlesubject/bpmod_singlesubject.h>
#include <testing/catch.h>

//----------------------------------------------------------------------
// mh_tests.cpp
//     Test MH and child classes
//
// Note: More tests needed.
//----------------------------------------------------------------------

//
// Test first implementation of ModifiedMetropolisHastings in
// SS_DrawFixedEffects child class
//
context( "MMH - SS_DrawFixedEffects" ) {
  
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
  SS_DrawFixedEffects draw_fixed_effects_mass(0.5, 500, 25000, 0.35, false, false, 5000);
  SS_DrawFixedEffects draw_fixed_effects_width(30, 500, 25000, 0.35, true, false, 5000);
  
  
  //
  // Now, time for tests
  //
  
  test_that( "Check sub-functions" ) {
    
    expect_true( draw_fixed_effects_mass.pv.getpsd() == sqrt(0.5)  );
    expect_true( draw_fixed_effects_mass.pv.getpv() == Approx(0.5) );
    expect_true( draw_fixed_effects_width.pv.getpsd() == sqrt(30)  );
    expect_true( draw_fixed_effects_width.pv.getpv() == Approx(30) );
    
  }
  
  test_that( "Check tracking iterations and adjusting pv/psd" ) {
    
    double initial_pv, adjusted_pv, final_pv;
    int iter = 0;
    initial_pv = draw_fixed_effects_mass.pv.getpv();
    
    while (iter < 501) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean, iter);
      iter++;
    }
    
    adjusted_pv = draw_fixed_effects_mass.pv.getpv();
    expect_true( adjusted_pv == (initial_pv * 1.1) );
    
    while (iter < 25000) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean, iter);
      iter++;
    }
    
    adjusted_pv = draw_fixed_effects_mass.pv.getpv();
    
    // Test before and after the final change
    adjusted_pv = draw_fixed_effects_mass.pv.getpv();
    expect_true( draw_fixed_effects_mass.pv.getpv() == adjusted_pv );
    draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean, iter);
    iter++;
    expect_true( draw_fixed_effects_mass.pv.getpv() == adjusted_pv );
    expect_true( iter == 25001 );
    
    final_pv = draw_fixed_effects_mass.pv.getpv();
    while (iter < 50000) {
      draw_fixed_effects_mass.sample(patient, &patient->estimates.mass_mean, iter);
      iter++;
    }
    expect_true( draw_fixed_effects_mass.pv.getpv() == final_pv );
    expect_true( iter == 50000 );
    
    
  }
  
}


context( "MMH - SS_DrawLocationsStrauss" ) {
  
  //
  // Testing setup
  //
  // Create patient object -- using the default test dataset/specs
  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  Patient * patient = &pat;
  patient = utils.add_default_pulses(patient);
  
  // Create sampler object 
  SS_DrawLocationsStrauss draw_pulse_locations_strauss(10, 500, 25000, 0.35, false, 5000);
  
  //
  // Now, time for tests
  //
  
  test_that( "Check sub-functions" ) {
    
    expect_true( draw_pulse_locations_strauss.pv.getpsd() == sqrt(10)    );
    expect_true( draw_pulse_locations_strauss.pv.getpv() == Approx(10.0) );
    
  }
  
  test_that( "Check tracking iterations and adjusting pv/psd" ) {
    
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
    expect_true( adjusted_pv == (initial_pv * 1.1) );
    expect_true( adjusted_psd == Approx(sqrt(initial_pv * 1.1)) );
    
    while (iter < 25000) {
      draw_pulse_locations_strauss.sample_pulses(patient, iter);
      ++iter;
    }
    
    // Test before and after the final change
    adjusted_pv = draw_pulse_locations_strauss.pv.getpv();
    expect_true( draw_pulse_locations_strauss.pv.getpv() == adjusted_pv );
    expect_true( draw_pulse_locations_strauss.pv.getpsd() == Approx(sqrt(adjusted_pv)) );
    draw_pulse_locations_strauss.sample_pulses(patient, iter);
    ++iter;
    // Note: failing this test -- believe it's because the decision to adjust pv
    // is based on success rate of the mh algo and the new Cauchy prior on sds
    // chages the decision basis/likelihood/estimates at this point.  look into
    // more.
    expect_true( draw_pulse_locations_strauss.pv.getpv() != adjusted_pv );
    
    // Test final psd change
    final_pv = draw_pulse_locations_strauss.pv.getpv();
    while (iter < 50000) {
      draw_pulse_locations_strauss.sample_pulses(patient, iter);
      ++iter;
    }
    
    expect_true( draw_pulse_locations_strauss.pv.getpv() == final_pv );
    expect_true( draw_pulse_locations_strauss.pv.getpsd() == Approx(sqrt(final_pv)) );
    
  }
  
  
}

context( "MMH - Temp/Partial tests" ) {
  
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
  SS_DrawFixedEffects draw_fixed_effects(1.1, 500, 25000, 0.35, false, false, 5000);
  SS_DrawSDRandomEffects draw_sd_pulse_masses(2, 500, 25000, 0.35, false, false, 5000);
  SS_DrawBaselineHalflife draw_baselinehalflife(bhl_pv, 500, 25000, 0.25, false, 5000);
  SS_DrawLocationsStrauss draw_pulse_locations_strauss(10, 500, 25000, 0.35, false, 5000);
  SS_DrawRandomEffects draw_pulse_masses(1.1, 500, 25000, 0.35, false, false, 5000);
  SS_DrawTVarScale draw_pulse_tvarscale(1.01, 500, 25000, 0.35, false, false, 5000);
  
  SS_DrawError draw_error;
  
  
  //
  // Now, time for testing
  //
  
  test_that( "Check sub-functions" ) {
    
    expect_true( draw_fixed_effects.pv.getpsd() == sqrt(1.1)  );
    expect_true( draw_fixed_effects.pv.getpv() == Approx(1.1) );
    
    expect_true( draw_sd_pulse_masses.pv.getpsd() == sqrt(2)  );
    expect_true( draw_sd_pulse_masses.pv.getpv() == Approx(2) );
    
    double x, y, xy;
    x = 0.5; y = 45; xy = -0.9 * sqrt(x * y);
    arma::mat checkpv = { { x, xy }, { xy, y } };
    arma::mat checkchol = arma::chol(checkpv);
    expect_true( arma::approx_equal(draw_baselinehalflife.pv.getpsd(), checkchol,
                                "absdiff", 0.0000001) );
    expect_true( arma::approx_equal(draw_baselinehalflife.pv.getpv(), checkpv,
                                "absdiff", 0.0000001) );
    
    expect_true( draw_pulse_locations_strauss.pv.getpsd() == sqrt(10) );
    expect_true( draw_pulse_locations_strauss.pv.getpv() == Approx(10.0) );
    
    expect_true( draw_pulse_masses.pv.getpsd() == sqrt(1.1) );
    expect_true( draw_pulse_masses.pv.getpv() == Approx(1.1) );
    
    expect_true( draw_pulse_tvarscale.pv.getpsd() == sqrt(1.01) );
    expect_true( draw_pulse_tvarscale.pv.getpv() == Approx(1.01) );
    
  }
  
  test_that( "Temporary test section - Run sampler for other MMH objects" ) {
    
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
    
    expect_true( draw_fixed_effects.pv.getpv()           != pvfe );
    expect_true( draw_sd_pulse_masses.pv.getpv()         != pvsd );
    expect_true( draw_pulse_locations_strauss.pv.getpv() != pvloc );
    expect_true( draw_pulse_masses.pv.getpv()            != pvpmass );
    expect_true( draw_pulse_tvarscale.pv.getpv()         != pvpscale );
    expect_true( !arma::approx_equal(draw_baselinehalflife.pv.getpv(), checkpv,
                                 "absdiff", 0.0000001) );
    
  }
  
}