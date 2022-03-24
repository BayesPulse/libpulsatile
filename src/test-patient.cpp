#include <testthat.h>
#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bpmod_singlesubject/bpmod_singlesubject.h>
// #include <testing/catch.h>

context("Patient - Single Patient Constructor") {
  
  NumericVector time(144);
  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  
  test_that( "Estimates can be accessed" ) {
    expect_true(pat.estimates.baseline_halflife(0) == 2.6);
    expect_true(pat.estimates.mass_mean == 3.5);
    expect_true(pat.estimates.mass_sd == 1/10);
  }
  
  test_that( "Estimates can be updated" ) {
    pat.estimates.baseline_halflife(0) = 10;
    expect_true(pat.estimates.baseline_halflife(0) == 10);
    pat.estimates.mass_mean = 5.0;
    expect_true(pat.estimates.mass_mean == 5.0);
    pat.estimates.mass_prec = 1/50^2;
    expect_true(pat.estimates.mass_prec == 1/50^2);
    expect_true(pat.estimates.get_decay() == (log(2)/45));
  }
  
  test_that( "Priors can be accessed" ) {
    expect_true(pat.priors.baseline_mean == 2.6);
    expect_true(pat.priors.mass_prec_param == 1/5^2);
    expect_true(pat.priors.error_alpha == 0.0001);
    expect_true(pat.priors.num_orderstat == 3);
    expect_true(pat.priors.strauss_repulsion == 0);
  }
  
  test_that( "Priors can be updated" ) {
    pat.priors.baseline_mean     = 2.25;
    pat.priors.mass_prec_param      = 0.001;
    pat.priors.mass_prec_param_rate = 0.001;
    pat.priors.error_alpha       = 700;
    pat.priors.num_orderstat     = 4;
    pat.priors.strauss_repulsion = 0.75;
    expect_true(pat.priors.baseline_mean     == 2.25);
    expect_true(pat.priors.mass_prec_param   == 0.001);
    expect_true(pat.priors.mass_prec_param_rate == 0.001);
    expect_true(pat.priors.error_alpha       == 700);
    expect_true(pat.priors.num_orderstat     == 4);
    expect_true(pat.priors.strauss_repulsion == 0.75);
  }
  
  test_that( "Data can be accessed" ) {
    expect_true(pat.data.time(1)              == 20);
    expect_true(pat.data.time(143)            == 1440);
    expect_true(pat.data.concentration(1)     == Approx(log(3.619304)));
    expect_true(pat.data.concentration(143)   == Approx(log(3.079000)));
    expect_true(pat.data.time.size()          == 144);
    expect_true(pat.data.concentration.size() == 144);
    expect_true(pat.data.response_concentration.size() == 0);
    expect_true(pat.data.avg_period_of_obs    == 10);
    expect_true(pat.data.duration_of_obs      == 1430);
    expect_true(pat.data.number_of_obs        == 144);
    expect_true(pat.data.fitstart             == -40);
    expect_true(pat.data.fitend               == 1460);
  }
  
  // TODO: This test_that is causing trouble (initial pulse values), not giving
  // sensible results...
  test_that( "Has initial pulse" ) {
    std::list<PulseEstimates>::const_iterator this_piter = pat.pulses.begin();
    expect_true(pat.get_pulsecount()        == 1);
    expect_true(this_piter->time            == pat.data.fitstart);
    expect_true(this_piter->mass            == 1);
    expect_true(this_piter->width           == 1);
    expect_true(this_piter->tvarscale_mass  == 1);
    expect_true(this_piter->tvarscale_width == 1);
  }
  
  // Add pulses
  Patient * patient = &pat;
  patient = utils.add_default_pulses(patient);
  arma::vec location { 26.54152, 174.63993, 298.62117, 360.55329, 494.61155,
                       689.09242, 763.89017, 839.80027, 925.80251, 975.47320, 
                       1199.00866, 1322.82471 };
  
  test_that( "Can add pulses and iterate with iterators" ) {
    ++pat.piter;
    expect_true(pat.get_pulsecount() == 12);
    expect_true(pat.piter->time == Approx(location(0)));
    ++pat.piter;
    expect_true(pat.piter->time == Approx(location(1)));
    ++pat.piter;
    expect_true(pat.piter->time == Approx(location(2)));
    ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter;
    ++pat.piter; ++pat.piter; ++pat.piter;
    expect_true(pat.piter->time == Approx(location(10)));
  }
  
  test_that( "Can remove a pulse" ) {
    ++pat.piter; ++pat.piter;
    expect_true(pat.piter->time == Approx(location(1)));
    pat.piter = pat.pulses.erase(pat.piter); // this is how you delete and keep iter correct
    expect_true(pat.get_pulsecount() == 11);
    expect_true(pat.piter->time == Approx(location(2)));
  }
  
  // Get mean concentration of all pulses
  arma::vec mconc = pat.mean_concentration(false);
  
  test_that( "Can get mean_concentration" ) {
    
    expect_true(mconc.n_elem == 144);
    expect_true(mconc.n_elem == 144);
    expect_true(mconc.min() > 0);
    expect_true(mconc.max() < 10);
    
  }
  
  test_that( "Can get mean_concentration, after removing one pulse" ) {
    
    ++pat.piter;
    expect_true( pat.piter->time == Approx( location(0) ));
    arma::vec mconc_excl1 = pat.mean_concentration(false, pat.piter);
    expect_true( mconc_excl1.n_elem == 144 );
    
    expect_true( !arma::approx_equal(mconc_excl1, mconc, "absdiff", 0.0000001) );
    expect_true( arma::all(mconc_excl1 <= mconc) );
    expect_true( arma::all(mconc_excl1 >= 0) );
    expect_true( arma::all(mconc_excl1 < 10) );
    
  }
  
  test_that( "Can get mean_concentration, after removing one pulse" ) {
    
    // move to another pulse for excluding
    ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter;
    expect_true( pat.piter->time == Approx( location(3) ));
    
    // calc mean_conc excluding one pulse
    arma::vec mconc_excl3 = pat.mean_concentration(false, pat.piter);
    expect_true( mconc_excl3.n_elem == 144 );
    
    // expect_true mconc excluding 3 (#4) to be <= and not all == to full mconc
    expect_true( !arma::approx_equal(mconc_excl3, mconc, "absdiff", 0.0000001) );
    expect_true( arma::all(mconc_excl3 <= mconc) );
    expect_true( arma::all(mconc_excl3 >= 0) );
    expect_true( arma::all(mconc_excl3 < 10) );
    
  }
  
  test_that( "Can get likelihood" ) {
    
    //expect_true(pat.likelihood(false) == Approx(71.8522693142)); changed w/ 
    expect_true(pat.likelihood(false) == Approx(162.4336751778));
    
  }
  
}

context("Patient - Population Constructor") {
  NumericVector time(144);
  NumericVector conc = rnorm(144, 3, 0.1);
  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;
  
  PatientData data(time, conc);
  PatientEstimates estimates(2.6, 45, 0.05, 3.5, 30, 10, 10, false); // population constructor
  //PatientData * data = &pd;
  //PatientEstimates * estimates = &pep;
  
  //Patient pat(data, estimates);
  Patient pat(data, estimates);
  
  test_that( "Estimates can be accessed" ) {
    expect_true(pat.estimates.baseline_halflife(0) == 2.6);
    expect_true(pat.estimates.mass_mean == 3.5);
  }
  
  test_that( "Estimates can be updated" ) {
    pat.estimates.baseline_halflife(0) = 10;
    expect_true(pat.estimates.baseline_halflife(0) == 10);
    pat.estimates.mass_mean = 5.0;
    expect_true(pat.estimates.mass_mean == 5.0);
  }
  
  test_that( "Data can be accessed" ) {
    expect_true(pat.data.time(1) == 20);
    expect_true(pat.data.time(143) == 1440);
    expect_true(pat.data.concentration(1) < 20);
    expect_true(pat.data.concentration(1) > 0);
    expect_true(pat.data.concentration(143) < 20);
    expect_true(pat.data.concentration(143) > 0);
    expect_true(pat.data.time.size() == 144);
    expect_true(pat.data.concentration.size() == 144);
    expect_true(pat.data.response_concentration.size() == 0);
  }
}
