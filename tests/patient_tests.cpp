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
#include <catch.h>


//
// patient_tests.cpp
//     Patient class tests
//
// note: haven't added tests of response_concentration data
//



//
// Single-subject patient (OnlyPatient)
//

TEST_CASE( "Patient class constructor for single-subject works", "[patient]" ) {

  NumericVector time(144);
  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();

  SECTION( "Estimates can be accessed" ) {
    REQUIRE(pat.estimates.baseline_halflife(0) == 2.6);
    REQUIRE(pat.estimates.mass_mean == 3.5);
    REQUIRE(pat.estimates.mass_sd == 10);
  }

  SECTION( "Estimates can be updated" ) {
    pat.estimates.baseline_halflife(0) = 10;
    REQUIRE(pat.estimates.baseline_halflife(0) == 10);
    pat.estimates.mass_mean = 5.0;
    REQUIRE(pat.estimates.mass_mean == 5.0);
    pat.estimates.mass_sd = 50;
    REQUIRE(pat.estimates.mass_sd == 50);
    REQUIRE(pat.estimates.get_decay() == (log(2)/45));
  }

  SECTION( "Priors can be accessed" ) {
    REQUIRE(pat.priors.baseline_mean == 1.5);
    REQUIRE(pat.priors.mass_sd_max == 10);
    REQUIRE(pat.priors.error_alpha == 1000);
    REQUIRE(pat.priors.num_orderstat == 3);
    REQUIRE(pat.priors.strauss_repulsion == 0);
  }

  SECTION( "Priors can be updated" ) {
    pat.priors.baseline_mean     = 2.25;
    pat.priors.mass_sd_max       = 90;
    pat.priors.error_alpha       = 700;
    pat.priors.num_orderstat     = 4;
    pat.priors.strauss_repulsion = 0.75;
    REQUIRE(pat.priors.baseline_mean     == 2.25);
    REQUIRE(pat.priors.mass_sd_max       == 90);
    REQUIRE(pat.priors.error_alpha       == 700);
    REQUIRE(pat.priors.num_orderstat     == 4);
    REQUIRE(pat.priors.strauss_repulsion == 0.75);
  }

  SECTION( "Data can be accessed" ) {
    REQUIRE(pat.data.time(1)              == 20);
    REQUIRE(pat.data.time(143)            == 1440);
    REQUIRE(pat.data.concentration(1)     == log(5.156800));
    REQUIRE(pat.data.concentration(143)   == log(8.155648));
    REQUIRE(pat.data.time.size()          == 144);
    REQUIRE(pat.data.concentration.size() == 144);
    REQUIRE(pat.data.response_concentration.size() == 0);
    REQUIRE(pat.data.avg_period_of_obs    == 10);
    REQUIRE(pat.data.duration_of_obs      == 1430);
    REQUIRE(pat.data.number_of_obs        == 144);
    REQUIRE(pat.data.fitstart             == -40);
    REQUIRE(pat.data.fitend               == 1460);
  }

  // TODO: This section is causing trouble (initial pulse values), not giving
  // sensible results...
  SECTION( "Has initial pulse" ) {
    std::list<PulseEstimates>::const_iterator this_piter = pat.pulses.begin();
    REQUIRE(pat.get_pulsecount()        == 1);
    REQUIRE(this_piter->time            == pat.data.fitstart);
    REQUIRE(this_piter->mass            == 1);
    REQUIRE(this_piter->width           == 1);
    REQUIRE(this_piter->tvarscale_mass  == 1);
    REQUIRE(this_piter->tvarscale_width == 1);
  }

  // Add pulses
  Patient * patient = &pat;
  patient = utils.add_default_pulses(patient);
  arma::vec location { -37.65204, 41.04917, 125.58297, 306.32966, 461.65469,
    616.95867, 835.64466, 1000.29703, 1149.08360, 1319.19888, 1414.31830 };

  SECTION( "Can add pulses and iterate with iterators" ) {
    ++pat.piter;
    REQUIRE(pat.get_pulsecount() == 11);
    REQUIRE(pat.piter->time == location(0));
    ++pat.piter;
    REQUIRE(pat.piter->time == location(1));
    ++pat.piter;
    REQUIRE(pat.piter->time == location(2));
    ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter;
    ++pat.piter; ++pat.piter; ++pat.piter;
    REQUIRE(pat.piter->time == location(10));
  }

  SECTION( "Can remove a pulse" ) {
    ++pat.piter; ++pat.piter;
    REQUIRE(pat.piter->time == location(1));
    pat.piter = pat.pulses.erase(pat.piter); // this is how you delete and keep iter correct
    REQUIRE(pat.get_pulsecount() == 10);
    REQUIRE(pat.piter->time == location(2));
  }

  // Get mean concentration of all pulses
  arma::vec mconc = pat.mean_concentration(false);

  SECTION( "Can get mean_concentration" ) {

    REQUIRE(mconc.n_elem == 144);
    REQUIRE(mconc.n_elem == 144);
    REQUIRE(mconc.min() > 0);
    REQUIRE(mconc.max() < 10);

  }

  SECTION( "Can get mean_concentration, after removing one pulse" ) {

    ++pat.piter;
    REQUIRE( pat.piter->time == location(0) );
    arma::vec mconc_excl1 = pat.mean_concentration(false, pat.piter);
    REQUIRE( mconc_excl1.n_elem == 144 );

    REQUIRE( !arma::approx_equal(mconc_excl1, mconc, "absdiff", 0.0000001) );
    REQUIRE( arma::all(mconc_excl1 <= mconc) );
    REQUIRE( arma::all(mconc_excl1 >= 0) );
    REQUIRE( arma::all(mconc_excl1 < 10) );

  }

  SECTION( "Can get mean_concentration, after removing one pulse" ) {

    // move to another pulse for excluding
    ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter;
    REQUIRE( pat.piter->time == location(3) );

    // calc mean_conc excluding one pulse
    arma::vec mconc_excl3 = pat.mean_concentration(false, pat.piter);
    REQUIRE( mconc_excl3.n_elem == 144 );

    // Require mconc excluding 3 (#4) to be <= and not all == to full mconc
    REQUIRE( !arma::approx_equal(mconc_excl3, mconc, "absdiff", 0.0000001) );
    REQUIRE( arma::all(mconc_excl3 <= mconc) );
    REQUIRE( arma::all(mconc_excl3 >= 0) );
    REQUIRE( arma::all(mconc_excl3 < 10) );

  }

  SECTION( "Can get likelihood" ) {

    REQUIRE(pat.likelihood(false) == Approx(83.36087));

  }

}



//
// Population patient
//

TEST_CASE( "Patient class constructor for population model works", "[patient]" ) {

  NumericVector time(144);
  NumericVector conc = rnorm(144, 3, 0.1);
  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;

  PatientData data(time, conc);
  PatientEstimates estimates(2.6, 45, 0.05, 3.5, 30); // population constructor
  //PatientData * data = &pd;
  //PatientEstimates * estimates = &pep;

  //Patient pat(data, estimates);
  Patient pat(data, estimates);

  SECTION( "Estimates can be accessed" ) {
    REQUIRE(pat.estimates.baseline_halflife(0) == 2.6);
    REQUIRE(pat.estimates.mass_mean == 3.5);
  }

  SECTION( "Estimates can be updated" ) {
    pat.estimates.baseline_halflife(0) = 10;
    REQUIRE(pat.estimates.baseline_halflife(0) == 10);
    pat.estimates.mass_mean = 5.0;
    REQUIRE(pat.estimates.mass_mean == 5.0);
  }

  SECTION( "Data can be accessed" ) {
    REQUIRE(pat.data.time(1) == 20);
    REQUIRE(pat.data.time(143) == 1440);
    REQUIRE(pat.data.concentration(1) < 20);
    REQUIRE(pat.data.concentration(1) > 0);
    REQUIRE(pat.data.concentration(143) < 20);
    REQUIRE(pat.data.concentration(143) > 0);
    REQUIRE(pat.data.time.size() == 144);
    REQUIRE(pat.data.concentration.size() == 144);
    REQUIRE(pat.data.response_concentration.size() == 0);
  }

  // TODO: look at what is necessary for the pop constructor version
  //SECTION( "Can add a pulse" ) {
  //}
  //SECTION( "Can remove a pulse" ) {
  //}

}








