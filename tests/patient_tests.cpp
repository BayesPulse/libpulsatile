#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/patient.h>
#include <bp_datastructures/patientdata.h>
#include <bp_datastructures/patientpriors.h>
#include <bp_datastructures/patientestimates.h>
#include <bp_datastructures/pulseestimates.h>
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
  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;
  NumericVector conc =
    { 5.237937,  5.156800,  5.773434,  6.364134,  9.872640,  8.782145, 8.633523,
      7.166886,  6.063834,  6.393283,  5.547755,  5.609554, 17.809125,
      14.864233, 14.488007, 11.012978, 10.074941,  8.579306,  7.925953,
      7.505312,  7.408293, 6.296038,  6.138372,  6.114483,  4.589039,  4.787283,
      4.125096,  4.560712, 3.838159,  3.611227, 10.149188,  9.114352,  8.475730,
      8.166254,  7.115865, 6.364551,  6.076659,  4.936654,  4.815760,  5.000790,
      4.405922,  4.621510, 3.612532,  3.861247,  3.480776,  4.382934,  6.908732,
      5.677212,  5.697563, 5.049410,  5.033481,  5.036818,  4.766425,  4.108461,
      4.493277,  3.970064, 3.479031,  3.105488,  3.182802,  2.832481,  2.904569,
      13.393956, 12.160363, 12.604881,  9.282701,  8.971699,  7.872894,
      7.341750,  6.129950,  6.754877, 6.236324,  5.221481,  4.712401,  4.251457,
      4.107650,  3.642404,  3.409108, 3.673586,  3.573048,  3.006202,  3.574850,
      2.984541,  2.942794,  6.142817, 6.239127,  5.607103,  5.781841,  5.365100,
      4.483500,  4.318313,  3.729648, 4.186974,  3.661007,  3.587046,  3.420915,
      3.554909,  2.951518,  3.197826, 2.768718,  4.050486,  6.375893,  5.984611,
      4.975292,  4.758750,  5.006260, 4.137620,  3.676703,  4.133922,  3.445431,
      3.813060,  3.304948,  3.308399, 3.362024,  3.239059,  5.554801,  5.157033,
      5.049831,  5.366527,  4.410822, 4.798887,  4.386467,  4.013058,  3.749040,
      3.580090,  3.434933,  3.086843, 3.471470,  3.334008,  2.913723,  2.980340,
      2.862694,  5.891940,  6.136203, 6.251388,  6.170612,  5.191264,  4.810252,
      4.989489,  3.951583,  3.793401, 4.320815, 10.128100,  9.961510,  8.155648 };

  PatientData pdone(time, conc);
  PatientPriors ppsingle(1.5, 100, 45, 100, 3.5, 100, 30, 100,
                         10, 100, 1000, 1000, 12, 0, 40);
  PatientEstimates pesingle(2.6, 45, 0.05, 3.5, 30, 12, 10, 10);
  PatientData * data = &pdone;
  PatientPriors * priors = &ppsingle;
  PatientEstimates * estimates = &pesingle;

  Patient pat(data, priors, estimates);

  SECTION( "Estimates can be accessed" ) {
    REQUIRE(pat.estimates->baseline_halflife(0) == 2.6);
    REQUIRE(pat.estimates->mass_mean == 3.5);
    REQUIRE(pat.estimates->pulse_count == 1);
    REQUIRE(pat.estimates->mass_sd == 10);
  }

  SECTION( "Estimates can be updated" ) {
    pat.estimates->baseline_halflife(0) = 10;
    REQUIRE(pat.estimates->baseline_halflife(0) == 10);
    pat.estimates->mass_mean = 5.0;
    REQUIRE(pat.estimates->mass_mean == 5.0);
    pat.estimates->pulse_count = 6;
    REQUIRE(pat.estimates->pulse_count == 6);
    pat.estimates->mass_sd = 50;
    REQUIRE(pat.estimates->mass_sd == 50);
    REQUIRE(pat.estimates->get_decay() == (log(2)/45));
  }

  SECTION( "Priors can be accessed" ) {
    REQUIRE(pat.priors->baseline_mean == 1.5);
    REQUIRE(pat.priors->mass_sd_max == 10);
    REQUIRE(pat.priors->error_alpha == 1000);
    REQUIRE(pat.priors->num_orderstat == 3);
    REQUIRE(pat.priors->strauss_repulsion == 0);
  }

  SECTION( "Priors can be updated" ) {
    pat.priors->baseline_mean     = 2.25;
    pat.priors->mass_sd_max       = 90;
    pat.priors->error_alpha       = 700;
    pat.priors->num_orderstat     = 4;
    pat.priors->strauss_repulsion = 0.75;
    REQUIRE(pat.priors->baseline_mean     == 2.25);
    REQUIRE(pat.priors->mass_sd_max       == 90);
    REQUIRE(pat.priors->error_alpha       == 700);
    REQUIRE(pat.priors->num_orderstat     == 4);
    REQUIRE(pat.priors->strauss_repulsion == 0.75);
  }

  SECTION( "Data can be accessed" ) {
    REQUIRE(pat.data->time(1)              == 20);
    REQUIRE(pat.data->time(143)            == 1440);
    REQUIRE(pat.data->concentration(1)     == log(5.156800));
    REQUIRE(pat.data->concentration(143)   == log(8.155648));
    REQUIRE(pat.data->time.size()          == 144);
    REQUIRE(pat.data->concentration.size() == 144);
    REQUIRE(pat.data->response_concentration.size() == 0);
    REQUIRE(pat.data->avg_period_of_obs    == 10);
    REQUIRE(pat.data->duration_of_obs      == 1430);
    REQUIRE(pat.data->number_of_obs        == 144);
    REQUIRE(pat.data->fitstart             == -40);
    REQUIRE(pat.data->fitend               == 1460);
  }

  // TODO: This section is causing trouble (initial pulse values), not giving
  // sensible results...
  SECTION( "Has initial pulse" ) {
    std::list<PulseEstimates>::const_iterator this_piter = pat.pulses.begin();
    REQUIRE(pat.get_pulsecount()        == 1);
    REQUIRE(this_piter->time            == pat.data->fitstart);
    REQUIRE(this_piter->mass            == 1);
    REQUIRE(this_piter->width           == 1);
    REQUIRE(this_piter->tvarscale_mass  == 1);
    REQUIRE(this_piter->tvarscale_width == 1);
  }

  // Delete initial pulse for these tests
  ++pat.piter;
  pat.piter = pat.pulses.erase(pat.piter);

  // Now add a bunch of pulses (from matt's sim data)
  arma::vec location        { -37.65204, 41.04917, 125.58297, 306.32966,
                              461.65469, 616.95867, 835.64466, 1000.29703,
                              1149.08360, 1319.19888, 1414.31830 };
  arma::vec mass            { 6.635627,  5.940251, 13.463474,  7.649251,
                              4.309497, 12.285200,  4.241401,  3.734942,
                              3.553504,  4.589623, 7.364310 };
  arma::vec width           { 4.596894, 6.740636, 5.986173, 4.476816, 4.701983,
                              3.829487, 5.917652, 4.932739, 6.367862, 4.915543,
                              5.692852 };
  arma::vec tvarscale_mass  { 1.3362258, 0.9683636, 0.2660084, 0.1090139,
                              0.5613391, 0.4339317, 0.3184480, 0.2060616,
                              1.4164371, 2.1759781, 1.0744991 };
  arma::vec tvarscale_width { 1.0301410, 0.7512098, 1.9517880, 1.9834742,
                              1.2743407, 0.4760933, 1.7499186, 0.6057495,
                              1.5219036, 0.8946568, 0.8632652 };

  for (int i = 0; i < location.n_elem; i++) {
    PulseEstimates pulse(location(i), mass(i), width(i), tvarscale_mass(i), tvarscale_width(i),
                        pat.estimates->get_decay(), pat.data->time);
    pat.pulses.push_back(pulse);
  }

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

  PatientData pd(time, conc);
  PatientEstimates pep(2.6, 45, 0.05, 3.5, 30, 12); // population constructor
  PatientData * data = &pd;
  PatientEstimates * estimates = &pep;

  Patient pat(data, estimates);

  SECTION( "Estimates can be accessed" ) {
    REQUIRE(pat.estimates->baseline_halflife(0) == 2.6);
    REQUIRE(pat.estimates->mass_mean == 3.5);
    REQUIRE(pat.estimates->pulse_count == 1);
  }

  SECTION( "Estimates can be updated" ) {
    pat.estimates->baseline_halflife(0) = 10;
    REQUIRE(pat.estimates->baseline_halflife(0) == 10);
    pat.estimates->mass_mean = 5.0;
    REQUIRE(pat.estimates->mass_mean == 5.0);
    pat.estimates->pulse_count = 6;
    REQUIRE(pat.estimates->pulse_count == 6);
  }

  SECTION( "Data can be accessed" ) {
    REQUIRE(pat.data->time(1) == 20);
    REQUIRE(pat.data->time(143) == 1440);
    REQUIRE(pat.data->concentration(1) < 20);
    REQUIRE(pat.data->concentration(1) > 0);
    REQUIRE(pat.data->concentration(143) < 20);
    REQUIRE(pat.data->concentration(143) > 0);
    REQUIRE(pat.data->time.size() == 144);
    REQUIRE(pat.data->concentration.size() == 144);
    REQUIRE(pat.data->response_concentration.size() == 0);
  }

  // TODO: look at what is necessary for the pop constructor version
  //SECTION( "Can add a pulse" ) {
  //}
  //SECTION( "Can remove a pulse" ) {
  //}

}








