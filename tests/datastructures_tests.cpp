#include <RcppArmadillo.h>
#include <RInside.h>
#include "patient.h"
#include "datastructures.h"
#include "utils.h"
#include "proposalvariance.h"
#include "catch.h"

//
// datastructures_tests.cpp
//  For now, also contains tests of utils.h (first set of tests)
//


//
// Test utility functions
//

TEST_CASE( "orderstat_default", "[utils]" ) {

  SECTION( "equal to 3" ) {

    REQUIRE(pulseutils::orderstat_default() == 3);

  }

}

TEST_CASE( "rmvnorm Function", "[utils]" ) {

  RInside R;

  arma::vec initial_means = { 2, 3 };
  arma::vec initial_pvs = { 0.7, 0.1 };
  ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);

  SECTION( "can generate random value" ) {

    arma::vec answer = { 1.3264, 3.1485 };

    pulseutils::set_seed(171227);
    REQUIRE(
            arma::approx_equal(pulseutils::rmvnorm(initial_means, pv.getpsd()),
                               answer, "absdiff", 0.0001) 
            );
  }

}




//
// Test datastructures structs
//


//
// PatientPriors based objects
//
TEST_CASE( "PatientPriors base constructor works (never otherwise used directly)", 
           "[datastructures]" ) {

  PatientPriors ppriors(1.5, 100, 45, 100, 3.5, 100, 30, 100);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    // Test use pointers just for fun.
    //PatientPriors * pppointer = &ppriors;

    REQUIRE(ppriors.baseline_mean     == 1.5);
    REQUIRE(ppriors.baseline_variance == 100.0);
    REQUIRE(ppriors.halflife_mean     == 45.0);
    REQUIRE(ppriors.halflife_variance == 100.0);
    REQUIRE(ppriors.mass_mean         == 3.5);
    REQUIRE(ppriors.mass_variance     == 100.0);
    REQUIRE(ppriors.width_mean        == 30.0);
    REQUIRE(ppriors.width_variance    == 100.0);

  }

}


TEST_CASE( "PopulationEstimates (PatientPriors Population) constructor works", 
           "[datastructures]" ) {

  PopulationEstimates pppop(1.5, 100, 45, 100, 3.5, 100, 30, 100, 5, 10);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(pppop.baseline_mean     == 1.5);
    REQUIRE(pppop.baseline_variance == 100.0);
    REQUIRE(pppop.halflife_mean     == 45.0);
    REQUIRE(pppop.halflife_variance == 100.0);
    REQUIRE(pppop.mass_mean         == 3.5);
    REQUIRE(pppop.mass_variance     == 100.0);
    REQUIRE(pppop.width_mean        == 30.0);
    REQUIRE(pppop.width_variance    == 100.0);
    REQUIRE(pppop.mass_mean_sd      == 5.0);
    REQUIRE(pppop.width_mean_sd     == 10.0);

  }

}


TEST_CASE( "PatientPriors single-subject constructor works",
           "[datastructures]" ) {

  PatientPriors_Single ppsingle(1.5, 100, 45, 100, 3.5, 100, 30, 100,
                                10, 100, 1000, 1000, 12, 0, 40);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(ppsingle.baseline_mean     == 1.5);
    REQUIRE(ppsingle.baseline_variance == 100.0);
    REQUIRE(ppsingle.halflife_mean     == 45.0);
    REQUIRE(ppsingle.halflife_variance == 100.0);
    REQUIRE(ppsingle.mass_mean         == 3.5);
    REQUIRE(ppsingle.mass_variance     == 100.0);
    REQUIRE(ppsingle.width_mean        == 30.0);
    REQUIRE(ppsingle.width_variance    == 100.0);
    REQUIRE(ppsingle.mass_sd_max       == 10);
    REQUIRE(ppsingle.width_sd_max      == 100);
    REQUIRE(ppsingle.error_alpha       == 1000);
    REQUIRE(ppsingle.error_beta        == 1000);
    REQUIRE(ppsingle.pulse_count       == 12);
    REQUIRE(ppsingle.strauss_repulsion == 0);
    REQUIRE(ppsingle.strauss_repulsion_range == 40);

  }

}


//
// PatientEstimates based objects
//
TEST_CASE( "PatientEstimates_Pop population constructor works",
           "[datastructures]" ) {

  PatientEstimates_Pop pepop(3, 45, 0.05, 3.5, 30, 12);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(pepop.baseline    == 3);
    REQUIRE(pepop.halflife    == 45);
    REQUIRE(pepop.errorsq     == 0.05);
    REQUIRE(pepop.mass_mean   == 3.5);
    REQUIRE(pepop.width_mean  == 30);
    REQUIRE(pepop.pulse_count == 12);

  }

  SECTION( "Accessor methods do correct calculations." ) {

    // Change values used in the calculations
    pepop.halflife = 75;
    pepop.errorsq = 0.01;

    REQUIRE(pepop.get_decay() == (log(2) / pepop.halflife));
    REQUIRE(pepop.get_logerrorsq() == log(pepop.errorsq));

  }

}


TEST_CASE( "PatientEstimates_Single single-subject constructor works",
           "[datastructures]" ) {

  PatientEstimates_Single pesingle(3, 45, 0.05, 3.5, 30, 12, 10, 10);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(pesingle.baseline    == 3);
    REQUIRE(pesingle.halflife    == 45);
    REQUIRE(pesingle.errorsq     == 0.05);
    REQUIRE(pesingle.mass_mean   == 3.5);
    REQUIRE(pesingle.width_mean  == 30);
    REQUIRE(pesingle.pulse_count == 12);
    REQUIRE(pesingle.mass_sd  == 10);
    REQUIRE(pesingle.width_sd == 10);

  }

  SECTION( "Accessor methods do correct calculations." ) {

    // Change values used in the calculations
    pesingle.halflife = 75;
    pesingle.errorsq = 0.01;

    REQUIRE(pesingle.get_decay() == (log(2) / pesingle.halflife));
    REQUIRE(pesingle.get_logerrorsq() == log(pesingle.errorsq));

  }

}


//
// PatientData object
//
TEST_CASE( "PatientData single hormone constructor works",
           "[datastructures]" ) {

  NumericVector time(144);
  NumericVector conc = rnorm(144, 3, 0.1);
  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;

  PatientData pdone(time, conc);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(arma::approx_equal(pdone.time, as<arma::vec>(time),
                               "absdiff", 0.00001));
    REQUIRE(arma::approx_equal(pdone.concentration, as<arma::vec>(conc),
                               "absdiff", 0.00001));
  }

  SECTION( "Check calculated variables" ) {

    REQUIRE(pdone.number_of_obs == time.size());
    REQUIRE(pdone.number_of_obs == 144);
    REQUIRE(pdone.duration_of_obs == 1430);
    REQUIRE(pdone.avg_period_of_obs == 10);

  }

}


TEST_CASE( "PatientData two-hormone constructor works",
           "[datastructures]" ) {

  NumericVector time(144);
  NumericVector conc = rnorm(144, 3, 0.1);
  NumericVector responseconc = conc + rnorm(144, 3, 0.1);
  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;

  PatientData pdtwo(time, conc, responseconc);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(arma::approx_equal(pdtwo.time, as<arma::vec>(time),
                               "absdiff", 0.00001));
    REQUIRE(arma::approx_equal(pdtwo.concentration, as<arma::vec>(conc),
                               "absdiff", 0.00001));
    REQUIRE(arma::approx_equal(pdtwo.response_concentration,
                               as<arma::vec>(responseconc),
                               "absdiff", 0.00001));
  }

  SECTION( "Check calculated variables" ) {

    REQUIRE(pdtwo.number_of_obs == time.size());
    REQUIRE(pdtwo.number_of_obs == 144);
    REQUIRE(pdtwo.duration_of_obs == 1430);
    REQUIRE(pdtwo.avg_period_of_obs == 10);

  }

}


//
// PulseEstimate object
//
TEST_CASE( "PulseEstimate works" , "[datastructures]" ) {

  arma::vec conc =
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
  double time = 10.7;
  double mass = 5.1;
  double width = 32.3;
  double tvarscale_mass = 0.174;
  double tvarscale_width = 0.764;
  //double lambda = 0.5;
  double decay_rate = 0.5;
  PulseEstimate pulse(time, mass, width, tvarscale_mass, tvarscale_width,
                      decay_rate, conc);

  SECTION( "member variables can be access" ) {
    REQUIRE(pulse.time == time);
    REQUIRE(pulse.mass == mass);
    REQUIRE(pulse.width == width);
    REQUIRE(pulse.tvarscale_mass == tvarscale_mass);
    REQUIRE(pulse.tvarscale_width == tvarscale_width);
  }

  arma::vec  mc = 
  { 0.3175083, 0.3121017, 0.3538967, 0.3950185, 0.6154131, 0.5570913, 0.5481228,
    0.4512207, 0.3740260, 0.3970621, 0.3384259, 0.3426448, 0.4615028, 0.6297887,
    0.6442114, 0.6584116, 0.6245282, 0.5448004, 0.5029293, 0.4745658, 0.4679122,
    0.3902477, 0.3792236, 0.3775567, 0.2752629, 0.2879139, 0.2466672, 0.2734754,
    0.2297571, 0.2168345, 0.6277222, 0.5763390, 0.5383817, 0.5186881, 0.4476725,
    0.3950478, 0.3749196, 0.2976010, 0.2897507, 0.3017990, 0.2638002, 0.2773183,
    0.2169076, 0.2310946, 0.2095932, 0.2623770, 0.4332139, 0.3472793, 0.3486763,
    0.3049962, 0.3039473, 0.3041669, 0.2865715, 0.2456700, 0.2692408, 0.2374534,
    0.2094972, 0.1895487, 0.1935809, 0.1757215, 0.1793099, 0.6728469, 0.6787910,
    0.6799961, 0.5856320, 0.5682158, 0.4994018, 0.4633291, 0.3786359, 0.4224346,
    0.3860683, 0.3164091, 0.2831069, 0.2543074, 0.2456214, 0.2185856, 0.2056743,
    0.2203447, 0.2147008, 0.1844456, 0.2148012, 0.1833435, 0.1812310, 0.3795339,
    0.3862643, 0.3424772, 0.3544762, 0.3260439, 0.2686294, 0.2583958, 0.2235270,
    0.2503941, 0.2196341, 0.2154817, 0.2063170, 0.2136913, 0.1816712, 0.1943704,
    0.1725852, 0.2422106, 0.3958428, 0.3685141, 0.3001274, 0.2860782, 0.3021581,
    0.2474193, 0.2205210, 0.2471971, 0.2076552, 0.2283078, 0.2000547, 0.2002395,
    0.2031228, 0.1965469, 0.3389062, 0.3121172, 0.3050239, 0.3261401, 0.2641040,
    0.2886618, 0.2625955, 0.2399906, 0.2246336, 0.2150935, 0.2070815, 0.1885839,
    0.2090819, 0.2016134, 0.1797688, 0.1831302, 0.1772199, 0.3620849, 0.3790722,
    0.3871222, 0.3814748, 0.3143943, 0.2893951, 0.3010577, 0.2363670, 0.2271762,
    0.2585494, 0.6268234, 0.6194906, 0.5180005 };

  SECTION( "mean_contribution is working on initialization" ) {
    arma::vec pulsemc = pulse.get_mean_contribution(conc, decay_rate);
    REQUIRE(conc.n_elem == 144);
    REQUIRE(pulsemc.n_elem == 144);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001));
  }

  SECTION( "mean_contribution changes with new decay rate" ) {
    arma::vec pulsemc = pulse.get_mean_contribution(conc, 0.1);
    REQUIRE(pulsemc.n_elem == 144);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == false);
  }

  SECTION( "mean_contribution changes with new pulse time" ) {
    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(conc, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == true);

    // now it should be updated.
    pulse.time = 12.1;
    pulsemc    = pulse.get_mean_contribution(conc, decay_rate);
    REQUIRE(pulse.time == 12.1);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.01) == false);
  }

  SECTION( "mean_contribution changes with new pulse mass" ) {
    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(conc, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == true);

    // now it should be updated.
    pulse.mass = 2.1;
    pulsemc    = pulse.get_mean_contribution(conc, decay_rate);
    REQUIRE(pulse.time == time);
    REQUIRE(pulse.mass == 2.1);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == false);
  }

  SECTION( "mean_contribution changes with new pulse width" ) {
    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(conc, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == true);

    // now it should be updated.
    pulse.width = 42.1;
    pulsemc     = pulse.get_mean_contribution(conc, decay_rate);
    REQUIRE(pulse.mass == mass);
    REQUIRE(pulse.width == 42.1);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == false);
  }


  //SECTION( "mean_contribution changes with new decay rate" ) {
  //  arma::vec pulsemc = pulse.get_mean_contribution(conc, decay_rate);
  //  REQUIRE(conc.n_elem == 144);
  //  REQUIRE(pulsemc.n_elem == 144);
  //  REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == false);
  //}

}




//
// Patient class tests
//
// note: haven't added tests of response_concentration data
//


// Single-subject patient (OnlyPatient)
TEST_CASE( "OnlyPatient class constructor works", "[patient]" ) {

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
  PatientPriors_Single ppsingle(1.5, 100, 45, 100, 3.5, 100, 30, 100,
                                10, 100, 1000, 1000, 12, 0, 40);
  PatientEstimates_Single pesingle(3, 45, 0.05, 3.5, 30, 12, 10, 10);
  PatientData * data = &pdone;
  PatientPriors_Single * priors = &ppsingle;
  PatientEstimates_Single * estimates = &pesingle;

  OnlyPatient pat(data, priors, estimates);

  SECTION( "Estimates can be accessed" ) {
    REQUIRE(pat.estimates->baseline == 3);
    REQUIRE(pat.estimates->mass_mean == 3.5);
    REQUIRE(pat.estimates->pulse_count == 12);
    REQUIRE(pat.estimates->mass_sd == 10);
  }

  SECTION( "Estimates can be updated" ) {
    pat.estimates->baseline = 10;
    REQUIRE(pat.estimates->baseline == 10);
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
    REQUIRE(pat.data->time(1) == 20);
    REQUIRE(pat.data->time(143) == 1440);
    REQUIRE(pat.data->concentration(1) == 5.156800);
    REQUIRE(pat.data->concentration(143) == 8.155648);
    REQUIRE(pat.data->time.size() == 144);
    REQUIRE(pat.data->concentration.size() == 144);
    REQUIRE(pat.data->response_concentration.size() == 0);
  }

  // TODO: look at how the first one is created in the pulsatile() pkg code

  for (int i = 0; i < 12; i++) {
    double time = 10.7 + (i*100);
    double mass = 5.1 + Rf_rnorm(0, 1);
    double width = 32.3 + Rf_rnorm(0, 5);
    double tvarscale_mass = 0.174;
    double tvarscale_width = 0.764;
    //double lambda = 0.5;
    PulseEstimate pulse(time, mass, width, tvarscale_mass, tvarscale_width,
                        pat.estimates->get_decay(), pat.data->concentration);
    pat.pulses.push_back(pulse);
  }

  SECTION( "Can add pulses and iterate with iterators" ) {
    ++pat.piter;
    REQUIRE(pat.get_pulsecount() == 12);
    REQUIRE(pat.piter->time == 10.7);
    ++pat.piter;
    REQUIRE(pat.piter->time == 110.7);
    ++pat.piter;
    REQUIRE(pat.piter->time == 210.7);
    ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter;
    ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter;
    REQUIRE(pat.piter->time == 1110.7);
  }

  SECTION( "Can remove a pulse" ) {
    ++pat.piter; ++pat.piter;
    pat.piter = pat.pulses.erase(pat.piter); // this is how you delete and keep iter correct
    REQUIRE(pat.get_pulsecount() == 11);
    REQUIRE(pat.piter->time == 210.7);
  }

}


// Population patient (OnePatient)
TEST_CASE( "OnePatient class constructor works", "[patient]" ) {

  NumericVector time(144);
  NumericVector conc = rnorm(144, 3, 0.1);
  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;

  PatientData pd(time, conc);
  PatientEstimates_Pop pep(3, 45, 0.05, 3.5, 30, 12);
  PatientData * data = &pd;
  PatientEstimates_Pop * estimates = &pep;

  OnePatient pat(data, estimates);

  SECTION( "Estimates can be accessed" ) {
    REQUIRE(pat.estimates->baseline == 3);
    REQUIRE(pat.estimates->mass_mean == 3.5);
    REQUIRE(pat.estimates->pulse_count == 12);
  }

  SECTION( "Estimates can be updated" ) {
    pat.estimates->baseline = 10;
    REQUIRE(pat.estimates->baseline == 10);
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

  // TODO: look at how the first one is created in the pulsatile() pkg code
  // TODO: Look at iterators for this
  SECTION( "Can add a pulse" ) {
  }

  SECTION( "Can remove a pulse" ) {
  }

}








