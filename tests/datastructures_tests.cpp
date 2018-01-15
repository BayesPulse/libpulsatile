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

TEST_CASE( "PopulationEstimates (PatientPriors Population) constructor works", 
           "[datastructures]" ) {

  PatientPriors pppop(1.5, 100, 45, 100, 3.5, 100, 30, 100, 5, 10);

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

  PatientPriors ppsingle(1.5, 100, 45, 100, 3.5, 100, 30, 100,
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
TEST_CASE( "PatientEstimates population constructor works",
           "[datastructures]" ) {

  PatientEstimates pepop(2.6, 45, 0.05, 3.5, 30, 12);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(pepop.baseline    == 2.6);
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


TEST_CASE( "PatientEstimates single-subject constructor works",
           "[datastructures]" ) {

  PatientEstimates pesingle(2.6, 45, 0.05, 3.5, 30, 12, 10, 10);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(pesingle.baseline    == 2.6);
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
    REQUIRE(arma::approx_equal(pdone.concentration, log(as<arma::vec>(conc)),
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
    REQUIRE(arma::approx_equal(pdtwo.concentration, 
                               log(as<arma::vec>(conc)),
                               "absdiff", 0.00001));
    REQUIRE(arma::approx_equal(pdtwo.response_concentration,
                               log(as<arma::vec>(responseconc)),
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
////
//TEST_CASE( "PulseEstimate works" , "[datastructures]" ) {
//
//  arma::vec data_time(144);
//  for (int i = 0; i < 144; i++) data_time(i) = 10 * i + 10;
//  double time = 126;
//  double mass = 13.5;
//  double width = 5.99;
//  double tvarscale_mass = 0.174;
//  double tvarscale_width = 0.764;
//  //double lambda = 0.5;
//  double decay_rate = 0.015;
//  PulseEstimate pulse(time, mass, width, tvarscale_mass, tvarscale_width,
//                      decay_rate, data_time);
//
//  SECTION( "member variables can be access" ) {
//    REQUIRE(pulse.time == time);
//    REQUIRE(pulse.mass == mass);
//    REQUIRE(pulse.width == width);
//    REQUIRE(pulse.tvarscale_mass == tvarscale_mass);
//    REQUIRE(pulse.tvarscale_width == tvarscale_width);
//  }
//
//  arma::vec  mc =
//    { 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
//      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
//      0.0000000004, 0.0948776609, 12.0218825294, 10.9502638717, 9.4249795205,
//      8.1121550510, 6.9821965584, 6.0096322708, 5.1725384308, 4.4520450858,
//      3.8319107167, 3.2981561188, 2.8387492790, 2.4433341475, 2.1029971898,
//      1.8100664557, 1.5579386363, 1.3409302110, 1.1541493284, 0.9933855330,
//      0.8550148519, 0.7359181030, 0.6334105812, 0.5451815396, 0.4692420998,
//      0.4038804181, 0.3476230974, 0.2992019727, 0.2575255245, 0.2216542731,
//      0.1907796008, 0.1642055242, 0.1413330044, 0.1216464442, 0.1047020649,
//      0.0901179024, 0.0775651974, 0.0667609841, 0.0574617115, 0.0494577534,
//      0.0425686829, 0.0366392049, 0.0315356559, 0.0271429906, 0.0233621885,
//      0.0201080220, 0.0173071349, 0.0148963891, 0.0128214409, 0.0110355164,
//      0.0094983570, 0.0081753117, 0.0070365560, 0.0060564198, 0.0052128089,
//      0.0044867062, 0.0038617438, 0.0033238337, 0.0028608502, 0.0024623566,
//      0.0021193699, 0.0018241586, 0.0015700679, 0.0013513699, 0.0011631349,
//      0.0010011195, 0.0008616715, 0.0007416475, 0.0006383420, 0.0005494260,
//      0.0004728954, 0.0004070248, 0.0003503295, 0.0003015314, 0.0002595305,
//      0.0002233799, 0.0001922649, 0.0001654839, 0.0001424333, 0.0001225935,
//      0.0001055172, 0.0000908195, 0.0000781691, 0.0000672807, 0.0000579091,
//      0.0000498428, 0.0000429001, 0.0000369245, 0.0000317812, 0.0000273543,
//      0.0000235441, 0.0000202646, 0.0000174419, 0.0000150124, 0.0000129213,
//      0.0000111214, 0.0000095723, 0.0000082390, 0.0000070913, 0.0000061036,
//      0.0000052534, 0.0000045216, 0.0000038918, 0.0000033497, 0.0000028831,
//      0.0000024815, 0.0000021359, 0.0000018384, 0.0000015823, 0.0000013619,
//      0.0000011722, 0.0000010089, 0.0000008684, 0.0000007474, 0.0000006433,
//      0.0000005537, 0.0000004766, 0.0000004102, 0.0000003531, 0.0000003039,
//      0.0000002616, 0.0000002251, 0.0000001938, 0.0000001668, 0.0000001435,
//      0.0000001235, 0.0000001063, 0.0000000915, 0.0000000788, 0.0000000678,
//      0.0000000584, 0.0000000502, 0.0000000432, 0.0000000372 };
//
//  SECTION( "mean_contribution is working on initialization" ) {
//
//    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
//    REQUIRE(data_time.n_elem == 144);
//    REQUIRE(pulsemc.n_elem == 144);
//    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001));
//
//  }
//
//  SECTION( "mean_contribution changes with new decay rate" ) {
//
//    arma::vec pulsemc = pulse.get_mean_contribution(data_time, 0.1);
//    REQUIRE(pulsemc.n_elem == 144);
//    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == false);
//
//  }
//
//  SECTION( "mean_contribution changes with new pulse time" ) {
//
//    // mean contrib still the same
//    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
//    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == true);
//
//    // now it should be updated.
//    pulse.time = 12.1;
//    pulsemc    = pulse.get_mean_contribution(data_time, decay_rate);
//    REQUIRE(pulse.time == 12.1);
//    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == false);
//
//  }
//
//  SECTION( "mean_contribution changes with new pulse mass" ) {
//
//    // mean contrib still the same
//    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
//    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == true);
//
//    // now it should be updated.
//    pulse.mass = 5;
//    pulsemc    = pulse.get_mean_contribution(data_time, decay_rate);
//    REQUIRE(pulse.time == time);
//    REQUIRE(pulse.mass == 5);
//    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == false);
//
//  }
//
//  SECTION( "mean_contribution changes with new pulse width" ) {
//
//    // mean contrib still the same
//    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
//    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == true);
//
//    // now it should be updated.
//    pulse.width = 10.;
//    pulsemc     = pulse.get_mean_contribution(data_time, decay_rate);
//    REQUIRE(pulse.mass == mass);
//    REQUIRE(pulse.width == 10.);
//    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == false);
//
//  }
//
//}
//



//
// Patient class tests
//
// note: haven't added tests of response_concentration data
//


// Single-subject patient (OnlyPatient)
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
    REQUIRE(pat.estimates->baseline == 2.6);
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
    std::list<PulseEstimate>::const_iterator this_piter = pat.pulses.begin();
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
    PulseEstimate pulse(location(i), mass(i), width(i), tvarscale_mass(i), tvarscale_width(i),
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


// Population patient
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
    REQUIRE(pat.estimates->baseline == 2.6);
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

  // TODO: look at what is necessary for the pop constructor version
  //SECTION( "Can add a pulse" ) {
  //}
  //SECTION( "Can remove a pulse" ) {
  //}

}








