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

  PatientEstimates pepop(3, 45, 0.05, 3.5, 30, 12);

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


TEST_CASE( "PatientEstimates single-subject constructor works",
           "[datastructures]" ) {

  PatientEstimates pesingle(3, 45, 0.05, 3.5, 30, 12, 10, 10);

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
//
TEST_CASE( "PulseEstimate works" , "[datastructures]" ) {

  arma::vec data_time(144);
  for (int i = 0; i < 144; i++) data_time(i) = 10 * i + 10;
  double time = 126;
  double mass = 13.5;
  double width = 5.99;
  double tvarscale_mass = 0.174;
  double tvarscale_width = 0.764;
  //double lambda = 0.5;
  double decay_rate = 0.015;
  PulseEstimate pulse(time, mass, width, tvarscale_mass, tvarscale_width,
                      decay_rate, data_time);

  SECTION( "member variables can be access" ) {
    REQUIRE(pulse.time == time);
    REQUIRE(pulse.mass == mass);
    REQUIRE(pulse.width == width);
    REQUIRE(pulse.tvarscale_mass == tvarscale_mass);
    REQUIRE(pulse.tvarscale_width == tvarscale_width);
  }

  arma::vec  mc =
    { 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000001, 0.0182916377, 2.3177207155, 2.1111213950, 1.8170590358,
      1.5639572058, 1.3461104418, 1.1586079944, 0.9972231423, 0.8583179129,
      0.7387610739, 0.6358575490, 0.5472876643, 0.4710548581, 0.4054406737,
      0.3489660218, 0.3003578385, 0.2585203873, 0.2225105595, 0.1915166134,
      0.1648398767, 0.1418789967, 0.1221163842, 0.1051065459, 0.0904660424,
      0.0778648443, 0.0670188926, 0.0576836954, 0.0496488168, 0.0427331326,
      0.0367807481, 0.0316574833, 0.0272478484, 0.0234524404, 0.0201857025,
      0.0173739952, 0.0149539362, 0.0128709722, 0.0110781484, 0.0095350507,
      0.0082068942, 0.0070637393, 0.0060798168, 0.0052329468, 0.0045040390,
      0.0038766623, 0.0033366742, 0.0028719021, 0.0024718690, 0.0021275574,
      0.0018312056, 0.0015761333, 0.0013565905, 0.0011676283, 0.0010049870,
      0.0008650003, 0.0007445126, 0.0006408080, 0.0005515485, 0.0004747222,
      0.0004085972, 0.0003516829, 0.0003026963, 0.0002605331, 0.0002242429,
      0.0001930077, 0.0001661232, 0.0001429836, 0.0001230671, 0.0001059248,
      0.0000911704, 0.0000784711, 0.0000675407, 0.0000581328, 0.0000500354,
      0.0000430658, 0.0000370671, 0.0000319040, 0.0000274600, 0.0000236350,
      0.0000203429, 0.0000175093, 0.0000150704, 0.0000129712, 0.0000111644,
      0.0000096093, 0.0000082708, 0.0000071187, 0.0000061272, 0.0000052737,
      0.0000045391, 0.0000039068, 0.0000033627, 0.0000028943, 0.0000024911,
      0.0000021441, 0.0000018455, 0.0000015884, 0.0000013672, 0.0000011767,
      0.0000010128, 0.0000008717, 0.0000007503, 0.0000006458, 0.0000005558,
      0.0000004784, 0.0000004118, 0.0000003544, 0.0000003051, 0.0000002626,
      0.0000002260, 0.0000001945, 0.0000001674, 0.0000001441, 0.0000001240,
      0.0000001067, 0.0000000919, 0.0000000791, 0.0000000681, 0.0000000586,
      0.0000000504, 0.0000000434, 0.0000000374, 0.0000000322, 0.0000000277,
      0.0000000238, 0.0000000205, 0.0000000176, 0.0000000152, 0.0000000131,
      0.0000000113, 0.0000000097, 0.0000000083, 0.0000000072 };

  SECTION( "mean_contribution is working on initialization" ) {
    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(data_time.n_elem == 144);
    REQUIRE(pulsemc.n_elem == 144);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0001));

  }

  SECTION( "mean_contribution changes with new decay rate" ) {
    arma::vec pulsemc = pulse.get_mean_contribution(data_time, 0.1);
    REQUIRE(pulsemc.n_elem == 144);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == false);
    std::cout << "Internal mc = " << pulsemc << "\n";
    std::cout << "External mc = " << mc << "\n";
  }

  SECTION( "mean_contribution changes with new pulse time" ) {
    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == true);

    // now it should be updated.
    pulse.time = 12.1;
    pulsemc    = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(pulse.time == 12.1);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.01) == false);
  }

  SECTION( "mean_contribution changes with new pulse mass" ) {
    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == true);

    // now it should be updated.
    pulse.mass = 5;
    pulsemc    = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(pulse.time == time);
    REQUIRE(pulse.mass == 5);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == false);
  }

  SECTION( "mean_contribution changes with new pulse width" ) {
    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == true);

    // now it should be updated.
    pulse.width = 10.;
    pulsemc     = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(pulse.mass == mass);
    REQUIRE(pulse.width == 10.);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000001) == false);
  }

}




////
//// Patient class tests
////
//// note: haven't added tests of response_concentration data
////
//
//
//// Single-subject patient (OnlyPatient)
//TEST_CASE( "Patient class constructor for single-subject works", "[patient]" ) {
//
//  NumericVector time(144);
//  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;
//  NumericVector conc =
//    { 5.237937,  5.156800,  5.773434,  6.364134,  9.872640,  8.782145, 8.633523,
//      7.166886,  6.063834,  6.393283,  5.547755,  5.609554, 17.809125,
//      14.864233, 14.488007, 11.012978, 10.074941,  8.579306,  7.925953,
//      7.505312,  7.408293, 6.296038,  6.138372,  6.114483,  4.589039,  4.787283,
//      4.125096,  4.560712, 3.838159,  3.611227, 10.149188,  9.114352,  8.475730,
//      8.166254,  7.115865, 6.364551,  6.076659,  4.936654,  4.815760,  5.000790,
//      4.405922,  4.621510, 3.612532,  3.861247,  3.480776,  4.382934,  6.908732,
//      5.677212,  5.697563, 5.049410,  5.033481,  5.036818,  4.766425,  4.108461,
//      4.493277,  3.970064, 3.479031,  3.105488,  3.182802,  2.832481,  2.904569,
//      13.393956, 12.160363, 12.604881,  9.282701,  8.971699,  7.872894,
//      7.341750,  6.129950,  6.754877, 6.236324,  5.221481,  4.712401,  4.251457,
//      4.107650,  3.642404,  3.409108, 3.673586,  3.573048,  3.006202,  3.574850,
//      2.984541,  2.942794,  6.142817, 6.239127,  5.607103,  5.781841,  5.365100,
//      4.483500,  4.318313,  3.729648, 4.186974,  3.661007,  3.587046,  3.420915,
//      3.554909,  2.951518,  3.197826, 2.768718,  4.050486,  6.375893,  5.984611,
//      4.975292,  4.758750,  5.006260, 4.137620,  3.676703,  4.133922,  3.445431,
//      3.813060,  3.304948,  3.308399, 3.362024,  3.239059,  5.554801,  5.157033,
//      5.049831,  5.366527,  4.410822, 4.798887,  4.386467,  4.013058,  3.749040,
//      3.580090,  3.434933,  3.086843, 3.471470,  3.334008,  2.913723,  2.980340,
//      2.862694,  5.891940,  6.136203, 6.251388,  6.170612,  5.191264,  4.810252,
//      4.989489,  3.951583,  3.793401, 4.320815, 10.128100,  9.961510,  8.155648 };
//
//  PatientData pdone(time, conc);
//  PatientPriors ppsingle(1.5, 100, 45, 100, 3.5, 100, 30, 100,
//                                10, 100, 1000, 1000, 12, 0, 40);
//  PatientEstimates pesingle(3, 45, 0.05, 3.5, 30, 12, 10, 10);
//  PatientData * data = &pdone;
//  PatientPriors * priors = &ppsingle;
//  PatientEstimates * estimates = &pesingle;
//
//  Patient pat(data, priors, estimates);
//
//  SECTION( "Estimates can be accessed" ) {
//    REQUIRE(pat.estimates->baseline == 3);
//    REQUIRE(pat.estimates->mass_mean == 3.5);
//    REQUIRE(pat.estimates->pulse_count == 12);
//    REQUIRE(pat.estimates->mass_sd == 10);
//  }
//
//  SECTION( "Estimates can be updated" ) {
//    pat.estimates->baseline = 10;
//    REQUIRE(pat.estimates->baseline == 10);
//    pat.estimates->mass_mean = 5.0;
//    REQUIRE(pat.estimates->mass_mean == 5.0);
//    pat.estimates->pulse_count = 6;
//    REQUIRE(pat.estimates->pulse_count == 6);
//    pat.estimates->mass_sd = 50;
//    REQUIRE(pat.estimates->mass_sd == 50);
//    REQUIRE(pat.estimates->get_decay() == (log(2)/45));
//  }
//
//  SECTION( "Priors can be accessed" ) {
//    REQUIRE(pat.priors->baseline_mean == 1.5);
//    REQUIRE(pat.priors->mass_sd_max == 10);
//    REQUIRE(pat.priors->error_alpha == 1000);
//    REQUIRE(pat.priors->num_orderstat == 3);
//    REQUIRE(pat.priors->strauss_repulsion == 0);
//  }
//
//  SECTION( "Priors can be updated" ) {
//    pat.priors->baseline_mean     = 2.25;
//    pat.priors->mass_sd_max       = 90;
//    pat.priors->error_alpha       = 700;
//    pat.priors->num_orderstat     = 4;
//    pat.priors->strauss_repulsion = 0.75;
//    REQUIRE(pat.priors->baseline_mean     == 2.25);
//    REQUIRE(pat.priors->mass_sd_max       == 90);
//    REQUIRE(pat.priors->error_alpha       == 700);
//    REQUIRE(pat.priors->num_orderstat     == 4);
//    REQUIRE(pat.priors->strauss_repulsion == 0.75);
//  }
//
//  SECTION( "Data can be accessed" ) {
//    REQUIRE(pat.data->time(1)              == 20);
//    REQUIRE(pat.data->time(143)            == 1440);
//    REQUIRE(pat.data->concentration(1)     == log(5.156800));
//    REQUIRE(pat.data->concentration(143)   == log(8.155648));
//    REQUIRE(pat.data->time.size()          == 144);
//    REQUIRE(pat.data->concentration.size() == 144);
//    REQUIRE(pat.data->response_concentration.size() == 0);
//    REQUIRE(pat.data->avg_period_of_obs    == 10);
//    REQUIRE(pat.data->duration_of_obs      == 1430);
//    REQUIRE(pat.data->number_of_obs        == 144);
//    REQUIRE(pat.data->fitstart             == -40);
//    REQUIRE(pat.data->fitend               == 1460);
//  }
//
//  // TODO: This section is causing trouble (initial pulse values), not giving
//  // sensible results...
//  SECTION( "Has initial pulse" ) {
//    std::list<PulseEstimate>::const_iterator this_piter = pat.pulses.begin();
//    //std::cout << "Pulse count = " << this_piter->time << "\n"; 
//    //std::cout << "time = " << this_piter->time << "\n"; 
//    //std::cout << "mass = " << this_piter->mass << "\n"; 
//    //std::cout << "width = " << this_piter->width << "\n"; 
//    REQUIRE(pat.get_pulsecount()        == 1);
//    REQUIRE(this_piter->time            == pat.data->fitstart);
//    REQUIRE(this_piter->mass            == 1);
//    REQUIRE(this_piter->width           == 1);
//    REQUIRE(this_piter->tvarscale_mass  == 1);
//    REQUIRE(this_piter->tvarscale_width == 1);
//  }
//
//  // Delete initial pulse for these tests
//  ++pat.piter;
//  pat.piter = pat.pulses.erase(pat.piter);
//
//  for (int i = 0; i < 12; i++) {
//    double time = 10.7 + (i*100);
//    double mass = 5.1 + Rf_rnorm(0, 1);
//    double width = 32.3 + Rf_rnorm(0, 5);
//    double tvarscale_mass = 0.174;
//    double tvarscale_width = 0.764;
//    //double lambda = 0.5;
//    PulseEstimate pulse(time, mass, width, tvarscale_mass, tvarscale_width,
//                        pat.estimates->get_decay(), pat.data->concentration);
//    pat.pulses.push_back(pulse);
//    //std::cout << "Pulse #" << i + 1 << "/12 mean contribution\n";
//    //std::cout << pulse.get_mean_contribution(pat.data->concentration, pat.estimates->get_decay()) << "\n";
//  }
//
//  SECTION( "Can add pulses and iterate with iterators" ) {
//    ++pat.piter;
//    REQUIRE(pat.get_pulsecount() == 12);
//    REQUIRE(pat.piter->time == 10.7);
//    ++pat.piter;
//    REQUIRE(pat.piter->time == 110.7);
//    ++pat.piter;
//    REQUIRE(pat.piter->time == 210.7);
//    ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter;
//    ++pat.piter; ++pat.piter; ++pat.piter; ++pat.piter;
//    REQUIRE(pat.piter->time == 1110.7);
//  }
//
//  SECTION( "Can remove a pulse" ) {
//    ++pat.piter; ++pat.piter;
//    pat.piter = pat.pulses.erase(pat.piter); // this is how you delete and keep iter correct
//    REQUIRE(pat.get_pulsecount() == 11);
//    REQUIRE(pat.piter->time == 210.7);
//  }
//
//  SECTION( "Can get mean_concentration" ) {
//
//    //arma::vec mconc_all, mconc_lessone;
//    arma::vec mconc = pat.mean_concentration(false);
//    REQUIRE(mconc.n_elem == 144);
//    //std::cout << "Conc = " << mconc << "\n";
//    REQUIRE(mconc.n_elem == 144);
//    //REQUIRE(mconc > 0);
//    //REQUIRE(mconc < 10);
//
//  }
//
//}
//
//
//// Population patient
//TEST_CASE( "Patient class constructor for population model works", "[patient]" ) {
//
//  NumericVector time(144);
//  NumericVector conc = rnorm(144, 3, 0.1);
//  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;
//
//  PatientData pd(time, conc);
//  PatientEstimates pep(3, 45, 0.05, 3.5, 30, 12); // population constructor
//  PatientData * data = &pd;
//  PatientEstimates * estimates = &pep;
//
//  Patient pat(data, estimates);
//
//  SECTION( "Estimates can be accessed" ) {
//    REQUIRE(pat.estimates->baseline == 3);
//    REQUIRE(pat.estimates->mass_mean == 3.5);
//    REQUIRE(pat.estimates->pulse_count == 12);
//  }
//
//  SECTION( "Estimates can be updated" ) {
//    pat.estimates->baseline = 10;
//    REQUIRE(pat.estimates->baseline == 10);
//    pat.estimates->mass_mean = 5.0;
//    REQUIRE(pat.estimates->mass_mean == 5.0);
//    pat.estimates->pulse_count = 6;
//    REQUIRE(pat.estimates->pulse_count == 6);
//  }
//
//  SECTION( "Data can be accessed" ) {
//    REQUIRE(pat.data->time(1) == 20);
//    REQUIRE(pat.data->time(143) == 1440);
//    REQUIRE(pat.data->concentration(1) < 20);
//    REQUIRE(pat.data->concentration(1) > 0);
//    REQUIRE(pat.data->concentration(143) < 20);
//    REQUIRE(pat.data->concentration(143) > 0);
//    REQUIRE(pat.data->time.size() == 144);
//    REQUIRE(pat.data->concentration.size() == 144);
//    REQUIRE(pat.data->response_concentration.size() == 0);
//  }
//
//  // TODO: look at how the first one is created in the pulsatile() pkg code
//  // TODO: Look at iterators for this
//  SECTION( "Can add a pulse" ) {
//  }
//
//  SECTION( "Can remove a pulse" ) {
//  }
//
//}
//
//
//
//




