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

  arma::vec mc(100);
  PulseEstimate pulse { 10.7, 5.1, 32.3, 0.174, 0.764, mc };

  SECTION( "member variables can be access" ) {
    REQUIRE(pulse.time == 10.7);
    REQUIRE(pulse.mass == 5.1);
    REQUIRE(pulse.width == 32.3);
    REQUIRE(pulse.tvarscale_mass == 0.174);
    REQUIRE(pulse.tvarscale_width == 0.764);
    REQUIRE(approx_equal(pulse.mean_contribution, mc, "absdiff", 0.0001));
    REQUIRE(pulse.mean_contribution(0) == 0.0);
    REQUIRE(pulse.mean_contribution(99) == 0.0);
    REQUIRE(pulse.mean_contribution.size() == 100);
  }

}


//
// patient_tests.cpp
//
// note: haven't added tests of response_concentration data
//

// Single-subject patient (OnlyPatient)
TEST_CASE( "OnlyPatient class constructor works", "[patient]" ) {

  NumericVector time(144);
  NumericVector conc = rnorm(144, 3, 0.1);
  for (int i = 0; i < time.size(); i++)  time(i) = (i + 1) * 10;

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


