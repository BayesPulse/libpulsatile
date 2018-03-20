#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/patientdata.h>
#include <bp_datastructures/patientpriors.h>
#include <bp_datastructures/patientestimates.h>
#include <bp_datastructures/pulseestimates.h>
#include <bp_mcmc/utils.h>
#include <bp_mcmc/proposalvariance.h>
#include <testing/catch.h>


//
// datastructures_tests.cpp
//     Test datastructures structs
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

    REQUIRE(pepop.baseline_halflife(0) == 2.6);
    REQUIRE(pepop.baseline_halflife(1) == 45);
    REQUIRE(pepop.errorsq     == 0.05);
    REQUIRE(pepop.mass_mean   == 3.5);
    REQUIRE(pepop.width_mean  == 30);
    REQUIRE(pepop.pulse_count == 1);

  }

  SECTION( "Accessor methods do correct calculations." ) {

    // Change values used in the calculations
    pepop.baseline_halflife(1) = 75;
    pepop.errorsq = 0.01;

    REQUIRE(pepop.get_decay() == (log(2) / pepop.baseline_halflife(1)));
    REQUIRE(pepop.get_logerrorsq() == log(pepop.errorsq));

  }

}


TEST_CASE( "PatientEstimates single-subject constructor works",
           "[datastructures]" ) {

  PatientEstimates pesingle(2.6, 45, 0.05, 3.5, 30, 12, 10, 10);

  SECTION( "Variables included in constructor are initialized as expected." ) {

    REQUIRE(pesingle.baseline_halflife(0) == 2.6);
    REQUIRE(pesingle.baseline_halflife(1) == 45);
    REQUIRE(pesingle.errorsq     == 0.05);
    REQUIRE(pesingle.mass_mean   == 3.5);
    REQUIRE(pesingle.width_mean  == 30);
    REQUIRE(pesingle.pulse_count == 1);
    REQUIRE(pesingle.mass_sd  == 10);
    REQUIRE(pesingle.width_sd == 10);

  }

  SECTION( "Accessor methods do correct calculations." ) {

    // Change values used in the calculations
    pesingle.baseline_halflife(1) = 75;
    pesingle.errorsq = 0.01;

    REQUIRE(pesingle.get_decay() == (log(2) / pesingle.baseline_halflife(1)));
    REQUIRE(pesingle.get_logerrorsq() == log(pesingle.errorsq));

  }

}



//
// PatientData object
//

TEST_CASE( "PatientData single hormone constructor works",
           "[datastructures]" ) {

  // R sesh for calling R functions
  //RInside R;

  Rcpp::NumericVector time(144);
  Rcpp::NumericVector conc = rnorm(144, 3, 0.1);
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

  Rcpp::NumericVector time(144);
  Rcpp::NumericVector conc = rnorm(144, 3, 0.1);
  Rcpp::NumericVector responseconc = conc + rnorm(144, 3, 0.1);
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
// PulseEstimates object
//

TEST_CASE( "PulseEstimates works" , "[datastructures]" ) {

  arma::vec data_time(144);
  for (int i = 0; i < 144; i++) data_time(i) = 10 * i + 10;
  double time = 126;
  double mass = 13.5;
  double width = 5.99;
  double tvarscale_mass = 0.174;
  double tvarscale_width = 0.764;
  //double lambda = 0.5;
  double decay_rate = 0.015;
  PulseEstimates pulse(time, mass, width, tvarscale_mass, tvarscale_width,
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
      0.0000000004, 0.0948776609, 12.0218825294, 10.9502638717, 9.4249795205,
      8.1121550510, 6.9821965584, 6.0096322708, 5.1725384308, 4.4520450858,
      3.8319107167, 3.2981561188, 2.8387492790, 2.4433341475, 2.1029971898,
      1.8100664557, 1.5579386363, 1.3409302110, 1.1541493284, 0.9933855330,
      0.8550148519, 0.7359181030, 0.6334105812, 0.5451815396, 0.4692420998,
      0.4038804181, 0.3476230974, 0.2992019727, 0.2575255245, 0.2216542731,
      0.1907796008, 0.1642055242, 0.1413330044, 0.1216464442, 0.1047020649,
      0.0901179024, 0.0775651974, 0.0667609841, 0.0574617115, 0.0494577534,
      0.0425686829, 0.0366392049, 0.0315356559, 0.0271429906, 0.0233621885,
      0.0201080220, 0.0173071349, 0.0148963891, 0.0128214409, 0.0110355164,
      0.0094983570, 0.0081753117, 0.0070365560, 0.0060564198, 0.0052128089,
      0.0044867062, 0.0038617438, 0.0033238337, 0.0028608502, 0.0024623566,
      0.0021193699, 0.0018241586, 0.0015700679, 0.0013513699, 0.0011631349,
      0.0010011195, 0.0008616715, 0.0007416475, 0.0006383420, 0.0005494260,
      0.0004728954, 0.0004070248, 0.0003503295, 0.0003015314, 0.0002595305,
      0.0002233799, 0.0001922649, 0.0001654839, 0.0001424333, 0.0001225935,
      0.0001055172, 0.0000908195, 0.0000781691, 0.0000672807, 0.0000579091,
      0.0000498428, 0.0000429001, 0.0000369245, 0.0000317812, 0.0000273543,
      0.0000235441, 0.0000202646, 0.0000174419, 0.0000150124, 0.0000129213,
      0.0000111214, 0.0000095723, 0.0000082390, 0.0000070913, 0.0000061036,
      0.0000052534, 0.0000045216, 0.0000038918, 0.0000033497, 0.0000028831,
      0.0000024815, 0.0000021359, 0.0000018384, 0.0000015823, 0.0000013619,
      0.0000011722, 0.0000010089, 0.0000008684, 0.0000007474, 0.0000006433,
      0.0000005537, 0.0000004766, 0.0000004102, 0.0000003531, 0.0000003039,
      0.0000002616, 0.0000002251, 0.0000001938, 0.0000001668, 0.0000001435,
      0.0000001235, 0.0000001063, 0.0000000915, 0.0000000788, 0.0000000678,
      0.0000000584, 0.0000000502, 0.0000000432, 0.0000000372 };

  SECTION( "mean_contribution is working on initialization" ) {

    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(data_time.n_elem == 144);
    REQUIRE(pulsemc.n_elem == 144);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001));

  }

  SECTION( "mean_contribution changes with new decay rate" ) {

    arma::vec pulsemc = pulse.get_mean_contribution(data_time, 0.1);
    REQUIRE(pulsemc.n_elem == 144);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == false);

  }

  SECTION( "mean_contribution changes with new pulse time" ) {

    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == true);

    // now it should be updated.
    pulse.time = 12.1;
    pulsemc    = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(pulse.time == 12.1);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == false);

  }

  SECTION( "mean_contribution changes with new pulse mass" ) {

    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == true);

    // now it should be updated.
    pulse.mass = 5;
    pulsemc    = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(pulse.time == time);
    REQUIRE(pulse.mass == 5);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == false);

  }

  SECTION( "mean_contribution changes with new pulse width" ) {

    // mean contrib still the same
    arma::vec pulsemc = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == true);

    // now it should be updated.
    pulse.width = 10.;
    pulsemc     = pulse.get_mean_contribution(data_time, decay_rate);
    REQUIRE(pulse.mass == mass);
    REQUIRE(pulse.width == 10.);
    REQUIRE(approx_equal(pulsemc, mc, "absdiff", 0.0000000001) == false);

  }

}


