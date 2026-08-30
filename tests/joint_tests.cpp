#include <RcppArmadillo.h>
#include <cmath>
#include <bp_datastructures/bp_datastructures.h>
#include <bp_datastructures/associationparameters.h>
#include <bp_mcmc/utils.h>
#include <bpmod_joint/joint_update_lambda.h>
#include <bpmod_joint/joint_birthdeath.h>
#include <bpmod_joint/joint_draw_association.h>
// Note: joint_mcmc_iteration.h not included to avoid linker errors from
// bpmod_singlesubject.h having non-inline functions in headers
// #include <bpmod_joint/joint_mcmc_iteration.h>
#include <testing/catch.h>

using namespace Rcpp;


//
// Test suite for joint hormone model components
//


TEST_CASE("AssociationEstimates structure", "[joint][association]") {

  SECTION("Default constructor") {
    AssociationEstimates assoc;
    REQUIRE(assoc.get_rho() == 1.0);
    REQUIRE(assoc.get_nu() == 1.0);
    REQUIRE(assoc.get_log_rho() == 0.0);
    REQUIRE(assoc.get_log_nu() == 0.0);
  }

  SECTION("Constructor with values") {
    AssociationEstimates assoc(2.0, 0.5);
    REQUIRE(assoc.get_rho() == 2.0);
    REQUIRE(assoc.get_nu() == 0.5);
    REQUIRE(assoc.get_log_rho() == Approx(std::log(2.0)));
    REQUIRE(assoc.get_log_nu() == Approx(std::log(0.5)));
  }

  SECTION("Setters maintain consistency") {
    AssociationEstimates assoc;

    assoc.set_rho(3.0);
    REQUIRE(assoc.get_rho() == 3.0);
    REQUIRE(assoc.get_log_rho() == Approx(std::log(3.0)));

    assoc.set_nu(0.25);
    REQUIRE(assoc.get_nu() == 0.25);
    REQUIRE(assoc.get_log_nu() == Approx(std::log(0.25)));

    assoc.set_log_rho(1.5);
    REQUIRE(assoc.get_log_rho() == 1.5);
    REQUIRE(assoc.get_rho() == Approx(std::exp(1.5)));

    assoc.set_log_nu(-0.5);
    REQUIRE(assoc.get_log_nu() == -0.5);
    REQUIRE(assoc.get_nu() == Approx(std::exp(-0.5)));
  }

  SECTION("Kernel function") {
    AssociationEstimates assoc(2.0, 1.0);

    // Kernel at same location should be maximum
    double k_same = assoc.kernel(5.0, 5.0);
    REQUIRE(k_same == Approx(2.0 / std::sqrt(2 * M_PI * 1.0)));

    // Kernel should decay with distance
    double k_near = assoc.kernel(5.0, 5.5);
    double k_far = assoc.kernel(5.0, 7.0);
    REQUIRE(k_near > k_far);
    REQUIRE(k_far > 0.0);

    // Kernel should be symmetric
    double k1 = assoc.kernel(5.0, 7.0);
    double k2 = assoc.kernel(7.0, 5.0);
    REQUIRE(k1 == Approx(k2));
  }
}


TEST_CASE("AssociationPriors structure", "[joint][association]") {

  SECTION("Default constructor") {
    AssociationPriors priors;
    // Just check it compiles and constructs
  }

  SECTION("Full constructor") {
    AssociationPriors priors(0.5, 1.0, -0.5, 0.8);
    REQUIRE(priors.log_rho_mean == 0.5);
    REQUIRE(priors.log_rho_var == 1.0);
    REQUIRE(priors.log_nu_mean == -0.5);
    REQUIRE(priors.log_nu_var == 0.8);
  }

  SECTION("Constructor from natural scale") {
    // Create priors from natural scale means and SDs
    AssociationPriors priors = AssociationPriors::from_natural_scale(
      2.0, 0.5,  // rho: mean=2, sd=0.5
      1.0, 0.3   // nu: mean=1, sd=0.3
    );

    // Check that priors were created (exact values depend on log-normal approximation)
    REQUIRE(priors.log_rho_mean < std::log(2.0));  // Should be adjusted downward
    REQUIRE(priors.log_rho_var > 0.0);
    REQUIRE(priors.log_nu_mean < std::log(1.0));
    REQUIRE(priors.log_nu_var > 0.0);
  }
}


TEST_CASE("PulseEstimates with lambda field", "[joint][pulse]") {

  arma::vec data_time = arma::linspace(0, 24, 25);

  SECTION("Constructor with default lambda") {
    PulseEstimates pulse(5.0, 10.0, 1.5, 1.0, 1.0, 0.1, data_time);
    REQUIRE(pulse.lambda == 0.0);
  }

  SECTION("Constructor with specified lambda") {
    PulseEstimates pulse(5.0, 10.0, 1.5, 1.0, 1.0, 0.1, data_time, 2.5);
    REQUIRE(pulse.lambda == 2.5);
  }

  SECTION("get_vector_of_values_with_lambda") {
    PulseEstimates pulse(5.0, 10.0, 1.5, 1.0, 1.0, 0.1, data_time, 3.0);
    arma::rowvec vals = pulse.get_vector_of_values_with_lambda();
    REQUIRE(vals.n_elem == 6);
    REQUIRE(vals(0) == 5.0);   // time
    REQUIRE(vals(1) == 10.0);  // mass
    REQUIRE(vals(2) == 1.5);   // width
    REQUIRE(vals(3) == 1.0);   // tvarscale_mass
    REQUIRE(vals(4) == 1.0);   // tvarscale_width
    REQUIRE(vals(5) == 3.0);   // lambda
  }

  SECTION("Empty pulse constructor") {
    PulseEstimates pulse;
    REQUIRE(pulse.lambda == 0.0);
  }
}


TEST_CASE("Lambda update mechanism", "[joint][lambda]") {

  arma::vec data_time = arma::linspace(0, 24, 25);

  SECTION("Single driver, single response") {
    // Create driver pulse at t=5
    std::list<PulseEstimates> driver_pulses;
    driver_pulses.push_back(PulseEstimates(5.0, 10.0, 1.0, 1.0, 1.0, 0.1, data_time));

    // Create response pulse at t=6
    std::list<PulseEstimates> response_pulses;
    response_pulses.push_back(PulseEstimates(6.0, 8.0, 1.2, 1.0, 1.0, 0.1, data_time));

    // Association parameters
    AssociationEstimates assoc(2.0, 1.0);

    // Update lambda
    update_lambda_values(driver_pulses, response_pulses, assoc);

    // Check lambda was calculated
    double expected_lambda = assoc.kernel(6.0, 5.0);
    REQUIRE(response_pulses.front().lambda == Approx(expected_lambda));
    REQUIRE(response_pulses.front().lambda > 0.0);
  }

  SECTION("Multiple drivers, single response") {
    std::list<PulseEstimates> driver_pulses;
    driver_pulses.push_back(PulseEstimates(3.0, 10.0, 1.0, 1.0, 1.0, 0.1, data_time));
    driver_pulses.push_back(PulseEstimates(7.0, 10.0, 1.0, 1.0, 1.0, 0.1, data_time));

    std::list<PulseEstimates> response_pulses;
    response_pulses.push_back(PulseEstimates(5.0, 8.0, 1.2, 1.0, 1.0, 0.1, data_time));

    AssociationEstimates assoc(2.0, 1.0);
    update_lambda_values(driver_pulses, response_pulses, assoc);

    // Lambda should be sum of contributions from both drivers
    double expected_lambda = assoc.kernel(5.0, 3.0) + assoc.kernel(5.0, 7.0);
    REQUIRE(response_pulses.front().lambda == Approx(expected_lambda));
  }

  SECTION("Single driver, multiple responses") {
    std::list<PulseEstimates> driver_pulses;
    driver_pulses.push_back(PulseEstimates(5.0, 10.0, 1.0, 1.0, 1.0, 0.1, data_time));

    std::list<PulseEstimates> response_pulses;
    response_pulses.push_back(PulseEstimates(4.0, 8.0, 1.2, 1.0, 1.0, 0.1, data_time));
    response_pulses.push_back(PulseEstimates(6.0, 8.0, 1.2, 1.0, 1.0, 0.1, data_time));

    AssociationEstimates assoc(2.0, 1.0);
    update_lambda_values(driver_pulses, response_pulses, assoc);

    // Both responses should have lambda values
    auto it = response_pulses.begin();
    REQUIRE(it->lambda == Approx(assoc.kernel(4.0, 5.0)));
    ++it;
    REQUIRE(it->lambda == Approx(assoc.kernel(6.0, 5.0)));
  }

  SECTION("No driver pulses") {
    std::list<PulseEstimates> driver_pulses;  // Empty

    std::list<PulseEstimates> response_pulses;
    response_pulses.push_back(PulseEstimates(5.0, 8.0, 1.2, 1.0, 1.0, 0.1, data_time));

    AssociationEstimates assoc(2.0, 1.0);
    update_lambda_values(driver_pulses, response_pulses, assoc);

    // Lambda should be zero with no drivers
    REQUIRE(response_pulses.front().lambda == 0.0);
  }

  SECTION("update_single_lambda") {
    std::list<PulseEstimates> driver_pulses;
    driver_pulses.push_back(PulseEstimates(3.0, 10.0, 1.0, 1.0, 1.0, 0.1, data_time));
    driver_pulses.push_back(PulseEstimates(7.0, 10.0, 1.0, 1.0, 1.0, 0.1, data_time));

    PulseEstimates response_pulse(5.0, 8.0, 1.2, 1.0, 1.0, 0.1, data_time);

    AssociationEstimates assoc(2.0, 1.0);
    update_single_lambda(driver_pulses, response_pulse, assoc);

    double expected_lambda = assoc.kernel(5.0, 3.0) + assoc.kernel(5.0, 7.0);
    REQUIRE(response_pulse.lambda == Approx(expected_lambda));
  }
}


TEST_CASE("JointBirthDeathProcess instantiation", "[joint][birthdeath]") {

  SECTION("Can create joint birth-death object") {
    JointBirthDeathProcess joint_bd;
    // Just verify it compiles and constructs
    REQUIRE(true);
  }
}


//
// JointBirthDeathProcess::add_new_response_pulse -- log-normal draw path
//
// The response-hormone birth branch draws a new pulse's mass/width from the
// response patient's fixed effects. Under response_patient->lognormal_pulses the
// draw is exp(rnorm(mu, sigma/sqrt(kappa))); with gaussian_random_effects the
// kappa is fixed at 1 (no rgamma draw), so add_new_response_pulse consumes
// exactly two rnorm draws, which a separately-seeded reference reproduces
// exactly. This verifies (a) the log-normal branch draws exp() of a normal and
// is strictly positive, and (b) the new pulse's lambda is the driver-coupling
// kernel evaluated at its position.
//
TEST_CASE("joint add_new_response_pulse: log-normal draw is exp(normal) and positive",
          "[joint][birthdeath][lognormal]") {

  PulseUtils pu;

  // Data grid shared by the driver and response patients.
  NumericVector time = NumericVector::create(0, 2, 4, 6, 8);
  NumericVector conc = NumericVector::create(5, 7, 6, 8, 5);
  PatientData data(time, conc);

  PatientPriors priors(3.0, 1.0, 50.0, 100.0, 10.0, 9.0, 1.5, 0.25,
                       3.0, 0.5, 0.001, 0.001, 4, 0.5, 1.0);

  // Single-subject estimates constructor carries mass_sd/width_sd. The mean/SD
  // values here are LOG-scale (the log-normal parameterization).
  PatientEstimates est(2.6, 45.0, 0.005, 1.2, 3.0, 0.5, 0.7);

  Patient driver(data, priors, est);
  Patient response(data, priors, est);

  response.lognormal_pulses        = true;   // draw on the log scale
  response.gaussian_random_effects = true;   // kappa fixed at 1 -> two rnorms
  response.estimates.mass_mean  = 1.2;       // mean of LOG mass
  response.estimates.mass_sd    = 0.5;
  response.estimates.width_mean = 3.0;       // mean of LOG width
  response.estimates.width_sd   = 0.7;

  // A single driver pulse at t = 4 so the new response pulse gets a positive,
  // known coupling lambda from the kernel.
  const double driver_time = 4.0;
  driver.pulses.clear();
  driver.pulses.push_back(
      PulseEstimates(driver_time, 10.0, 1.0, 1.0, 1.0, 0.1, data.time));

  AssociationEstimates assoc(2.0, 1.0);
  JointBirthDeathProcess jbd;

  const double position = 4.0;
  const int n0 = response.get_pulsecount();

  // Reference draws matching the gaussian-kappa log-normal branch.
  double n_mass = 0.0, n_width = 0.0;
  {
    pu.set_seed(13572468);
    Rcpp::RNGScope scope;
    n_mass  = Rf_rnorm(1.2, 0.5);
    n_width = Rf_rnorm(3.0, 0.7);
  }

  pu.set_seed(13572468);
  jbd.add_new_response_pulse(&driver, &response, assoc, position);

  REQUIRE(response.get_pulsecount() == n0 + 1);
  PulseEstimates & p = response.pulses.back();

  SECTION("mass/width are exp() of the underlying normal draws and positive") {
    REQUIRE(p.mass  == Approx(std::exp(n_mass)));
    REQUIRE(p.width == Approx(std::exp(n_width)));
    REQUIRE(p.mass  > 0.0);
    REQUIRE(p.width > 0.0);
    REQUIRE(p.tvarscale_mass  == 1.0);   // gaussian kappa
    REQUIRE(p.tvarscale_width == 1.0);
  }

  SECTION("new pulse carries the driver-coupling lambda at its position") {
    REQUIRE(p.lambda == Approx(assoc.kernel(position, driver_time)));
    REQUIRE(p.lambda > 0.0);
  }
}


TEST_CASE("Association parameter samplers", "[joint][samplers]") {

  // Note: Joint_DrawRho and Joint_DrawNu are abstract classes that can't be
  // instantiated directly - they're only used through the sample() method in
  // the joint MCMC iteration. Testing them directly would require creating
  // concrete derived classes, which is beyond the scope of unit tests.

  SECTION("JointData structure") {
    NumericVector time = NumericVector::create(0, 2, 4, 6, 8);
    NumericVector conc = NumericVector::create(5, 7, 6, 8, 5);

    PatientData data(time, conc);

    // Create priors using the single-subject constructor
    PatientPriors priors(
      3.0,    // baseline_mean
      1.0,    // baseline_variance
      50.0,   // halflife_mean
      100.0,  // halflife_variance
      10.0,   // mass_mean
      9.0,    // mass_variance
      1.5,    // width_mean
      0.25,   // width_variance
      3.0,    // mass_sd_param
      0.5,    // width_sd_param
      0.001,  // error_alpha
      0.001,  // error_beta
      4,      // pulse_count
      0.5,    // strauss_repulsion
      1.0     // strauss_range
    );

    // Create estimates using the constructor
    PatientEstimates est1(3.0, 50.0, 10.0, 1.5, 1.0);

    Patient driver(data, priors, est1);
    Patient response(data, priors, est1);

    AssociationPriors assoc_priors(0.0, 1.0, 0.0, 1.0);

    JointData joint_data(&driver, &response, &assoc_priors);

    REQUIRE(joint_data.driver_patient == &driver);
    REQUIRE(joint_data.response_patient == &response);
    REQUIRE(joint_data.assoc_priors == &assoc_priors);
  }
}

// Note: JointSamplers construction test commented out due to linker errors
// from including bpmod_singlesubject headers (non-inline functions in headers).
// The joint MCMC orchestration is tested through the R package integration tests.
