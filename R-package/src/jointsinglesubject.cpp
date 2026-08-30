#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/bp_datastructures.h>
#include <bp_datastructures/associationparameters.h>
#include <bp_datastructures/chains.h>
#include <bpmod_joint/joint_mcmc_iteration.h>

using namespace Rcpp;


//
// jointsinglesubject.cpp
//   Rcpp export function for joint driver-response hormone model
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// This function implements the joint model for two coupled hormones:
//   - Driver hormone (e.g., LH): Standard pulsatile model
//   - Response hormone (e.g., FSH): Pulsatile model with coupling to driver
//
// The coupling is modeled through lambda (λ), which modulates the birth rate
// of response pulses based on driver pulse locations.
//


//' Joint Single Subject Model
//'
//' Fit a joint Bayesian model for two coupled pulsatile hormones measured
//' from a single subject. Driver pulses influence response pulse occurrence
//' through a coupling kernel parameterized by rho and nu.
//'
//' @param driver_concentration Numeric vector of driver hormone concentrations
//' @param driver_time Numeric vector of driver hormone measurement times
//' @param response_concentration Numeric vector of response hormone concentrations
//' @param response_time Numeric vector of response hormone measurement times
//' @param location_prior Character string: "strauss" (only option currently)
//' @param driver_priors List of prior parameters for driver hormone
//' @param response_priors List of prior parameters for response hormone
//' @param association_priors List of prior parameters for coupling (rho, nu)
//' @param proposalvars List of proposal variances for all parameters
//' @param driver_startingvals List of starting values for driver parameters
//' @param response_startingvals List of starting values for response parameters
//' @param association_startingvals List of starting values for coupling parameters
//' @param mcmc_iterations Integer: total MCMC iterations to run
//' @param thin Integer: thinning interval for saving samples
//' @param burnin Integer: number of initial iterations to discard
//' @param verbose Logical: print progress messages
//' @param pv_adjust_iter Integer: interval for proposal variance adaptation
//' @param pv_adjust_max_iter Integer: iteration to stop adaptation
//' @param bivariate_pv_target_ratio Numeric: target acceptance rate for bivariate proposals
//' @param univariate_pv_target_ratio Numeric: target acceptance rate for univariate proposals
//'
//' @return List containing MCMC chains for driver, response, and association parameters
//' @export
// [[Rcpp::export]]
Rcpp::List jointsinglesubject_(Rcpp::NumericVector driver_concentration,
                                Rcpp::NumericVector driver_time,
                                Rcpp::NumericVector response_concentration,
                                Rcpp::NumericVector response_time,
                                Rcpp::CharacterVector location_prior,
                                Rcpp::List driver_priors,
                                Rcpp::List response_priors,
                                Rcpp::List association_priors,
                                Rcpp::List proposalvars,
                                Rcpp::List driver_startingvals,
                                Rcpp::List response_startingvals,
                                Rcpp::List association_startingvals,
                                int mcmc_iterations,
                                int thin,
                                int burnin,
                                bool verbose,
                                int pv_adjust_iter,
                                int pv_adjust_max_iter,
                                double bivariate_pv_target_ratio,
                                double univariate_pv_target_ratio) {

  // Note: no RNGScope here -- RcppExports.cpp already opens one around this
  // .Call, so an inner scope would be a reference-counted no-op. set.seed()
  // reproducibility is handled by that outer scope.

  //
  // Initialize driver hormone patient
  //
  PatientData driver_data(driver_time, driver_concentration);
  PatientPriors driver_patient_priors(
    driver_priors["baseline_mean"],
    driver_priors["baseline_variance"],
    driver_priors["halflife_mean"],
    driver_priors["halflife_variance"],
    driver_priors["mass_mean"],
    driver_priors["mass_variance"],
    driver_priors["width_mean"],
    driver_priors["width_variance"],
    driver_priors["mass_sd_param"],
    driver_priors["width_sd_param"],
    driver_priors["error_alpha"],
    driver_priors["error_beta"],
    driver_priors["pulse_count"],
    driver_priors["strauss_repulsion"],
    driver_priors["strauss_repulsion_range"]);
  PatientEstimates driver_estimates(
    driver_startingvals["baseline"],
    driver_startingvals["halflife"],
    driver_startingvals["errorsq"],
    driver_startingvals["mass_mean"],
    driver_startingvals["width_mean"],
    driver_startingvals["mass_sd"],
    driver_startingvals["width_sd"]);
  Patient driver_patient(driver_data, driver_patient_priors, driver_estimates);

  // Random-effects distribution for the driver: Student-t (default, per-pulse
  // t-scale kappa) or Gaussian (kappa fixed at 1). Optional element of the
  // priors list; absent -> Student-t. Read as a real so logical/numeric both work.
  {
    bool driver_student_t = true;
    if (driver_priors.containsElementNamed("student_t_pulses")) {
      driver_student_t = (Rf_asReal(driver_priors["student_t_pulses"]) != 0.0);
    }
    driver_patient.gaussian_random_effects = !driver_student_t;
  }

  // Pulse random-effects scale (log-normal vs natural-scale truncated-normal) and
  // SD prior (uniform vs half-Cauchy) for the driver. Optional priors-list
  // elements; absent -> false (legacy natural-scale / half-Cauchy). The
  // SS_DrawSDRandomEffects and birth-death read these flags DIRECTLY from the
  // Patient; the mass/width MEAN draws take the log-normal flag at construction
  // (threaded into JointSamplers below).
  bool driver_lognormal = false;
  if (driver_priors.containsElementNamed("lognormal_pulses")) {
    driver_lognormal = (Rf_asReal(driver_priors["lognormal_pulses"]) != 0.0);
  }
  bool driver_uniform_sd = false;
  if (driver_priors.containsElementNamed("uniform_sd_prior")) {
    driver_uniform_sd = (Rf_asReal(driver_priors["uniform_sd_prior"]) != 0.0);
  }
  driver_patient.lognormal_pulses = driver_lognormal;
  driver_patient.uniform_sd_prior = driver_uniform_sd;

  // Optional Uniform(0, .) SD upper bounds. Only consulted when
  // uniform_sd_prior is true; absent -> the PatientPriors default (10.0). Mirror
  // singlesubject.cpp / population.cpp so the user's prior_driver_sd_* bounds are
  // honored instead of the hardcoded constructor default.
  if (driver_priors.containsElementNamed("mass_sd_max")) {
    driver_patient.priors.mass_sd_max = Rf_asReal(driver_priors["mass_sd_max"]);
  }
  if (driver_priors.containsElementNamed("width_sd_max")) {
    driver_patient.priors.width_sd_max = Rf_asReal(driver_priors["width_sd_max"]);
  }

  //
  // Initialize response hormone patient
  //
  PatientData response_data(response_time, response_concentration);
  PatientPriors response_patient_priors(
    response_priors["baseline_mean"],
    response_priors["baseline_variance"],
    response_priors["halflife_mean"],
    response_priors["halflife_variance"],
    response_priors["mass_mean"],
    response_priors["mass_variance"],
    response_priors["width_mean"],
    response_priors["width_variance"],
    response_priors["mass_sd_param"],
    response_priors["width_sd_param"],
    response_priors["error_alpha"],
    response_priors["error_beta"],
    response_priors["pulse_count"],
    response_priors["strauss_repulsion"],
    response_priors["strauss_repulsion_range"]);
  PatientEstimates response_estimates(
    response_startingvals["baseline"],
    response_startingvals["halflife"],
    response_startingvals["errorsq"],
    response_startingvals["mass_mean"],
    response_startingvals["width_mean"],
    response_startingvals["mass_sd"],
    response_startingvals["width_sd"]);
  Patient response_patient(response_data, response_patient_priors, response_estimates);

  // Random-effects distribution for the response (see driver note above).
  {
    bool response_student_t = true;
    if (response_priors.containsElementNamed("student_t_pulses")) {
      response_student_t = (Rf_asReal(response_priors["student_t_pulses"]) != 0.0);
    }
    response_patient.gaussian_random_effects = !response_student_t;
  }

  // Pulse random-effects scale + SD prior for the response (see driver note above).
  bool response_lognormal = false;
  if (response_priors.containsElementNamed("lognormal_pulses")) {
    response_lognormal = (Rf_asReal(response_priors["lognormal_pulses"]) != 0.0);
  }
  bool response_uniform_sd = false;
  if (response_priors.containsElementNamed("uniform_sd_prior")) {
    response_uniform_sd = (Rf_asReal(response_priors["uniform_sd_prior"]) != 0.0);
  }
  response_patient.lognormal_pulses = response_lognormal;
  response_patient.uniform_sd_prior = response_uniform_sd;

  // Optional Uniform(0, .) SD upper bounds for the response (see driver note).
  if (response_priors.containsElementNamed("mass_sd_max")) {
    response_patient.priors.mass_sd_max = Rf_asReal(response_priors["mass_sd_max"]);
  }
  if (response_priors.containsElementNamed("width_sd_max")) {
    response_patient.priors.width_sd_max = Rf_asReal(response_priors["width_sd_max"]);
  }

  //
  // Initialize association parameters
  //

  // Validate association_priors list structure
  if (association_priors.size() != 4) {
    Rcpp::stop("association_priors must have exactly 4 elements (log_rho_mean, log_rho_var, log_nu_mean, log_nu_var)");
  }
  if (!association_priors.containsElementNamed("log_rho_mean") ||
      !association_priors.containsElementNamed("log_rho_var") ||
      !association_priors.containsElementNamed("log_nu_mean") ||
      !association_priors.containsElementNamed("log_nu_var")) {
    Rcpp::stop("association_priors missing required elements: log_rho_mean, log_rho_var, log_nu_mean, log_nu_var");
  }

  // Use name-based access for clarity and safety
  double log_rho_mean = Rcpp::as<double>(association_priors["log_rho_mean"]);
  double log_rho_var = Rcpp::as<double>(association_priors["log_rho_var"]);
  double log_nu_mean = Rcpp::as<double>(association_priors["log_nu_mean"]);
  double log_nu_var = Rcpp::as<double>(association_priors["log_nu_var"]);

  // Validate association_startingvals list structure
  if (association_startingvals.size() != 2) {
    Rcpp::stop("association_startingvals must have exactly 2 elements (rho, nu)");
  }
  if (!association_startingvals.containsElementNamed("rho") ||
      !association_startingvals.containsElementNamed("nu")) {
    Rcpp::stop("association_startingvals missing required elements: rho, nu");
  }

  double sv_rho = Rcpp::as<double>(association_startingvals["rho"]);
  double sv_nu = Rcpp::as<double>(association_startingvals["nu"]);

  AssociationPriors assoc_priors(log_rho_mean, log_rho_var, log_nu_mean, log_nu_var);
  AssociationEstimates assoc_est(sv_rho, sv_nu);

  //
  // Initialize MCMC samplers
  //
  std::string loc_prior_str = Rcpp::as<std::string>(location_prior[0]);
  JointSamplers samplers(proposalvars, pv_adjust_iter, pv_adjust_max_iter,
                        bivariate_pv_target_ratio, univariate_pv_target_ratio,
                        verbose, mcmc_iterations, loc_prior_str,
                        driver_lognormal, response_lognormal);

  //
  // Initialize output chains
  //
  int num_saved = (mcmc_iterations - burnin) / thin;

  // Driver chains
  arma::mat driver_fixed_effects_chain(num_saved, 7, arma::fill::zeros);  // baseline, halflife, mass_mean, width_mean, errorsq, mass_sd, width_sd
  MatrixVector driver_pulse_chains;

  // Response chains
  arma::mat response_fixed_effects_chain(num_saved, 7, arma::fill::zeros);  // baseline, halflife, mass_mean, width_mean, errorsq, mass_sd, width_sd
  MatrixVector response_pulse_chains;

  // Association chains
  arma::mat association_chain(num_saved, 2, arma::fill::zeros);  // rho, nu

  //
  // Initialize lambda values for response pulses
  //
  update_lambda_values(driver_patient.pulses, response_patient.pulses, assoc_est);

  //
  // Main MCMC loop
  //
  if (verbose) Rcpp::Rcout << "Starting MCMC iterations..." << std::endl;

  int save_index = 0;

  for (int iter = 0; iter < mcmc_iterations; iter++) {

    // Run one MCMC iteration
    joint_mcmc_iteration(&driver_patient, &response_patient,
                        &assoc_est, &assoc_priors,
                        samplers, iter);

    // Save samples after burnin, applying thinning
    if (iter >= burnin && (iter - burnin) % thin == 0) {

      // Driver fixed effects (baseline, halflife, mass_mean, width_mean, errorsq)
      driver_fixed_effects_chain(save_index, 0) = driver_patient.estimates.baseline_halflife(0);
      driver_fixed_effects_chain(save_index, 1) = driver_patient.estimates.baseline_halflife(1);
      driver_fixed_effects_chain(save_index, 2) = driver_patient.estimates.mass_mean;
      driver_fixed_effects_chain(save_index, 3) = driver_patient.estimates.width_mean;
      driver_fixed_effects_chain(save_index, 4) = driver_patient.estimates.errorsq;
      driver_fixed_effects_chain(save_index, 5) = driver_patient.estimates.mass_sd;
      driver_fixed_effects_chain(save_index, 6) = driver_patient.estimates.width_sd;

      // Driver pulses
      arma::mat driver_pulses_matrix(driver_patient.get_pulsecount(), 5);
      int pulse_idx = 0;
      for (const auto& pulse : driver_patient.pulses) {
        driver_pulses_matrix.row(pulse_idx++) = pulse.get_vector_of_values();
      }
      driver_pulse_chains.push_back(driver_pulses_matrix);

      // Response fixed effects (baseline, halflife, mass_mean, width_mean, errorsq)
      response_fixed_effects_chain(save_index, 0) = response_patient.estimates.baseline_halflife(0);
      response_fixed_effects_chain(save_index, 1) = response_patient.estimates.baseline_halflife(1);
      response_fixed_effects_chain(save_index, 2) = response_patient.estimates.mass_mean;
      response_fixed_effects_chain(save_index, 3) = response_patient.estimates.width_mean;
      response_fixed_effects_chain(save_index, 4) = response_patient.estimates.errorsq;
      response_fixed_effects_chain(save_index, 5) = response_patient.estimates.mass_sd;
      response_fixed_effects_chain(save_index, 6) = response_patient.estimates.width_sd;

      // Response pulses (including lambda)
      arma::mat response_pulses_matrix(response_patient.get_pulsecount(), 6);
      pulse_idx = 0;
      for (const auto& pulse : response_patient.pulses) {
        response_pulses_matrix.row(pulse_idx++) = pulse.get_vector_of_values_with_lambda();
      }
      response_pulse_chains.push_back(response_pulses_matrix);

      // Association parameters
      association_chain(save_index, 0) = assoc_est.get_rho();
      association_chain(save_index, 1) = assoc_est.get_nu();

      save_index++;
    }

    // Progress reporting
    if (verbose && (iter + 1) % 1000 == 0) {
      Rcpp::Rcout << "Iteration " << (iter + 1) << " / " << mcmc_iterations << std::endl;
    }
  }

  if (verbose) Rcpp::Rcout << "MCMC complete!" << std::endl;

  //
  // Convert MatrixVector to R List (each arma::mat needs to be converted to NumericMatrix)
  //
  Rcpp::List driver_pulse_list(driver_pulse_chains.size());
  for (size_t i = 0; i < driver_pulse_chains.size(); i++) {
    driver_pulse_list[i] = Rcpp::wrap(driver_pulse_chains[i]);
  }

  Rcpp::List response_pulse_list(response_pulse_chains.size());
  for (size_t i = 0; i < response_pulse_chains.size(); i++) {
    response_pulse_list[i] = Rcpp::wrap(response_pulse_chains[i]);
  }

  //
  // Return results
  //
  return Rcpp::List::create(
    Rcpp::Named("driver_fixed_effects") = driver_fixed_effects_chain,
    Rcpp::Named("driver_pulses") = driver_pulse_list,
    Rcpp::Named("response_fixed_effects") = response_fixed_effects_chain,
    Rcpp::Named("response_pulses") = response_pulse_list,
    Rcpp::Named("association") = association_chain,
    Rcpp::Named("driver_colnames") = Rcpp::CharacterVector::create(
      "baseline", "halflife", "mass_mean", "width_mean", "errorsq", "mass_sd", "width_sd"),
    Rcpp::Named("response_colnames") = Rcpp::CharacterVector::create(
      "baseline", "halflife", "mass_mean", "width_mean", "errorsq", "mass_sd", "width_sd"),
    Rcpp::Named("driver_pulse_colnames") = Rcpp::CharacterVector::create(
      "time", "mass", "width", "tvarscale_mass", "tvarscale_width"),
    Rcpp::Named("response_pulse_colnames") = Rcpp::CharacterVector::create(
      "time", "mass", "width", "tvarscale_mass", "tvarscale_width", "lambda"),
    Rcpp::Named("association_colnames") = Rcpp::CharacterVector::create(
      "rho", "nu")
  );
}
