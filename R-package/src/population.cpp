#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <bpmod_population/pop_mcmc_iteration.h>
#include <bp_datastructures/population.h>
#include <bp_datastructures/populationchains.h>
#include <iostream>
#ifndef NORINSIDE
#include <RInside.h>
#endif

using namespace Rcpp;

//
// population.cpp
//   Implementation of the population model. Creates Population object
//   (containing multiple subjects with data, priors, estimates), the samplers
//   (Metropolis Hastings, Gibbs, birth death process), and output routines
//   (chains, verbose output/diagnostics). Intended to be called by
//   fit_pulse_population() in the R package.
//
// Author: Matt Mulvahill
// Created: 12/30/24
//

// [[Rcpp::export]]
Rcpp::List population_(Rcpp::List subject_data_list,
                       Rcpp::CharacterVector location_prior,
                       Rcpp::List population_priors,
                       Rcpp::List proposalvars,
                       Rcpp::List population_startingvals,
                       Rcpp::List subject_startingvals_list,
                       int mcmc_iterations,
                       int thin,
                       int burnin,
                       bool verbose,
                       int pv_adjust_iter,
                       int pv_adjust_max_iter,
                       double bivariate_pv_target_ratio,
                       double univariate_pv_target_ratio)
{

  // Note: no RNGScope here -- RcppExports.cpp already opens one around this
  // .Call, so an inner scope would be a reference-counted no-op. set.seed()
  // reproducibility is handled by that outer scope.

  // Every nth iteration for printing verbose screen output
  int verbose_iter = 5000;

  // Input validation
  if (!population_priors.inherits("bp_population_priors"))
    stop("population_priors must be a bp_population_priors object");
  if (!proposalvars.inherits("bp_proposalvariance"))
    stop("proposalvars must be a bp_proposalvariance object");
  if (!population_startingvals.inherits("bp_population_startingvals"))
    stop("population_startingvals must be a bp_population_startingvals object");
  if (location_prior.size() != 1)
    stop("location_prior must be length 1");
  if (subject_data_list.size() != subject_startingvals_list.size())
    stop("Number of subjects in data and starting values must match");

  // Create shorter names for cleaner code
  int adj_iter       = pv_adjust_iter;
  int adj_max        = pv_adjust_max_iter;
  double biv_target  = bivariate_pv_target_ratio;
  double univ_target = univariate_pv_target_ratio;
  std::string loc_prior = Rcpp::as<std::string>(location_prior);

  // Random-effects distribution: Student-t (default, via per-pulse t-scale
  // kappa) or Gaussian (kappa fixed at 1). Optional element of the priors list
  // so the exported signature is unchanged; absent -> Student-t. Read as a real
  // so a logical or numeric both work.
  bool student_t_pulses = true;
  if (population_priors.containsElementNamed("student_t_pulses")) {
    student_t_pulses = (Rf_asReal(population_priors["student_t_pulses"]) != 0.0);
  }

  // Pulse random-effects scale: log-normal (papers/Horton) vs natural-scale
  // truncated-normal (research). And SD prior: uniform vs half-Cauchy. Both are
  // optional elements of the priors list so the exported signature is unchanged;
  // absent -> false (natural-scale / half-Cauchy) to preserve legacy behavior.
  // Read as reals so a logical or numeric both work.
  bool lognormal_pulses = false;
  if (population_priors.containsElementNamed("lognormal_pulses")) {
    lognormal_pulses = (Rf_asReal(population_priors["lognormal_pulses"]) != 0.0);
  }
  bool uniform_sd_prior = false;
  if (population_priors.containsElementNamed("uniform_sd_prior")) {
    uniform_sd_prior = (Rf_asReal(population_priors["uniform_sd_prior"]) != 0.0);
  }

  Rcpp::Rcout << "==================================================" << std::endl;
  Rcpp::Rcout << "Population Model MCMC" << std::endl;
  Rcpp::Rcout << "==================================================" << std::endl;
  Rcpp::Rcout << "Number of subjects: " << subject_data_list.size() << std::endl;
  Rcpp::Rcout << "Location prior: " << loc_prior << std::endl;
  Rcpp::Rcout << "MCMC iterations: " << mcmc_iterations << std::endl;
  Rcpp::Rcout << "Thin: " << thin << std::endl;
  Rcpp::Rcout << "Burnin: " << burnin << std::endl;
  Rcpp::Rcout << "==================================================" << std::endl;

  //
  // Create PopulationPriors object
  //
  PopulationPriors pop_priors(
    population_priors["mass_mean_mean"],
    population_priors["mass_mean_var"],
    population_priors["width_mean_mean"],
    population_priors["width_mean_var"],
    population_priors["baseline_mean_mean"],
    population_priors["baseline_mean_var"],
    population_priors["halflife_mean_mean"],
    population_priors["halflife_mean_var"],
    population_priors["mass_mean_sd_max"],
    population_priors["width_mean_sd_max"],
    population_priors["baseline_sd_max"],
    population_priors["halflife_sd_max"],
    population_priors["mass_sd_max"],
    population_priors["width_sd_max"],
    population_priors["error_alpha"],
    population_priors["error_beta"],
    population_priors["pulse_count"],
    population_priors["strauss_repulsion"],
    population_priors["strauss_repulsion_range"]
  );

  // Store the pulse parameterization flags on the priors so the population-level
  // SD draw (Pop_DrawPulseSDs) can reach them via the `population` pointer, and so
  // pop_mcmc_iteration can propagate them onto each subject every iteration.
  pop_priors.lognormal_pulses = lognormal_pulses;
  pop_priors.uniform_sd_prior = uniform_sd_prior;

  //
  // Create PopulationEstimates object (with starting values)
  //
  PopulationEstimates pop_estimates(
    population_startingvals["mass_mean"],
    population_startingvals["width_mean"],
    population_startingvals["baseline_mean"],
    population_startingvals["halflife_mean"],
    population_startingvals["mass_mean_sd"],
    population_startingvals["width_mean_sd"],
    population_startingvals["baseline_sd"],
    population_startingvals["halflife_sd"],
    population_startingvals["mass_sd"],
    population_startingvals["width_sd"],
    population_startingvals["errorsq"]
  );

  //
  // Create Patient objects for each subject
  //
  std::vector<Patient> subjects;
  int n_subjects = subject_data_list.size();

  for (int s = 0; s < n_subjects; s++) {

    // Extract data for this subject
    List subj_data = subject_data_list[s];

    // Validate subject data structure
    if (!subj_data.containsElementNamed("time") ||
        !subj_data.containsElementNamed("concentration")) {
      stop("Subject %d missing 'time' or 'concentration' field", s + 1);
    }

    NumericVector time = subj_data["time"];
    NumericVector concentration = subj_data["concentration"];

    // Validate vector dimensions
    if (time.size() != concentration.size()) {
      stop("Subject %d: time and concentration must have same length", s + 1);
    }
    if (time.size() == 0) {
      stop("Subject %d has no observations", s + 1);
    }

    PatientData data(time, concentration);

    // Extract starting values for this subject
    List subj_sv = subject_startingvals_list[s];
    PatientEstimates estimates(
      subj_sv["baseline"],
      subj_sv["halflife"],
      pop_estimates.errorsq,  // Error is shared across subjects
      subj_sv["mass_mean"],
      subj_sv["width_mean"]
    );  // Note: No mass_sd/width_sd - those come from population

    // Create PatientPriors from current population estimates
    // (will be updated each iteration)
    PatientPriors priors(
      pop_estimates.baseline_mean,
      pop_estimates.get_baseline_variance(),
      pop_estimates.halflife_mean,
      pop_estimates.get_halflife_variance(),
      pop_estimates.mass_mean,
      pop_estimates.get_mass_mean_variance(),
      pop_estimates.width_mean,
      pop_estimates.get_width_mean_variance(),
      pop_estimates.mass_mean_sd,
      pop_estimates.width_mean_sd
    );

    // Create Patient object and add to vector
    Patient patient(data, priors, estimates);
    patient.gaussian_random_effects = !student_t_pulses;
    patient.lognormal_pulses        = lognormal_pulses;
    patient.uniform_sd_prior        = uniform_sd_prior;
    subjects.push_back(patient);
  }

  //
  // Create Population object
  //
  Population population(subjects, pop_estimates, pop_priors);

  Rcpp::Rcout << "Population initialized with " << population.num_subjects
              << " subjects and " << population.get_total_observations()
              << " total observations" << std::endl;

  //
  // Create sampler objects
  //
  PopulationSamplers samplers(proposalvars, adj_iter, adj_max, biv_target,
                               univ_target, verbose, verbose_iter, loc_prior,
                               lognormal_pulses);

  //
  // Create output objects (chains)
  //
  PopulationChains chains(mcmc_iterations, thin, burnin, n_subjects,
                           verbose, verbose_iter);

  //
  // Run MCMC iterations
  //
  Rcpp::Rcout << "Starting MCMC..." << std::endl;

  for (int iteration = 0; iteration < mcmc_iterations; iteration++) {

    checkUserInterrupt();

    // Print diagnostic output periodically
    chains.print_diagnostic_output(&population, iteration);

    // Execute one MCMC iteration
    population_mcmc_iteration(&population, samplers, iteration);

    // Save to chains
    chains.save_sample(&population, iteration);
  }

  Rcpp::Rcout << "MCMC complete!" << std::endl;

  //
  // Return results
  //
  return chains.output(&population);
}
