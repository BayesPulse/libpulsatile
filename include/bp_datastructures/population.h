#ifndef GUARD_population_h
#define GUARD_population_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/patient.h>

//
// population.h
//   defining the population data structures for hierarchical Bayesian modeling
//
// Author: Matt Mulvahill
// Created: 12/30/24
// Notes:
//   These structures support population-level modeling where:
//   - Multiple subjects share common population-level parameters
//   - Subject-level parameters are drawn from population distributions
//   - Pulse-level variation is shared across all subjects
//

using namespace Rcpp;


//
// PopulationEstimates structure
//
// Holds population-level parameters that are updated via MCMC.
// These represent the hyperparameters that govern the distribution
// of subject-level parameters across the population.
//
struct PopulationEstimates {

  //
  // Population means (location parameters)
  //
  double mass_mean;        // μ_α: mean of subject mean pulse masses
  double width_mean;       // μ_ω: mean of subject mean pulse widths
  double baseline_mean;    // θ_b: population mean baseline
  double halflife_mean;    // θ_h: population mean half-life

  //
  // Subject-to-subject variation (between-subject SDs)
  //
  double mass_mean_sd;     // υ_α: SD of subject mean masses
  double width_mean_sd;    // υ_ω: SD of subject mean widths
  double baseline_sd;      // σ_b: SD of subject baselines
  double halflife_sd;      // σ_h: SD of subject half-lives

  //
  // Pulse-to-pulse variation (within-subject, shared across population)
  //
  double mass_sd;          // σ_α: pulse-to-pulse SD of mass
  double width_sd;         // σ_ω: pulse-to-pulse SD of width

  //
  // Model error
  //
  double errorsq;          // σ²_e: error variance

  // Default constructor
  PopulationEstimates() { }

  // Full constructor with starting values
  PopulationEstimates(double sv_mass_mean,
                      double sv_width_mean,
                      double sv_baseline_mean,
                      double sv_halflife_mean,
                      double sv_mass_mean_sd,
                      double sv_width_mean_sd,
                      double sv_baseline_sd,
                      double sv_halflife_sd,
                      double sv_mass_sd,
                      double sv_width_sd,
                      double sv_errorsq) {
    mass_mean      = sv_mass_mean;
    width_mean     = sv_width_mean;
    baseline_mean  = sv_baseline_mean;
    halflife_mean  = sv_halflife_mean;
    mass_mean_sd   = sv_mass_mean_sd;
    width_mean_sd  = sv_width_mean_sd;
    baseline_sd    = sv_baseline_sd;
    halflife_sd    = sv_halflife_sd;
    mass_sd        = sv_mass_sd;
    width_sd       = sv_width_sd;
    errorsq        = sv_errorsq;
  }

  // Utility functions
  double get_decay() const { return log(2) / halflife_mean; }
  double get_logerrorsq() const { return log(errorsq); }

  // Get variance instead of SD (for some calculations)
  double get_mass_mean_variance() const { return mass_mean_sd * mass_mean_sd; }
  double get_width_mean_variance() const { return width_mean_sd * width_mean_sd; }
  double get_baseline_variance() const { return baseline_sd * baseline_sd; }
  double get_halflife_variance() const { return halflife_sd * halflife_sd; }
  double get_mass_variance() const { return mass_sd * mass_sd; }
  double get_width_variance() const { return width_sd * width_sd; }
};


//
// PopulationPriors structure
//
// Holds user-specified prior distributions for population-level parameters.
// These are fixed throughout the MCMC and are not updated.
//
struct PopulationPriors {

  //
  // Priors on population means (normal priors)
  //
  double mass_mean_mean;
  double mass_mean_var;
  double width_mean_mean;
  double width_mean_var;
  double baseline_mean_mean;
  double baseline_mean_var;
  double halflife_mean_mean;
  double halflife_mean_var;

  //
  // Priors on subject-to-subject SDs (uniform priors, specify upper bounds)
  //
  double mass_mean_sd_max;
  double width_mean_sd_max;
  double baseline_sd_max;
  double halflife_sd_max;

  //
  // Priors on pulse-to-pulse SDs (uniform priors, specify upper bounds)
  //
  double mass_sd_max;
  double width_sd_max;

  //
  // Error prior (inverse-gamma)
  //
  double error_alpha;
  double error_beta;

  //
  // Pulse count prior (Poisson rate parameter)
  //
  double pulse_count_prior;

  //
  // Birth-death process priors (Strauss process)
  //
  double strauss_repulsion;       // gamma for interaction
  double strauss_repulsion_range; // range of interaction

  //
  // Pulse random-effects parameterization flags (population level)
  //
  // lognormal_pulses: when true (the published Horton parameterization) the pulse
  //   mass/width random effects are modeled on the LOG scale, log(theta) ~
  //   N(mu, sigma^2/kappa). The population-level pulse-to-pulse SD draw
  //   (Pop_DrawPulseSDs) then forms residuals as (log(theta) - subject_log_mean)
  //   and drops the truncated-normal normalizing constant. The subject means
  //   (mass_mean/width_mean) are means of the LOG random effects and are drawn by
  //   the per-subject SS_DrawFixedEffects with its own log-normal branch (threaded
  //   through PopulationSamplers). When false (natural-scale research option) the
  //   pulse SD draw keeps the natural-scale residuals and the truncation constant.
  // uniform_sd_prior: carried for parity with the single-subject/joint models. The
  //   population SD draws (Pop_DrawPulseSDs/Pop_DrawMeanSDs) ALREADY use a
  //   Uniform(0, max) prior unconditionally (matching Horton), so this flag does
  //   not alter population sampling; it is threaded onto each subject Patient so
  //   any per-subject SD machinery stays consistent.
  // Both default to false (C++ standalone default); the R specs flip them.
  bool lognormal_pulses = false;
  bool uniform_sd_prior = false;

  // Default constructor
  PopulationPriors() { }

  // Full constructor
  PopulationPriors(double prior_mass_mean_mean,
                   double prior_mass_mean_var,
                   double prior_width_mean_mean,
                   double prior_width_mean_var,
                   double prior_baseline_mean_mean,
                   double prior_baseline_mean_var,
                   double prior_halflife_mean_mean,
                   double prior_halflife_mean_var,
                   double prior_mass_mean_sd_max,
                   double prior_width_mean_sd_max,
                   double prior_baseline_sd_max,
                   double prior_halflife_sd_max,
                   double prior_mass_sd_max,
                   double prior_width_sd_max,
                   double prior_error_alpha,
                   double prior_error_beta,
                   double prior_pulse_count,
                   double prior_strauss_repulsion,
                   double prior_strauss_repulsion_range) {

    mass_mean_mean          = prior_mass_mean_mean;
    mass_mean_var           = prior_mass_mean_var;
    width_mean_mean         = prior_width_mean_mean;
    width_mean_var          = prior_width_mean_var;
    baseline_mean_mean      = prior_baseline_mean_mean;
    baseline_mean_var       = prior_baseline_mean_var;
    halflife_mean_mean      = prior_halflife_mean_mean;
    halflife_mean_var       = prior_halflife_mean_var;

    mass_mean_sd_max        = prior_mass_mean_sd_max;
    width_mean_sd_max       = prior_width_mean_sd_max;
    baseline_sd_max         = prior_baseline_sd_max;
    halflife_sd_max         = prior_halflife_sd_max;

    mass_sd_max             = prior_mass_sd_max;
    width_sd_max            = prior_width_sd_max;

    error_alpha             = prior_error_alpha;
    error_beta              = 1 / prior_error_beta; // Note: inverse as in original code
    pulse_count_prior       = prior_pulse_count;

    strauss_repulsion       = prior_strauss_repulsion;
    strauss_repulsion_range = prior_strauss_repulsion_range;
  }
};


//
// Population structure
//
// Main container for population-level modeling. Holds multiple subjects
// along with population-level parameters and priors.
//
struct Population {

  //
  // Data and parameters
  //
  std::vector<Patient> subjects;
  PopulationEstimates estimates;
  PopulationPriors priors;
  int num_subjects;

  // Default constructor
  Population() : num_subjects(0) { }

  // Constructor from vector of patients
  Population(std::vector<Patient> subj_list,
             PopulationEstimates est,
             PopulationPriors pri)
    : subjects(subj_list), estimates(est), priors(pri) {
    num_subjects = subjects.size();
  }

  //
  // Utility methods
  //

  // Get total number of observations across all subjects
  int get_total_observations() const {
    int total = 0;
    for (const auto& subj : subjects) {
      total += subj.data.time.n_elem;
    }
    return total;
  }

  // Get total number of pulses across all subjects
  int get_total_pulses() const {
    int total = 0;
    for (const auto& subj : subjects) {
      total += subj.get_pulsecount();
    }
    return total;
  }

  // Get all subject baselines as a vector
  arma::vec get_all_subject_baselines() const {
    arma::vec baselines(num_subjects);
    for (int i = 0; i < num_subjects; i++) {
      baselines(i) = subjects[i].estimates.baseline_halflife(0);
    }
    return baselines;
  }

  // Get all subject half-lives as a vector
  arma::vec get_all_subject_halflives() const {
    arma::vec halflives(num_subjects);
    for (int i = 0; i < num_subjects; i++) {
      halflives(i) = subjects[i].estimates.baseline_halflife(1);
    }
    return halflives;
  }

  // Get all subject mean masses as a vector
  arma::vec get_all_subject_mass_means() const {
    arma::vec mass_means(num_subjects);
    for (int i = 0; i < num_subjects; i++) {
      mass_means(i) = subjects[i].estimates.mass_mean;
    }
    return mass_means;
  }

  // Get all subject mean widths as a vector
  arma::vec get_all_subject_width_means() const {
    arma::vec width_means(num_subjects);
    for (int i = 0; i < num_subjects; i++) {
      width_means(i) = subjects[i].estimates.width_mean;
    }
    return width_means;
  }

  // Get all pulse masses across all subjects
  arma::vec get_all_pulse_masses() const {
    std::vector<double> masses;
    for (const auto& subj : subjects) {
      for (auto it = subj.pulses.begin(); it != subj.pulses.end(); ++it) {
        masses.push_back(it->mass);
      }
    }
    return arma::vec(masses);
  }

  // Get all pulse widths across all subjects
  arma::vec get_all_pulse_widths() const {
    std::vector<double> widths;
    for (const auto& subj : subjects) {
      for (auto it = subj.pulses.begin(); it != subj.pulses.end(); ++it) {
        widths.push_back(it->width);
      }
    }
    return arma::vec(widths);
  }

};


#endif
