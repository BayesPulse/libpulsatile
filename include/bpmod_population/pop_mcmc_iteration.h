#ifndef GUARD_bpmod_population_mcmc_iteration_h
#define GUARD_bpmod_population_mcmc_iteration_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/population.h>
#include <bp_datastructures/patient.h>
#include <bpmod_singlesubject/bpmod_singlesubject.h>
#include <bpmod_population/bpmod_population.h>

//
// pop_mcmc_iteration.h
//   Orchestrates a single MCMC iteration for the population model
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// The population MCMC iteration proceeds in two stages:
//
// Stage 1: Subject-level updates (for each subject s = 1, ..., n_subjects)
//   a. Birth-death process for pulse count
//   b. Subject-level parameters: baseline_s, halflife_s, mass_mean_s, width_mean_s
//   c. Pulse-level parameters: For each pulse k in subject s:
//      - mass_k,s, width_k,s, location_k,s
//      - tvarscale_mass_k,s, tvarscale_width_k,s
//
// Stage 2: Population-level updates
//   a. Population means: μ_α, μ_ω, θ_b, θ_h (Gibbs)
//   b. Subject-to-subject SDs: υ_α, υ_ω, σ_b, σ_h (MH)
//   c. Pulse-to-pulse SDs: σ_α, σ_ω (MH)
//   d. Error variance: σ²_e (Gibbs)
//

//
// PopulationSamplers struct
//   Container for all MCMC samplers needed for population model
//
struct PopulationSamplers {

  // Birth-death process (one per subject, but can reuse same object)
  BirthDeathProcess birth_death;

  // Subject-level samplers (reuse single-subject samplers)
  SS_DrawBaselineHalflife *draw_blhl;
  SS_DrawFixedEffects *draw_fixeff_mass;
  SS_DrawFixedEffects *draw_fixeff_width;

  // Pulse-level samplers (reuse single-subject samplers)
  SS_DrawLocations *draw_locations;
  SS_DrawRandomEffects *draw_masses;
  SS_DrawRandomEffects *draw_widths;
  SS_DrawTVarScale *draw_tvarscale_mass;
  SS_DrawTVarScale *draw_tvarscale_width;

  // Population-level samplers (Gibbs)
  Pop_DrawMeans draw_pop_means;
  Pop_DrawBaselineHalflifeMeans draw_pop_bh_means;
  Pop_DrawError draw_pop_error;

  // Population-level samplers (MH)
  Pop_DrawMeanSDs *draw_pop_mass_mean_sd;
  Pop_DrawMeanSDs *draw_pop_width_mean_sd;
  Pop_DrawBaselineHalflifeSD *draw_pop_baseline_sd;
  Pop_DrawBaselineHalflifeSD *draw_pop_halflife_sd;
  Pop_DrawPulseSDs *draw_pop_mass_sd;
  Pop_DrawPulseSDs *draw_pop_width_sd;

  // Constructor
  PopulationSamplers(Rcpp::List proposalvars,
                     int adj_iter,
                     int adj_max,
                     double biv_target,
                     double univ_target,
                     bool verbose,
                     int verbose_iter,
                     std::string loc_prior,
                     bool lognormal = false) {

    // Subject-level samplers
    arma::vec bhl_pv = { Rcpp::as<double>(proposalvars["baseline"]),
                         Rcpp::as<double>(proposalvars["halflife"]) };
    draw_blhl = new SS_DrawBaselineHalflife(bhl_pv, adj_iter, adj_max,
                                            biv_target, verbose, verbose_iter);

    // Thread the population's log-normal flag into the per-subject mean draws so
    // the subject mass/width MEAN full conditionals use log-scale residuals under
    // the papers' parameterization (the MMH container type here is `bool`, so the
    // flag cannot be read from the Patient and must be passed at construction).
    draw_fixeff_mass = new SS_DrawFixedEffects(
      Rcpp::as<double>(proposalvars["mass_mean"]),
      adj_iter, adj_max, univ_target, false, verbose, verbose_iter, lognormal);

    draw_fixeff_width = new SS_DrawFixedEffects(
      Rcpp::as<double>(proposalvars["width_mean"]),
      adj_iter, adj_max, univ_target, true, verbose, verbose_iter, lognormal);

    // Pulse-level samplers
    if (loc_prior == "strauss") {
      draw_locations = new SS_DrawLocationsStrauss(
        Rcpp::as<double>(proposalvars["location"]),
        adj_iter, adj_max, univ_target, verbose, verbose_iter);
    } else {
      Rcpp::stop("Order statistic prior is not yet implemented");
    }

    draw_masses = new SS_DrawRandomEffects(
      Rcpp::as<double>(proposalvars["pulse_mass"]),
      adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

    draw_widths = new SS_DrawRandomEffects(
      Rcpp::as<double>(proposalvars["pulse_width"]),
      adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

    draw_tvarscale_mass = new SS_DrawTVarScale(
      Rcpp::as<double>(proposalvars["sdscale_pulse_mass"]),
      adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

    draw_tvarscale_width = new SS_DrawTVarScale(
      Rcpp::as<double>(proposalvars["sdscale_pulse_width"]),
      adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

    // Population-level MH samplers
    draw_pop_mass_mean_sd = new Pop_DrawMeanSDs(
      Rcpp::as<double>(proposalvars["pop_mass_mean_sd"]),
      adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

    draw_pop_width_mean_sd = new Pop_DrawMeanSDs(
      Rcpp::as<double>(proposalvars["pop_width_mean_sd"]),
      adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

    draw_pop_baseline_sd = new Pop_DrawBaselineHalflifeSD(
      Rcpp::as<double>(proposalvars["pop_baseline_sd"]),
      adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

    draw_pop_halflife_sd = new Pop_DrawBaselineHalflifeSD(
      Rcpp::as<double>(proposalvars["pop_halflife_sd"]),
      adj_iter, adj_max, univ_target, true, verbose, verbose_iter);

    draw_pop_mass_sd = new Pop_DrawPulseSDs(
      Rcpp::as<double>(proposalvars["pop_mass_sd"]),
      adj_iter, adj_max, univ_target, false, verbose, verbose_iter);

    draw_pop_width_sd = new Pop_DrawPulseSDs(
      Rcpp::as<double>(proposalvars["pop_width_sd"]),
      adj_iter, adj_max, univ_target, true, verbose, verbose_iter);
  }

  // Destructor to clean up dynamically allocated samplers
  ~PopulationSamplers() {
    delete draw_blhl;
    delete draw_fixeff_mass;
    delete draw_fixeff_width;
    delete draw_locations;
    delete draw_masses;
    delete draw_widths;
    delete draw_tvarscale_mass;
    delete draw_tvarscale_width;
    delete draw_pop_mass_mean_sd;
    delete draw_pop_width_mean_sd;
    delete draw_pop_baseline_sd;
    delete draw_pop_halflife_sd;
    delete draw_pop_mass_sd;
    delete draw_pop_width_sd;
  }
};


//
// population_mcmc_iteration()
//   Execute one complete MCMC iteration for the population model
//
inline void population_mcmc_iteration(Population *population,
                                       PopulationSamplers &samplers,
                                       int iteration) {

  //
  // STAGE 1: Subject-level updates
  //
  for (auto& subject : population->subjects) {

    // Update subject priors from current population estimates
    // This is critical - subject priors come from population!
    subject.priors.baseline_mean     = population->estimates.baseline_mean;
    subject.priors.baseline_variance = population->estimates.get_baseline_variance();
    subject.priors.halflife_mean     = population->estimates.halflife_mean;
    subject.priors.halflife_variance = population->estimates.get_halflife_variance();
    subject.priors.mass_mean         = population->estimates.mass_mean;
    subject.priors.mass_variance     = population->estimates.get_mass_mean_variance();
    subject.priors.width_mean        = population->estimates.width_mean;
    subject.priors.width_variance    = population->estimates.get_width_mean_variance();

    // Note: mass_sd and width_sd used in pulse samplers come from population
    // We'll need to temporarily set them on the subject for compatibility
    subject.estimates.mass_sd  = population->estimates.mass_sd;
    subject.estimates.width_sd = population->estimates.width_sd;

    // Propagate the pulse random-effects parameterization flags from the
    // population onto each subject every iteration. SS_DrawRandomEffects and
    // SS_DrawTVarScale read subject.lognormal_pulses DIRECTLY, and the birth-death
    // process reads it to initialize new pulses, so the flag must be live on the
    // subject on every iteration (belt-and-suspenders against any per-iteration
    // subject rebuild). uniform_sd_prior is carried for parity; the population SD
    // draws are Uniform-only regardless.
    subject.lognormal_pulses = population->priors.lognormal_pulses;
    subject.uniform_sd_prior = population->priors.uniform_sd_prior;

    // Propagate the pulse-count (birth-rate) and Strauss priors from the
    // population to the subject. The per-subject PatientPriors are built with
    // the population constructor, which hardcodes these to 0; without this the
    // birth-death total_birth_rate is identically 0, no pulse can be born, and
    // every subject collapses to a single pulse - biasing baseline/mass high
    // and halflife low. (PopulationPriors carries these as pulse_count_prior.)
    subject.priors.pulse_count             = population->priors.pulse_count_prior;
    subject.priors.strauss_repulsion       = population->priors.strauss_repulsion;
    subject.priors.strauss_repulsion_range = population->priors.strauss_repulsion_range;

    // Birth-death process for pulse count
    samplers.birth_death.sample(&subject, false, iteration);

    // Subject-level parameters
    samplers.draw_blhl->sample(&subject, &subject.estimates.baseline_halflife, iteration);
    samplers.draw_fixeff_mass->sample(&subject, &subject.estimates.mass_mean, iteration);
    samplers.draw_fixeff_width->sample(&subject, &subject.estimates.width_mean, iteration);

    // Pulse-level parameters (for all pulses in this subject)
    samplers.draw_locations->sample_pulses(&subject, iteration);
    samplers.draw_masses->sample_pulses(&subject, iteration);
    samplers.draw_widths->sample_pulses(&subject, iteration);
    // Skip the t-scale (kappa) draws under Gaussian random effects; kappa stays
    // fixed at 1 for every pulse.
    if (!subject.gaussian_random_effects) {
      samplers.draw_tvarscale_mass->sample_pulses(&subject, iteration);
      samplers.draw_tvarscale_width->sample_pulses(&subject, iteration);
    }
  }

  //
  // STAGE 2: Population-level updates
  //

  // Population means (Gibbs samplers)
  samplers.draw_pop_means.sample(population);
  samplers.draw_pop_bh_means.sample(population);

  // Subject-to-subject SDs (MH samplers)
  // NOTE: the population is passed as BOTH sampling_unit and container. The
  // MH parameter_support() reads the uniform-prior upper bound (*_sd_max) from
  // the container; the 3-arg sample() overload would default-construct an empty
  // Population whose PopulationPriors leaves those bounds uninitialized, causing
  // every proposal to be rejected and the SD chains to freeze at their starting
  // values. Passing the real population as the container fixes this.
  samplers.draw_pop_mass_mean_sd->sample(population,
                                         &population->estimates.mass_mean_sd,
                                         population,
                                         iteration);
  samplers.draw_pop_width_mean_sd->sample(population,
                                          &population->estimates.width_mean_sd,
                                          population,
                                          iteration);
  samplers.draw_pop_baseline_sd->sample(population,
                                        &population->estimates.baseline_sd,
                                        population,
                                        iteration);
  samplers.draw_pop_halflife_sd->sample(population,
                                        &population->estimates.halflife_sd,
                                        population,
                                        iteration);

  // Pulse-to-pulse SDs (MH samplers)
  samplers.draw_pop_mass_sd->sample(population,
                                    &population->estimates.mass_sd,
                                    population,
                                    iteration);
  samplers.draw_pop_width_sd->sample(population,
                                     &population->estimates.width_sd,
                                     population,
                                     iteration);

  // Error variance (Gibbs sampler)
  samplers.draw_pop_error.sample(population);
}


#endif
