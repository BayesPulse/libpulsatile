#ifndef GUARD_bpmod_population_draw_pulse_sds_h
#define GUARD_bpmod_population_draw_pulse_sds_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/population.h>

//
// Pop_DrawPulseSDs
//   Metropolis-Hastings sampler for pulse-to-pulse SDs (shared across subjects)
//   
//   Updates:
//     - population.estimates.mass_sd (σ_α)
//     - population.estimates.width_sd (σ_ω)
//
//   Prior: σ ~ Uniform(0, σ_max)
//   
//   Likelihood: For all subjects s and pulses k:
//               α_k,s ~ TruncatedT(μ_α,s, σ_α / sqrt(κ_α,k,s))
//               where truncation is at 0
//
//   Note: This is shared across ALL subjects - affects every pulse in the population
//

class Pop_DrawPulseSDs : 
  public ModifiedMetropolisHastings<Population, Population, double, ProposalVariance>
{

  public:

    // Constructor
    Pop_DrawPulseSDs(double in_pv,
                     int in_adjust_iter,
                     int in_max_iter,
                     double in_target_ratio,
                     bool for_width,
                     bool verbose,
                     int verbose_iter) :
      ModifiedMetropolisHastings <Population, Population, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {

        // Choose which parameter to sample: width or mass
        if (for_width) {
          get_sd_          = &PopulationEstimates::width_sd;
          get_max_         = &PopulationPriors::width_sd_max;
          get_subj_mean_   = &PatientEstimates::width_mean;
          get_pulse_value_ = &PulseEstimates::width;
          get_tvarscale_   = &PulseEstimates::tvarscale_width;
          parameter_name   = "Pulse-to-pulse SD of width (σ_ω)";
        } else {
          get_sd_          = &PopulationEstimates::mass_sd;
          get_max_         = &PopulationPriors::mass_sd_max;
          get_subj_mean_   = &PatientEstimates::mass_mean;
          get_pulse_value_ = &PulseEstimates::mass;
          get_tvarscale_   = &PulseEstimates::tvarscale_mass;
          parameter_name   = "Pulse-to-pulse SD of mass (σ_α)";
        }

      }

  private:

    // Member pointers for flexibility
    double PopulationEstimates::*get_sd_;
    double PopulationPriors::*get_max_;
    double PatientEstimates::*get_subj_mean_;
    double PulseEstimates::*get_pulse_value_;
    double PulseEstimates::*get_tvarscale_;

    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; }

    bool parameter_support(double val, Population *population);

  public:

    // Exposed publicly so the closed-form log-ratio can be unit tested; the MH
    // base dispatches through its own private virtual, so this does not change
    // normal sampling behavior.
    double posterior_function(Population *population, double proposal, Population *notused);

};


//------------------------------------------------------------
// Defined functions for Pop_DrawPulseSDs
//------------------------------------------------------------

// parameter_support()
//   Check if proposal is within uniform prior support [0, max]
bool Pop_DrawPulseSDs::parameter_support(double val, Population *population) {
  double max_val = (population->priors).*get_max_;
  return (val > 0.0 && val < max_val);
}


// posterior_function()
//   Calculate acceptance ratio for MH sampler
//   
//   For truncated t-distribution likelihood:
//     f(x | μ, σ, κ) ∝ (1/σ) * exp(-0.5 * κ * (x-μ)²/σ²) / Φ(μ*sqrt(κ)/σ)
//   
//   where Φ is the standard normal CDF (normalizing constant for truncation)
//
double Pop_DrawPulseSDs::posterior_function(Population *population, 
                                            double proposal, 
                                            Population *notused) {

  double current_sd = (population->estimates).*get_sd_;
  bool lognormal    = population->priors.lognormal_pulses;

  double log_ratio = 0.0;
  int total_pulses = 0;
  double ssq_weighted = 0.0;
  double old_norm_const = 0.0;
  double new_norm_const = 0.0;

  // Loop over all subjects and their pulses
  for (auto& subject : population->subjects) {

    double subj_mean = (subject.estimates).*get_subj_mean_;

    for (auto it = subject.pulses.begin(); it != subject.pulses.end(); ++it) {

      double pulse_val = (*it).*get_pulse_value_;
      double tvarscale = (*it).*get_tvarscale_;

      total_pulses++;

      if (lognormal) {
        // Log-normal parameterization (Horton): log(theta) ~ N(mu, sigma^2/kappa),
        // where mu = subj_mean is the mean of the LOG random effects (the same
        // quantity SS_DrawFixedEffects targets under log-normal). Residuals are
        // formed on the log scale, (log(theta) - subj_mean), weighted by the
        // per-pulse t-scale kappa. There is no truncation at 0, so the
        // truncated-normal normalizing constants stay 0 and drop out below.
        double logre = log(pulse_val);
        ssq_weighted += tvarscale * (logre - subj_mean) * (logre - subj_mean);
      } else {
        // Natural-scale truncated-normal parameterization (research option).
        // Weighted sum of squared deviations
        ssq_weighted += tvarscale * (pulse_val - subj_mean) * (pulse_val - subj_mean);

        // Normalizing constants for truncated distribution
        // stdx = μ / (σ / sqrt(κ)) = μ * sqrt(κ) / σ
        double stdx_old = subj_mean * sqrt(tvarscale) / current_sd;
        double stdx_new = subj_mean * sqrt(tvarscale) / proposal;

        // Log of normalizing constant (using pnorm in log space)
        old_norm_const += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);  // log = TRUE
        new_norm_const += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
      }
    }
  }
  
  // Log likelihood ratio components:
  // 1. From the (1/σ) terms: n * log(σ_old / σ_new)
  double first_part = total_pulses * (log(current_sd) - log(proposal));
  
  // 2. From the exp(-0.5 * κ * (x-μ)²/σ²) terms
  double second_part = 0.5 * ssq_weighted * 
                       (1.0 / (current_sd * current_sd) - 1.0 / (proposal * proposal));
  
  // 3. From the normalizing constants
  double third_part = old_norm_const - new_norm_const;
  
  log_ratio = first_part + second_part + third_part;
  
  // Prior ratio (uniform prior cancels)
  
  return log_ratio;
}

#endif
