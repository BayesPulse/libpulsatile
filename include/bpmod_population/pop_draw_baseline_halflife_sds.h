#ifndef GUARD_bpmod_population_draw_baseline_halflife_sds_h
#define GUARD_bpmod_population_draw_baseline_halflife_sds_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/population.h>

//
// Pop_DrawBaselineHalflifeSD's
//   Metropolis-Hastings sampler for subject-to-subject SDs of baseline/halflife
//   
//   Updates:
//     - population.estimates.baseline_sd (σ_b)
//     - population.estimates.halflife_sd (σ_h)
//
//   Prior: σ ~ Uniform(0, σ_max)
//   
//   Likelihood: θ_i ~ N(θ, σ²) where θ_i are subject-level baseline/halflife
//

class Pop_DrawBaselineHalflifeSD : 
  public ModifiedMetropolisHastings<Population, Population, double, ProposalVariance>
{

  public:

    // Constructor
    Pop_DrawBaselineHalflifeSD(double in_pv,
                               int in_adjust_iter,
                               int in_max_iter,
                               double in_target_ratio,
                               bool for_halflife,
                               bool verbose,
                               int verbose_iter) :
      ModifiedMetropolisHastings <Population, Population, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {

        // Choose which parameter to sample: halflife or baseline
        if (for_halflife) {
          get_sd_        = &PopulationEstimates::halflife_sd;
          get_mean_      = &PopulationEstimates::halflife_mean;
          get_max_       = &PopulationPriors::halflife_sd_max;
          get_subj_vals_ = &Population::get_all_subject_halflives;
          parameter_name = "Subject SD of halflife (σ_h)";
        } else {
          get_sd_        = &PopulationEstimates::baseline_sd;
          get_mean_      = &PopulationEstimates::baseline_mean;
          get_max_       = &PopulationPriors::baseline_sd_max;
          get_subj_vals_ = &Population::get_all_subject_baselines;
          parameter_name = "Subject SD of baseline (σ_b)";
        }

      }

  private:

    // Member pointers for flexibility
    double PopulationEstimates::*get_sd_;
    double PopulationEstimates::*get_mean_;
    double PopulationPriors::*get_max_;
    arma::vec (Population::*get_subj_vals_)() const;

    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; }

    bool parameter_support(double val, Population *population);
    double posterior_function(Population *population, double proposal, Population *notused);

};


//------------------------------------------------------------
// Defined functions for Pop_DrawBaselineHalflifeSD
//------------------------------------------------------------

// parameter_support()
//   Check if proposal is within uniform prior support [0, max]
bool Pop_DrawBaselineHalflifeSD::parameter_support(double val, Population *population) {
  double max_val = (population->priors).*get_max_;
  return (val > 0.0 && val < max_val);
}


// posterior_function()
//   Calculate acceptance ratio for MH sampler
//   Same logic as Pop_DrawMeanSDs but for baseline/halflife
//
double Pop_DrawBaselineHalflifeSD::posterior_function(Population *population, 
                                                      double proposal, 
                                                      Population *notused) {

  double current_sd = (population->estimates).*get_sd_;
  double pop_mean   = (population->estimates).*get_mean_;
  int n_subjects    = population->num_subjects;
  
  // Get subject-level values
  arma::vec subj_vals = (population->*get_subj_vals_)();
  
  // Calculate sum of squared deviations from population mean
  double ssq = 0.0;
  for (int i = 0; i < n_subjects; i++) {
    double dev = subj_vals(i) - pop_mean;
    ssq += dev * dev;
  }
  
  // Log likelihood ratio
  double log_ratio = n_subjects * (log(current_sd) - log(proposal)) + 
                     0.5 * ssq * (1.0 / (current_sd * current_sd) - 
                                  1.0 / (proposal * proposal));
  
  return log_ratio;
}

#endif
