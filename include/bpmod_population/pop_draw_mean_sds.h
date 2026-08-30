#ifndef GUARD_bpmod_population_draw_mean_sds_h
#define GUARD_bpmod_population_draw_mean_sds_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/population.h>

//
// Pop_DrawMeanSDs
//   Metropolis-Hastings sampler for subject-to-subject SDs of mean parameters
//   
//   Updates:
//     - population.estimates.mass_mean_sd (υ_α)
//     - population.estimates.width_mean_sd (υ_ω)
//
//   Prior: υ ~ Uniform(0, υ_max)
//
//   Likelihood: θ_i ~ N(μ, υ²) where θ_i are subject-level means
//
//   Log-normal note: this draw is UNCHANGED under the log-normal pulse
//   parameterization. The subject-level means θ_i (mass_mean/width_mean) are the
//   same quantities SS_DrawFixedEffects targets; under log-normal they are means
//   of the LOG random effects and may be negative. The likelihood here is an
//   UNTRUNCATED normal, N(θ_i | μ, υ²), with NO pnorm truncation constant (see
//   posterior_function below), so it is already exactly the log-normal target --
//   residuals stay (θ_i - μ) with no re-logging. The natural-scale ("research")
//   behavior is likewise unchanged: it omits truncation too, which is the
//   existing (documented) research approximation and is deliberately preserved.
//

class Pop_DrawMeanSDs : 
  public ModifiedMetropolisHastings<Population, Population, double, ProposalVariance>
{

  public:

    // Constructor
    Pop_DrawMeanSDs(double in_pv,
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
          get_sd_        = &PopulationEstimates::width_mean_sd;
          get_mean_      = &PopulationEstimates::width_mean;
          get_max_       = &PopulationPriors::width_mean_sd_max;
          get_subj_vals_ = &Population::get_all_subject_width_means;
          parameter_name = "Subject SD of mean widths (υ_ω)";
        } else {
          get_sd_        = &PopulationEstimates::mass_mean_sd;
          get_mean_      = &PopulationEstimates::mass_mean;
          get_max_       = &PopulationPriors::mass_mean_sd_max;
          get_subj_vals_ = &Population::get_all_subject_mass_means;
          parameter_name = "Subject SD of mean masses (υ_α)";
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
// Defined functions for Pop_DrawMeanSDs
//------------------------------------------------------------

// parameter_support()
//   Check if proposal is within uniform prior support [0, max]
bool Pop_DrawMeanSDs::parameter_support(double val, Population *population) {
  double max_val = (population->priors).*get_max_;
  return (val > 0.0 && val < max_val);
}


// posterior_function()
//   Calculate acceptance ratio for MH sampler
//   
//   Log posterior ratio:
//     log(π(υ* | data) / π(υ | data))
//     = log(L(data | υ*) / L(data | υ)) + log(π(υ*) / π(υ))
//     = log(L(data | υ*) / L(data | υ))  (uniform prior cancels)
//     
//   where L(data | υ) = Π_i N(θ_i | μ, υ²)
//
double Pop_DrawMeanSDs::posterior_function(Population *population, 
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
  // log(L(new) / L(old)) = n * log(old_sd / new_sd) + 
  //                        0.5 * ssq * (1/old²  - 1/new²)
  double log_ratio = n_subjects * (log(current_sd) - log(proposal)) + 
                     0.5 * ssq * (1.0 / (current_sd * current_sd) - 
                                  1.0 / (proposal * proposal));
  
  // Prior ratio is 1 for uniform prior (cancels in log space)
  
  return log_ratio;
}

#endif
