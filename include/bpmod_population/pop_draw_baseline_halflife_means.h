#ifndef GUARD_bpmod_population_draw_baseline_halflife_means_h
#define GUARD_bpmod_population_draw_baseline_halflife_means_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_datastructures/population.h>

//
// Pop_DrawBaselineHalflifeMeans
//   Draw population mean parameters for baseline and half-life using Gibbs
//   
//   Updates:
//     - population.estimates.baseline_mean (θ_b)
//     - population.estimates.halflife_mean (θ_h)
//
//   For parameter θ (either θ_b or θ_h):
//     Prior: θ ~ N(θ_0, σ²_0)
//     Data:  θ_i ~ N(θ, σ²_θ) for i = 1, ..., n_subjects
//     Posterior: θ | θ_1, ..., θ_n ~ N(θ_n, σ²_n)
//   
//   Where:
//     σ²_n = 1 / (1/σ²_0 + n/σ²_θ)
//     θ_n  = σ²_n * (θ_0/σ²_0 + Σθ_i/σ²_θ)
//

class Pop_DrawBaselineHalflifeMeans 
{

  public:
    // Constructors
    Pop_DrawBaselineHalflifeMeans() { }

    //
    // Sample both baseline_mean and halflife_mean
    //
    void sample(Population *population) {
      sample_baseline_mean(population);
      sample_halflife_mean(population);
    }

    //
    // Sample baseline_mean (θ_b) only
    //
    void sample_baseline_mean(Population *population) {
      
      int n = population->num_subjects;
      
      // Prior parameters
      double prior_mean = population->priors.baseline_mean_mean;
      double prior_var  = population->priors.baseline_mean_var;
      
      // Subject-to-subject variance
      double subj_var = population->estimates.get_baseline_variance();
      
      // Get all subject baselines
      arma::vec subject_baselines = population->get_all_subject_baselines();
      double sum_baselines = arma::sum(subject_baselines);
      
      // Posterior parameters (conjugate normal)
      double post_prec = 1.0 / prior_var + n / subj_var;
      double post_var  = 1.0 / post_prec;
      double post_mean = post_var * (prior_mean / prior_var + sum_baselines / subj_var);
      
      // Draw from posterior
      population->estimates.baseline_mean = R::rnorm(post_mean, sqrt(post_var));
    }

    //
    // Sample halflife_mean (θ_h) only
    //
    void sample_halflife_mean(Population *population) {
      
      int n = population->num_subjects;
      
      // Prior parameters
      double prior_mean = population->priors.halflife_mean_mean;
      double prior_var  = population->priors.halflife_mean_var;
      
      // Subject-to-subject variance
      double subj_var = population->estimates.get_halflife_variance();
      
      // Get all subject half-lives
      arma::vec subject_halflives = population->get_all_subject_halflives();
      double sum_halflives = arma::sum(subject_halflives);
      
      // Posterior parameters (conjugate normal)
      double post_prec = 1.0 / prior_var + n / subj_var;
      double post_var  = 1.0 / post_prec;
      double post_mean = post_var * (prior_mean / prior_var + sum_halflives / subj_var);
      
      // Draw from posterior
      population->estimates.halflife_mean = R::rnorm(post_mean, sqrt(post_var));
    }

};

#endif
