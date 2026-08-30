#ifndef GUARD_bpmod_population_draw_means_h
#define GUARD_bpmod_population_draw_means_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_datastructures/population.h>

//
// Pop_DrawMeans
//   Draw population mean parameters for mass and width using Gibbs sampling
//   
//   Updates:
//     - population.estimates.mass_mean (μ_α)
//     - population.estimates.width_mean (μ_ω)
//
//   For parameter μ (either μ_α or μ_ω):
//     Prior: μ ~ N(μ_0, τ²_0)
//     Data:  θ_i ~ N(μ, υ²) for i = 1, ..., n_subjects
//     Posterior: μ | θ_1, ..., θ_n ~ N(μ_n, τ²_n)
//   
//   Where:
//     τ²_n = 1 / (1/τ²_0 + n/υ²)
//     μ_n  = τ²_n * (μ_0/τ²_0 + Σθ_i/υ²)
//
//   Log-normal note: this Gibbs draw is UNCHANGED under the log-normal pulse
//   parameterization. It is a normal-normal conjugate update that is
//   scale-agnostic: the subject means θ_i and the population mean μ live on
//   whatever scale the fixed-effects draw uses (log scale under log-normal,
//   natural scale otherwise), and the conjugate posterior N(μ_n, τ²_n) is exact
//   for an UNTRUNCATED normal -- which is exactly what the log-normal model
//   gives. The code makes no positivity/truncation assumption (μ is drawn with a
//   plain R::rnorm and may be negative, as required for a log-scale mean), so no
//   branch is needed.
//

class Pop_DrawMeans 
{

  public:
    // Constructors
    Pop_DrawMeans() { }

    //
    // Sample both mass_mean and width_mean
    //
    void sample(Population *population) {
      sample_mass_mean(population);
      sample_width_mean(population);
    }

    //
    // Sample mass_mean (μ_α) only
    //
    void sample_mass_mean(Population *population) {
      
      int n = population->num_subjects;
      
      // Prior parameters
      double prior_mean = population->priors.mass_mean_mean;
      double prior_var  = population->priors.mass_mean_var;
      
      // Subject-to-subject variance
      double subj_var = population->estimates.get_mass_mean_variance();
      
      // Get all subject mean masses
      arma::vec subject_means = population->get_all_subject_mass_means();
      double sum_means = arma::sum(subject_means);
      
      // Posterior parameters (conjugate normal)
      double post_prec = 1.0 / prior_var + n / subj_var;
      double post_var  = 1.0 / post_prec;
      double post_mean = post_var * (prior_mean / prior_var + sum_means / subj_var);
      
      // Draw from posterior
      population->estimates.mass_mean = R::rnorm(post_mean, sqrt(post_var));
    }

    //
    // Sample width_mean (μ_ω) only
    //
    void sample_width_mean(Population *population) {
      
      int n = population->num_subjects;
      
      // Prior parameters
      double prior_mean = population->priors.width_mean_mean;
      double prior_var  = population->priors.width_mean_var;
      
      // Subject-to-subject variance
      double subj_var = population->estimates.get_width_mean_variance();
      
      // Get all subject mean widths
      arma::vec subject_means = population->get_all_subject_width_means();
      double sum_means = arma::sum(subject_means);
      
      // Posterior parameters (conjugate normal)
      double post_prec = 1.0 / prior_var + n / subj_var;
      double post_var  = 1.0 / post_prec;
      double post_mean = post_var * (prior_mean / prior_var + sum_means / subj_var);
      
      // Draw from posterior
      population->estimates.width_mean = R::rnorm(post_mean, sqrt(post_var));
    }

};

#endif
