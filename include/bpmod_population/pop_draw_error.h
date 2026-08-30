#ifndef GUARD_bpmod_population_draw_error_h
#define GUARD_bpmod_population_draw_error_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_datastructures/population.h>

//
// Pop_DrawError
//   Draw the population-level model error variance using Gibbs sampling
//   
//   Updates: population.estimates.errorsq (σ²_e)
//
//   Prior: σ² ~ InvGamma(α, β)
//   Data:  Observations across all subjects
//   Posterior: σ² | data ~ InvGamma(α_n, β_n)
//   
//   Where:
//     α_n = α + N_total / 2
//     β_n = β + SSE_total / 2
//     N_total = total observations across all subjects
//     SSE_total = sum of squared errors across all subjects
//

class Pop_DrawError 
{

  public:
    // Constructors
    Pop_DrawError() { }

    void sample(Population *population) {

      // Get prior parameters
      double alpha = population->priors.error_alpha;
      double beta  = population->priors.error_beta;

      // Calculate total observations and sum of squared errors
      int N_total = population->get_total_observations();
      double ssq_total = 0.0;
      
      for (auto& subject : population->subjects) {
        ssq_total += subject.get_sumerrorsquared(false);
      }

      // Posterior parameters for inverse-gamma
      double post_alpha = alpha + N_total / 2.0;
      double post_beta  = 1.0 / (1.0 / beta + 0.5 * ssq_total);
      
      // Draw from inverse-gamma: σ² = 1 / Gamma(α, β)
      // Note: Rf_rgamma uses shape and scale parameterization
      population->estimates.errorsq = 1.0 / Rf_rgamma(post_alpha, post_beta);
    }

};

#endif
