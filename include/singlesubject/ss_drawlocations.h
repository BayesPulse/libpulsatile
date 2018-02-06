#ifndef GUARD_ss_drawlocations_h
#define GUARD_ss_drawlocations_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"
#include "utils.h"


// 
// SS_DrawLocations
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sampling the pulse locations
//
// TODO: Want to loop over pulses within the .sample() function for
// pulse-specific features
//

class SS_DrawLocationsStrauss : public ModifiedMetropolisHastings<PulseEstimate, double, ProposalVariance>
{

  public:
    // Constructors

  private:

    bool parameter_support(double val, double fitstart, double fitend) {
      return ((val <= fitend) && (val > fitstart));
    }

    //
    // posterior_function()
    //   for strauss location prior mmh
    //
    double posterior_function(PulseEstimate *pulse, double proposal) {

      //
      double acceptance_ratio;
      double prior_ratio;

      // Calculate sum_s_r for proposal value and current value
      int sum_s_r_proposal = calc_sr_strauss(proposal, list, pulse, priors);
      int sum_s_r_current  = calc_sr_strauss(pulse->time, list, pulse, priors);

      // Calculate prior ratio - 
      //      gamma^(sum(s_r*)_proposal - sum(s_r)_current)
      prior_ratio = pow(priors->gamma, sum_s_r_proposal - sum_s_r_current);

      // if prior ratio is 0 (~EPS), set acceptance_ratio to 0, 
      // else calculate it 
      // (note: necessary because log(0) = NaN, like_ratio on log scale)
      if (fabs(prior_ratio) < EPS) {

        acceptance_ratio = -1e300;

      } else {

        // Save current time and set its time to proposed value
        current_time = pulse->time;
        pulse->time  = proposal;

        // Save mean_contrib of this pulse 
        for (i = 0; i < N; i++) {
          tmp_mean_contrib[i] = pulse->mean_contrib[i];
        }

        // Recalculate mean_contrib of this pulse assuming proposed value 
        mean_contribution(pulse, ts, parms, N);

        // Calculate the likelihood under the proposed value 
        // returns log-likelihood 
        plikelihood = likelihood(list, ts, parms, N, list);

        // Calculate the likelihood ratio 
        acceptance_ratio = log(prior_ratio) + (plikelihood - *like);

      }

      return acceptance_ratio;

    }

}


#endif

