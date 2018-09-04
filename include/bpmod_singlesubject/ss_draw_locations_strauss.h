#ifndef GUARD_bpmod_singlesubject_draw_locations_strauss_h
#define GUARD_bpmod_singlesubject_draw_locations_strauss_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_datastructures/patient.h>
//#include <bp_mcmc/mh.h>
#include <bpmod_singlesubject/ss_draw_locations.h>
#include <bp_mcmc/utils.h>


// 
// SS_DrawLocationStrauss
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sampling the pulse locations
//

class SS_DrawLocationsStrauss : public SS_DrawLocations
{

  public:
    // Constructors
    SS_DrawLocationsStrauss(double in_pv, int in_adjust_iter, int in_max_iter,
                     double in_target_ratio, bool in_verbose, int in_verbose_iter) :
      SS_DrawLocations(in_pv, in_adjust_iter, in_max_iter, in_target_ratio,
                       in_verbose, in_verbose_iter) { };

    // Destructor
    ~SS_DrawLocationsStrauss() { }

  private:

    bool parameter_support(double val, Patient *patient) {
      return ((val <= patient->data.fitend) &&
              (val > patient->data.fitstart));
    }

    //
    // posterior_function()
    //   for strauss location prior mmh
    //
    double posterior_function(PulseIter *pulse,
                              double proposal,
                              Patient *patient) {

      // Internal variables
      PulseIter pulseiter = *pulse;
      double acceptance_ratio = 0;
      double prior_ratio = 0;
      // Extracted variables
      arma::vec curr_mean_contrib =
        pulseiter->get_mean_contribution(patient->data.time,
                                     patient->estimates.get_decay());
      double curr_likelihood = patient->likelihood(false); // would prefer get_likelihood()
      double plikelihood = 0;
      double gamma = patient->priors.strauss_repulsion;
      double current_time = 0;

      // Calculate sum_s_r for proposal value and current value
      // TODO: update arguments
      int sum_s_r_proposal = patient->calc_sr_strauss(proposal, &(*pulseiter));
      int sum_s_r_current  = patient->calc_sr_strauss(pulseiter->time, &(*pulseiter));

      // Calculate prior ratio - 
      //      gamma^(sum(s_r*)_proposal - sum(s_r)_current)
      // NOTE: updated gamma arg already
      prior_ratio = pow(gamma, sum_s_r_proposal - sum_s_r_current);

      // if prior ratio is 0 (~EPS), set acceptance_ratio to 0, 
      // else calculate it 
      // (note: necessary because log(0) = NaN, like_ratio on log scale)
      if (abs(prior_ratio) < 1.0e-42) {

        acceptance_ratio = -1e300;

      } else {

        // *** ***NOTE: TODO: this is not working*** ***
        // Save current time and set its time to proposed value
        current_time = pulseiter->time;
        pulseiter->time  = proposal;

        // Recalculate mean_contrib of this pulse assuming proposed value 
        // TODO: verify recalculation is done internally
        // It is -- likelihood() calls mean_concentration() which calls this
        // pulse->get_mean_contribution(patient->data->time,
        //                              patient->estimates->get_decay());

        // Calculate the likelihood under the proposed value 
        // returns log-likelihood  -- calculate with new pulse->time
        plikelihood = patient->likelihood(false); 

        // Calculate the likelihood ratio 
        acceptance_ratio = log(prior_ratio) + (plikelihood - curr_likelihood);

        // Reset pulse->time to current (sample() chooses whether to keep)
        //   and get_mean_contribution() will recalc that when requested.
        pulseiter->time = current_time;

      }

      return acceptance_ratio;

    }


};


#endif

