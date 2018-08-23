#ifndef GUARD_bpmod_singlesubject_draw_locations_h
#define GUARD_bpmod_singlesubject_draw_locations_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_datastructures/patient.h>
#include <bp_mcmc/mh.h>
#include <bp_mcmc/utils.h>


// 
// SS_DrawLocations
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sampling the pulse locations
//
// TODO: Want to loop over pulses within the .sample() function for
// pulse-specific features
//

class SS_DrawLocationsStrauss :
  public ModifiedMetropolisHastings<PulseEstimates, Patient, double, ProposalVariance>
{

  public:
    // Constructors
    SS_DrawLocationsStrauss(double in_pv, // double or arma::vec
                            int in_adjust_iter,
                            int in_max_iter,
                            double in_target_ratio) :
      ModifiedMetropolisHastings
      <PulseEstimates, Patient, double,
       ProposalVariance>::ModifiedMetropolisHastings(in_pv,
                                                     in_adjust_iter,
                                                     in_max_iter,
                                                     in_target_ratio) { };
    // Pulse level estimates need to be done at the pulse level
    void sample_pulses(Patient *patient, int iter) {

      //PulseIter pulse = patient->pulses.begin();
      //PulseConstIter pulse_end = patient->pulses.end();

      //while (pulse != pulse_end) {
      //  // Sample pulse,
      //  //   note: &(*pulse) derefs iter, then gets address of underlying obj
      //  sample(&(*pulse), &pulse->time, patient, iter);
      //  pulse++;
      //}
      for (auto &pulse : patient->pulses) {
        //sample(&(*pulse), &pulse->time, patient, iter);
        sample(&pulse, &(pulse.time), patient, iter);
      }

    }


  private:

    bool parameter_support(double val, Patient *patient) {
      return ((val <= patient->data.fitend) &&
              (val > patient->data.fitstart));
    }

    //
    // posterior_function()
    //   for strauss location prior mmh
    //
    double posterior_function(PulseEstimates *pulse,
                              double proposal,
                              Patient *patient) {

      // Calculated tmp variables
      double acceptance_ratio;
      double prior_ratio;
      // Extracted variables
      arma::vec curr_mean_contrib =
        pulse->get_mean_contribution(patient->data.time,
                                     patient->estimates.get_decay());
      double curr_likelihood = patient->likelihood(false); // would prefer get_likelihood()
      double plikelihood;
      double gamma = patient->priors.strauss_repulsion;
      double current_time;

      //std::cout << "Current mean contrib=" << curr_mean_contrib <<
      //  "; Current likelihood=" << curr_likelihood << std::endl;

      // Calculate sum_s_r for proposal value and current value
      // TODO: update arguments
      int sum_s_r_proposal = patient->calc_sr_strauss(proposal, pulse);
      int sum_s_r_current  = patient->calc_sr_strauss(pulse->time, pulse);

      // Calculate prior ratio - 
      //      gamma^(sum(s_r*)_proposal - sum(s_r)_current)
      // NOTE: updated gamma arg already
      prior_ratio = pow(gamma, sum_s_r_proposal - sum_s_r_current);

      // NOTE: prior ratio is working fine!
      std::cout << "prior ratio=" << prior_ratio <<
        "; gamma=" << gamma <<
        "; sum_s_r_proposal=" << sum_s_r_proposal <<
        "; sum_s_r_current=" << sum_s_r_current <<
        std::endl;

      // if prior ratio is 0 (~EPS), set acceptance_ratio to 0, 
      // else calculate it 
      // (note: necessary because log(0) = NaN, like_ratio on log scale)
      if (fabs(prior_ratio) < 1.0e-42) {

        acceptance_ratio = -1e300;

      } else {

        // *** ***NOTE: TODO: this is not working*** ***
        // Save current time and set its time to proposed value
        current_time = pulse->time;
        pulse->time  = proposal;

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
        pulse->time = current_time;

      }

      return acceptance_ratio;

    }


};


#endif

