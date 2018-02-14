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

class SS_DrawLocationsStrauss :
  public ModifiedMetropolisHastings<PulseEstimate, Patient, double, ProposalVariance>
{

  public:
    // Constructors
    SS_DrawLocationsStrauss(double in_pv, // double or arma::vec
                            int in_adjust_iter,
                            int in_max_iter,
                            double in_target_ratio) :
      ModifiedMetropolisHastings
      <PulseEstimate, Patient, double,
       ProposalVariance>::ModifiedMetropolisHastings(in_pv, in_adjust_iter,
                                                     in_max_iter,
                                                     in_target_ratio) { };
    // Pulse level estimates need to be done at the pulse level
    void sample_pulses(Patient *patient) {

      std::list<PulseEstimate>::iterator pulse = patient->pulses.begin();
      std::list<PulseEstimate>::const_iterator pulse_end = patient->pulses.end();
      while (pulse != pulse_end) {
        sample(&(*pulse), &pulse->time); // &(*pulse) derefs iter, then gets address of underlying obj
        pulse++;
      }

    }


  private:

    bool parameter_support(double val, Patient *patient) {
      return ((val <= patient->data->fitend) && (val > patient->data->fitstart));
    }

    //
    // posterior_function()
    //   for strauss location prior mmh
    //
    double posterior_function(PulseEstimate *pulse,
                              Patient *patient,
                              double proposal) {

      // Calculated tmp variables
      double acceptance_ratio;
      double prior_ratio;
      // Extracted variables
      arma::vec curr_mean_contrib = 
        pulse->get_mean_contribution(patient->data->time,
                                     patient->estimates->get_decay());
      double curr_likelihood = patient->likelihood(false); // would prefer get_likelihood()
      double plikelihood;
      double gamma = patient->priors->strauss_repulsion;
      double current_time;

      // Calculate sum_s_r for proposal value and current value
      // TODO: update arguments
      int sum_s_r_proposal = calc_sr_strauss(proposal, patient, pulse);
      int sum_s_r_current  = calc_sr_strauss(pulse->time, patient, pulse);

      // Calculate prior ratio - 
      //      gamma^(sum(s_r*)_proposal - sum(s_r)_current)
      // NOTE: updated gamma arg already
      prior_ratio = pow(gamma, sum_s_r_proposal - sum_s_r_current);

      // if prior ratio is 0 (~EPS), set acceptance_ratio to 0, 
      // else calculate it 
      // (note: necessary because log(0) = NaN, like_ratio on log scale)
      if (fabs(prior_ratio) < 1.0e-42) {

        acceptance_ratio = -1e300;

      } else {

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

    // calc_sr_strauss()
    //   Calculates sum(S(R)), the exponent on the gamma parameter in the Strauss
    //   process/prior for pulse location. Used for Strauss prior in birth_death
    //   and mh_time_strauss.
    int calc_sr_strauss(double location,     // location to test pulses against
                        Patient *patient,
                        PulseEstimate *pulse_excluded) {

      int s_r = 0;       // Sum of indicators where diff < 20
      double difference; // Time difference
      std::list<PulseEstimate>::iterator pulse = patient->pulses.begin();
      std::list<PulseEstimate>::const_iterator pulse_end = patient->pulses.end();

      while (pulse != pulse_end) {
        if (&(*pulse) != pulse_excluded) { // TODO: Test that pulse is actually excluded!
          // skip if node is same that location is from;
          difference = fabs(location - pulse->time);
          // increment by 1 if diff<R
          s_r = (difference < patient->priors->strauss_repulsion_range) ? s_r + 1 : s_r; 
        }
        pulse++;
      }

      // sum(S(R)) - scalar value for sum of # pulses too close to each other
      return(s_r); 

    }



};


#endif

