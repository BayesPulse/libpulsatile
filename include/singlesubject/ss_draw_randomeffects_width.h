#ifndef GUARD_ss_draw_randomeffectswidths_h
#define GUARD_ss_draw_randomeffectswidths_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"
//#include "utils.h"


// 
// SS_DrawRandomEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sampling the individual pulse mean and standard deviations
//

class SS_DrawRandomEffectsWidths :
  public ModifiedMetropolisHastings<PulseEstimate, Patient, double, ProposalVariance>
{

  public:
    // Constructors
    SS_DrawRandomEffectsWidths(double in_pv, // double or arma::vec
                         int in_adjust_iter,
                         int in_max_iter,
                         double in_target_ratio) :
      ModifiedMetropolisHastings <PulseEstimate, Patient, double, ProposalVariance>::
        ModifiedMetropolisHastings(in_pv,
                                   in_adjust_iter,
                                   in_max_iter,
                                   in_target_ratio) { };

    // Pulse-specific estimate -- this function samples for each pulse
    void sample_pulses(Patient *patient) { //, std::string measure) {

      PulseIter pulse = patient->pulses.begin();
      PulseConstIter pulse_end = patient->pulses.end();

      while (pulse != pulse_end) {
        // Sample pulse,
        //   note: &(*pulse) derefs iter, then gets address of underlying obj
        //if (measure == "width") {
            sample(&(*pulse), &pulse->width, patient);
        //} else if (measure == "width") {
        //    sample(&(*pulse), &pulse->width, patient);
        //}
        pulse++;
      }

    }


  private:

    bool parameter_support(double val, Patient *patient) {
      // NOTE: original was:
      //   width > 0.0 && width > 0.01 && width < 10
      return ( val > 0.0 );
    }

    //
    // posterior_function()
    //   for strauss location prior mmh
    //
    double posterior_function(PulseEstimate *pulse,
                              double proposal,
                              Patient *patient) {

        double prior_old, prior_new, prior_ratio, current_width, plikelihood;
        double patient_width_mean = patient->estimates->width_mean;
        double patient_width_sd   = patient->estimates->width_sd;
        double curr_likelihood   = patient->likelihood(false);

        // Compute the log of the ratio of the priors
        prior_old = pow(pulse->width - patient_width_mean, 2);
        prior_old *= 0.5 * prior_old;
        prior_new = proposal - patient_width_mean;
        prior_new *= 0.5 * prior_new;
        prior_ratio = prior_old - prior_new;
        prior_ratio /= patient_width_sd;
        prior_ratio /= patient_width_sd;

        // Save the current value of width/width and set to proposed value
        current_width = pulse->width;
        pulse->width  = proposal;

        // Calculate likelihood assuming proposed width/width 
        plikelihood = patient->likelihood(false);

        // Reset pulse->time to current (sample() chooses whether to keep)
        //   and get_mean_contribution() will recalc that when requested.
        pulse->width = current_width;

      return (prior_ratio + (plikelihood - curr_likelihood));

    }


};


#endif

