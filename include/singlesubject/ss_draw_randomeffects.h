#ifndef GUARD_ss_draw_randomeffects_h
#define GUARD_ss_draw_randomeffects_h

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

class SS_DrawRandomEffects :
  public ModifiedMetropolisHastings<PulseEstimate, Patient, double, ProposalVariance>
{

  public:
    // Constructors
    SS_DrawRandomEffects(double in_pv, // double or arma::vec
                         int in_adjust_iter,
                         int in_max_iter,
                         double in_target_ratio) :
      ModifiedMetropolisHastings <PulseEstimate, Patient, double, ProposalVariance>::
        ModifiedMetropolisHastings(in_pv,
                                   in_adjust_iter,
                                   in_max_iter,
                                   in_target_ratio) { };
    // Pulse level estimates need to be done at the pulse level
    void sample_pulses(Patient *patient, std::string measure) {

      std::list<PulseEstimate>::iterator pulse = patient->pulses.begin();
      std::list<PulseEstimate>::const_iterator pulse_end = patient->pulses.end();

      while (pulse != pulse_end) {
        // Sample pulse,
        //   note: &(*pulse) derefs iter, then gets address of underlying obj
        //if (measure == "mass") {
            sample(&(*pulse), &pulse->mass, patient);
        //} else if (measure == "width") {
        //    sample(&(*pulse), &pulse->width, patient);
        //}
        pulse++;
      }

    }


  private:

    bool parameter_support(double val) {
      // NOTE: original was:
      //   mass > 0.0 && width > 0.01 && width < 10
      return ( val > 0.0 );
    }

    //
    // posterior_function()
    //   for strauss location prior mmh
    //
    double posterior_function(PulseEstimate *pulse,
                              double proposal,
                              Patient *patient) {

        double patient_mass_mean = patient->estimates->mass;
        double patient_mass_sd   = patient->estimates->mass_sd;
        double curr_likelihood   = patient->likelihood(false);

        // Compute the log of the ratio of the priors
        prior_old = pow(pulse->mass - patient_mass_mean, 2);
        prior_old *= 0.5 * prior_old;
        prior_new = proposal - patient_mass_mean;
        prior_new *= 0.5 * prior_new;
        prior_ratio = prior_old - prior_new;
        prior_ratio /= patient_mass_sd;
        prior_ratio /= patient_mass_sd;

        // Save the current value of mass/width and set to proposed value
        current_mass = pulse->mass;
        pulse->mass = proposal;

        // Calculate likelihood assuming proposed mass/width 
        plikelihood      = patient->likelihood(false);

        // Reset pulse->time to current (sample() chooses whether to keep)
        //   and get_mean_contribution() will recalc that when requested.
        pulse->mass = current_mass;

      return (prior_ratio + (plikelihood - curr_likelihood));

    }


};


#endif

