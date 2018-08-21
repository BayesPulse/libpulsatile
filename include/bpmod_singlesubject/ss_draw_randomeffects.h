#ifndef GUARD_bpmod_singlesubject_draw_randomeffects_h
#define GUARD_bpmod_singlesubject_draw_randomeffects_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>
//#include "utils.h"


// 
// SS_DrawRandomEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sampling the individual pulse mean and standard deviations
//

class SS_DrawRandomEffects :
  public ModifiedMetropolisHastings<PulseEstimates, Patient, double, ProposalVariance>
{

  public:

    // Constructors
    SS_DrawRandomEffects(double in_pv, // double or arma::vec
                         int in_adjust_iter,
                         int in_max_iter,
                         double in_target_ratio,
                         bool for_width) :
      ModifiedMetropolisHastings <PulseEstimates, Patient, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv,
                                 in_adjust_iter,
                                 in_max_iter,
                                 in_target_ratio) { 

        // Choose which set of parameters to use: width or mass
        if (for_width) {
          est_mean_       = &PatientEstimates::width_mean;
          est_sd_         = &PatientEstimates::width_sd;
          randomeffect_   = &PulseEstimates::width;
        } else {
          est_mean_       = &PatientEstimates::mass_mean;
          est_sd_         = &PatientEstimates::mass_sd;
          randomeffect_   = &PulseEstimates::mass;
        }

      };

    // Pulse-specific estimate -- this function samples for each pulse
    void sample_pulses(Patient *patient, int iter) { //, std::string measure) {

      for (auto &pulse : patient->pulses) {
        sample(&pulse, &(pulse.*randomeffect_), patient, iter);
      }

    }


  private:

    double PatientEstimates::*est_mean_;
    double PatientEstimates::*est_sd_;
    double PulseEstimates::*randomeffect_; //pulse specific mass or width

    bool parameter_support(double val, Patient *patient) {
      // NOTE: original was:
      //   mass > 0.0 && width > 0.01 && width < 10
      return ( val > 0.0 );
    }

    //
    // posterior_function()
    //   for strauss location prior mmh
    //
    double posterior_function(PulseEstimates *pulse,
                              double proposal,
                              Patient *patient) {

        double prior_old, prior_new, prior_ratio, current_randomeffect, 
               plikelihood;
        PatientEstimates *est  = &patient->estimates;
        double patient_mean    = (*est).*est_mean_;
        double patient_sd      = (*est).*est_sd_;
        double curr_likelihood = patient->likelihood(false);

        // Compute the log of the ratio of the priors
        prior_old = pow((*pulse).*randomeffect_ - patient_mean, 2);
        prior_old *= 0.5 * prior_old;
        prior_new = proposal - patient_mean;
        prior_new *= 0.5 * prior_new;
        prior_ratio = prior_old - prior_new;
        prior_ratio /= patient_sd;
        prior_ratio /= patient_sd;

        // Save the current value of mass/width and set to proposed value
        current_randomeffect    = (*pulse).*randomeffect_;
        (*pulse).*randomeffect_ = proposal;

        // Calculate likelihood assuming proposed mass/width 
        plikelihood          = patient->likelihood(false);

        // Reset pulse->time to current (sample() chooses whether to keep)
        //   and get_mean_contribution() will recalc that when requested.
        (*pulse).*randomeffect_ = current_randomeffect;

      return (prior_ratio + (plikelihood - curr_likelihood));

    }


};


#endif

