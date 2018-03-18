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
                         double in_target_ratio,
                         bool for_width) :
      ModifiedMetropolisHastings <PulseEstimate, Patient, double, ProposalVariance>::
        ModifiedMetropolisHastings(in_pv,
                                   in_adjust_iter,
                                   in_max_iter,
                                   in_target_ratio) { 

         // Choose which set of parameters to use: width or mass
          if (for_width) {
            //prior_mean_     = &PatientPriors::width_mean;
            //prior_variance_ = &PatientPriors::width_variance;
            est_mean_       = &PatientEstimates::width_mean;
            est_sd_         = &PatientEstimates::width_sd;
            tvarscale_      = &PulseEstimate::tvarscale_width;
            randomeffect_   = &PulseEstimate::width;
          } else {
            //prior_mean_     = &PatientPriors::mass_mean;
            //prior_variance_ = &PatientPriors::mass_variance;
            est_mean_       = &PatientEstimates::mass_mean;
            est_sd_         = &PatientEstimates::mass_sd;
            tvarscale_      = &PulseEstimate::tvarscale_mass;
            randomeffect_   = &PulseEstimate::mass;
          }

        };

    // Pulse-specific estimate -- this function samples for each pulse
    void sample_pulses(Patient *patient) { //, std::string measure) {

      for (auto &pulse : patient->pulses) {
        sample(&(*pulse), (*pulse).*randomeffect_, patient);
      }

    }


  private:

    //double PatientPriors::*prior_mean_;
    //double PatientPriors::*prior_variance_;
    double PatientEstimates::*est_mean_;
    double PatientEstimates::*est_sd_;
    double PulseEstimate::*tvarscale_;
    double PulseEstimate::*randomeffect_; //pulse specific mass or width

    bool parameter_support(double val, Patient *patient) {
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

        double prior_old, prior_new, prior_ratio, current_mass, plikelihood;
        PatientPriors *priors = patient->priors;
        PatientEstimates *est = patient->estimates;
        double patient_mean = (*est).*est_mean_;
        double patient_sd   = (*est).*est_sd_;
        double curr_likelihood   = patient->likelihood(false);

        // Compute the log of the ratio of the priors
        prior_old = pow(pulse.*randomeffect_ - patient_mean, 2);
        prior_old *= 0.5 * prior_old;
        prior_new = proposal - patient_mass_mean;
        prior_new *= 0.5 * prior_new;
        prior_ratio = prior_old - prior_new;
        prior_ratio /= patient_mass_sd;
        prior_ratio /= patient_mass_sd;

        // Save the current value of mass/width and set to proposed value
        current_mass = pulse->mass;
        pulse->mass  = proposal;

        // Calculate likelihood assuming proposed mass/width 
        plikelihood = patient->likelihood(false);

        // Reset pulse->time to current (sample() chooses whether to keep)
        //   and get_mean_contribution() will recalc that when requested.
        pulse->mass = current_mass;

      return (prior_ratio + (plikelihood - curr_likelihood));

    }


};


#endif

