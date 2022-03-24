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
    SS_DrawRandomEffects(double in_pv,
                         int in_adjust_iter,
                         int in_max_iter,
                         double in_target_ratio,
                         bool for_width,
                         bool verbose,
                         int verbose_iter) :
      ModifiedMetropolisHastings <PulseEstimates, Patient, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {

        // Choose which set of parameters to use: width or mass
        if (for_width) {
          est_mean_       = &PatientEstimates::width_mean;
          est_prec_         = &PatientEstimates::width_prec;
          randomeffect_   = &PulseEstimates::width;
          parameter_name = "pulse width";
        } else {
          est_mean_       = &PatientEstimates::mass_mean;
          est_prec_         = &PatientEstimates::mass_prec;
          randomeffect_   = &PulseEstimates::mass;
          parameter_name = "pulse mass";
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
    double PatientEstimates::*est_prec_;
    double PulseEstimates::*randomeffect_; //pulse specific mass or width

    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; };

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
        double patient_prec      = (*est).*est_prec_;
        double curr_likelihood = patient->likelihood(false);

        //Rcpp::Rcout << "Patient mean: " << patient_mean <<
        //  "; Patient SD: " << patient_sd << 
        //  "; Current likelihood: " << curr_likelihood <<
        //  "; Proposal: " << proposal <<
        //  std::endl;

        // Compute the log of the ratio of the priors
        prior_old    = (*pulse).*randomeffect_ - patient_mean;
        prior_old   *= 0.5 * prior_old;
        prior_new    = proposal - patient_mean;
        prior_new   *= 0.5 * prior_new;
        prior_ratio  = prior_old - prior_new;
        prior_ratio /= 1/sqrt(patient_prec);
        prior_ratio /= 1/sqrt(patient_prec);

        // Save the current value of mass/width and set to proposed value
        //std::cout << "\n\nInitial random effect value: " << (*pulse).*randomeffect_ << std::endl;
        current_randomeffect    = (*pulse).*randomeffect_;
        //std::cout << "Saved initial random effect value: " << current_randomeffect << std::endl;
        (*pulse).*randomeffect_ = proposal;
        //std::cout << "New random effect value: " << proposal << std::endl;
        //std::cout << "New (in-place) random effect value: " << (*pulse).*randomeffect_ << std::endl;

        // Calculate likelihood assuming proposed mass/width 
        plikelihood          = patient->likelihood(false);

        //Rcpp::Rcout << "Patient mean: " << patient_mean << " Patient_sd: " << patient_sd << "\n"
        //            << "Current RE: " << current_randomeffect << " New RE: " << proposal
        //            << "\nPrior ratio: " << prior_ratio
        //            << "\npLikelihood: " << plikelihood << " currLikelihood: " << curr_likelihood
        //            << "\n";

        // Reset pulse->time to current (sample() chooses whether to keep)
        //   and get_mean_contribution() will recalc that when requested.
        (*pulse).*randomeffect_ = current_randomeffect;

      return (prior_ratio + (plikelihood - curr_likelihood));

    }


};


#endif

