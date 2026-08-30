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
          est_sd_         = &PatientEstimates::width_sd;
          randomeffect_   = &PulseEstimates::width;
          tvarscale_      = &PulseEstimates::tvarscale_width;
          parameter_name = "pulse width";
        } else {
          est_mean_       = &PatientEstimates::mass_mean;
          est_sd_         = &PatientEstimates::mass_sd;
          randomeffect_   = &PulseEstimates::mass;
          tvarscale_      = &PulseEstimates::tvarscale_mass;
          parameter_name = "pulse mass";
        }

      };

    // Pulse-specific estimate -- this function samples for each pulse
    void sample_pulses(Patient *patient, int iter) { //, std::string measure) {

      for (auto &pulse : patient->pulses) {
        sample(&pulse, &(pulse.*randomeffect_), patient, iter);
      }

    }

    //
    // posterior_function()
    //   Log acceptance ratio for the pulse-specific random effect (mass or
    //   width). The random effect is a scale-mixture-of-normals (Student-t):
    //     theta_k ~ N(mu, sigma^2 / kappa_k)  truncated at 0,
    //   where kappa_k = tvarscale is the per-pulse t-scale. The prior quadratic
    //   term must therefore be weighted by kappa_k -- exactly as the tvarscale
    //   draw (ss_draw_tvarscale.h) and the SD draws (ss_draw_sdrandomeffects.h,
    //   pop_draw_pulse_sds.h) already do. Omitting kappa here made the three
    //   full conditionals mutually inconsistent, so the sampler did not target a
    //   coherent posterior; the effect was worst for the weakly-identified pulse
    //   width, whose pulse-to-pulse SD then wandered/space-filled. The symmetric
    //   random-walk proposal means the truncation normalizing constant cancels
    //   between current and proposal, so it need not appear here.
    //
    //   Exposed publicly so the closed-form log-ratio can be unit tested; the MH
    //   base dispatches through its own private virtual, so this does not change
    //   normal sampling behavior.
    double posterior_function(PulseEstimates *pulse,
                              double proposal,
                              Patient *patient) {

        double prior_old, prior_new, prior_ratio, current_randomeffect,
               plikelihood;
        PatientEstimates *est  = &patient->estimates;
        double patient_mean    = (*est).*est_mean_;
        double patient_sd      = (*est).*est_sd_;
        double kappa           = (*pulse).*tvarscale_;
        double curr_likelihood = patient->likelihood(false);

        //Rcpp::Rcout << "Patient mean: " << patient_mean <<
        //  "; Patient SD: " << patient_sd <<
        //  "; Current likelihood: " << curr_likelihood <<
        //  "; Proposal: " << proposal <<
        //  std::endl;

        // Compute the log of the ratio of the priors. The per-pulse t-scale
        // kappa weights the quadratic deviation (prior variance is sigma^2/kappa).
        double current_re = (*pulse).*randomeffect_;
        if (patient->lognormal_pulses) {
          // Log-normal parameterization (papers): log(theta) ~ N(mu, sigma^2/kappa).
          // patient_mean/patient_sd are the mean/SD of the LOG random effect. The
          // MH proposal is a symmetric random walk on the natural-scale value, so
          // the target density of theta carries the 1/theta Jacobian of the
          // log-normal; its log-ratio contributes (log(theta) - log(theta')).
          // Both current_re and proposal are > 0 (parameter_support enforces it),
          // so the logs are well defined. No truncation constant is needed.
          double logold = log(current_re);
          double lognew = log(proposal);
          prior_old   = 0.5 * (logold - patient_mean) * (logold - patient_mean);
          prior_new   = 0.5 * (lognew - patient_mean) * (lognew - patient_mean);
          prior_ratio = kappa * (prior_old - prior_new) / (patient_sd * patient_sd);
          prior_ratio += (logold - lognew);  // log-normal Jacobian
        } else {
          // Natural-scale truncated-normal parameterization (research option).
          prior_old   = current_re - patient_mean;
          prior_old  *= 0.5 * prior_old;
          prior_new   = proposal - patient_mean;
          prior_new  *= 0.5 * prior_new;
          prior_ratio = kappa * (prior_old - prior_new) / (patient_sd * patient_sd);
        }

        // Save the current value of mass/width and set to proposed value
        //std::cout << "\n\nInitial random effect value: " << (*pulse).*randomeffect_ << std::endl;
        current_randomeffect    = (*pulse).*randomeffect_;
        //std::cout << "Saved initial random effect value: " << current_randomeffect << std::endl;
        (*pulse).*randomeffect_ = proposal;
        //std::cout << "New random effect value: " << proposal << std::endl;
        //std::cout << "New (in-place) random effect value: " << (*pulse).*randomeffect_ << std::endl;

        // Calculate likelihood assuming proposed mass/width 
        plikelihood          = patient->likelihood(false);

        // Reset pulse->time to current (sample() chooses whether to keep)
        //   and get_mean_contribution() will recalc that when requested.
        (*pulse).*randomeffect_ = current_randomeffect;

      return (prior_ratio + (plikelihood - curr_likelihood));

    }


  private:

    double PatientEstimates::*est_mean_;
    double PatientEstimates::*est_sd_;
    double PulseEstimates::*randomeffect_; //pulse specific mass or width
    double PulseEstimates::*tvarscale_;    //pulse specific t-scale (kappa)

    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; };

    bool parameter_support(double val, Patient *patient) {
      // NOTE: original was:
      //   mass > 0.0 && width > 0.01 && width < 10
      return ( val > 0.0 );
    }

};


#endif

