#ifndef GUARD_ss_draw_fixedeffects_h
#define GUARD_ss_draw_fixedeffects_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"


//
// SS_DrawFixedEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the mean mass & width
//

class SS_DrawFixedEffects : public ModifiedMetropolisHastings<Patient, bool, double, ProposalVariance>
{

  public:

    //
    // Constructor
    //   pass the proposal variance parameters to the constructor
    //
    SS_DrawFixedEffects(double in_pv, // double or arma::vec
                        int in_adjust_iter,
                        int in_max_iter,
                        double in_target_ratio,
                        bool for_width) :
      ModifiedMetropolisHastings
      <Patient, bool, double,
       ProposalVariance>::ModifiedMetropolisHastings(in_pv,
                                                     in_adjust_iter,
                                                     in_max_iter,
                                                     in_target_ratio) { 
         //for_width = in_for_width;
         if (for_width) {
           prior_mean_     = &PatientPriors::width_mean;
           prior_variance_ = &PatientPriors::width_variance;
           est_mean_       = &PatientEstimates::width_mean;
           est_sd_         = &PatientEstimates::width_sd;
           tvarscale_      = &PulseEstimate::tvarscale_width;
           randomeffect_   = &PulseEstimate::width;
         } else {
           prior_mean_     = &PatientPriors::mass_mean;
           prior_variance_ = &PatientPriors::mass_variance;
           est_mean_       = &PatientEstimates::mass_mean;
           est_sd_         = &PatientEstimates::mass_sd;
           tvarscale_      = &PulseEstimate::tvarscale_mass;
           randomeffect_   = &PulseEstimate::mass;
         }
       };

  private:

    double PatientPriors::*prior_mean_;
    double PatientPriors::*prior_variance_;
    double PatientEstimates::*est_mean_;
    double PatientEstimates::*est_sd_;
    double PulseEstimate::*tvarscale_;
    double PulseEstimate::*randomeffect_; //pulse specific mass or width

    bool parameter_support(double val, bool *notused) { return (val > 0.0); }

    double posterior_function(Patient *patient, double proposal, bool *notused) {

      double prior_ratio       = 0.0 ;
      double psum_old          = 0.0 ;
      double psum_new          = 0.0 ;
      double newint            = 0.0 ;
      double oldint            = 0.0 ;
      double normalizing_ratio = 0.0 ;
      double prop_ratio        = 0.0 ;
      PatientPriors *priors    = patient->priors;
      PatientEstimates *est    = patient->estimates;
      double prior_mass_mean   = (*priors).*prior_mean_;
      double prior_mass_var    = (*priors).*prior_variance_;
      double current           = (*est).*est_mean_;
      double stddev            = (*est).*est_sd_;
      double scale             = 0.0;
      double randomeffect      = 0.0;


      // Prior Ratio
      prior_ratio = (pow(current - prior_mass_mean, 2) -
                     pow(proposal - prior_mass_mean, 2)) /
                    (2 * prior_mass_var);

      // 'likelihood' ratio -- Ratio of p(alpha|mu, nu, kappa) for current and
      // proposed alpha
      for (auto &pulse : patient->pulses) {  // uses range based loop instead of iterators

        scale        = pulse.*tvarscale_;
        randomeffect = pulse.*randomeffect_;
        psum_old += pow(randomeffect - current, 2) * scale;
        psum_new += pow(randomeffect - proposal, 2) * scale;

        // Normalizing constants
        oldint += Rf_pnorm5(current * sqrt(scale) / stddev,
                            0, 1, 1.0, 1.0); // second 1.0 does the log xform for us
        newint += Rf_pnorm5(proposal * sqrt(scale) / stddev,
                            0, 1, 1.0, 1.0); // first 1.0 says to use lower tail

      }

      prop_ratio = 0.5 / pow(stddev, 2) * (psum_old - psum_new);
      normalizing_ratio = oldint - newint;

      return prior_ratio + prop_ratio + normalizing_ratio;

    }

};

#endif

