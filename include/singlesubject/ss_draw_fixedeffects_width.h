#ifndef GUARD_ss_draw_fixedeffectswidths_h
#define GUARD_ss_draw_fixedeffectswidths_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"


//
// SS_DrawFixedEffectsWidths
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the mean width
//

class SS_DrawFixedEffectsWidths : public ModifiedMetropolisHastings<Patient, bool, double, ProposalVariance>
{

  public:

    //
    // Constructor
    //   pass the proposal variance parameters to the constructor
    //
    SS_DrawFixedEffectsWidths(double in_pv, // double or arma::vec
                        int in_adjust_iter,
                        int in_max_iter,
                        double in_target_ratio) :
      ModifiedMetropolisHastings
      <Patient, bool, double,
       ProposalVariance>::ModifiedMetropolisHastings(in_pv,
                                                     in_adjust_iter,
                                                     in_max_iter,
                                                     in_target_ratio) { };

  private:

    bool parameter_support(double val, bool *notused) { return (val > 0.0); }

    double posterior_function(Patient *patient, double proposal, bool *notused) {

      double prior_ratio       = 0.0 ;
      double psum_old          = 0.0 ;
      double psum_new          = 0.0 ;
      double newint            = 0.0 ;
      double oldint            = 0.0 ;
      double normalizing_ratio = 0.0 ;
      double prop_ratio        = 0.0 ;
      double prior_width_mean   = patient->priors->width_mean;
      double prior_width_var    = patient->priors->width_variance;
      double curr_width         = patient->estimates->width_mean;
      double curr_width_sd      = patient->estimates->width_sd;
      std::list<PulseEstimate>::const_iterator pulse     = patient->pulses.begin();
      std::list<PulseEstimate>::const_iterator pulse_end = patient->pulses.end();

      // Prior Ratio
      prior_ratio = (pow(curr_width - prior_width_mean, 2) -
                     pow(proposal - prior_width_mean, 2)) /
                    (2 * prior_width_var);

      // 'likelihood' ratio -- Ratio of p(alpha|mu, nu, kappa)
      while (pulse != pulse_end) {
        psum_old += pow(pulse->width - curr_width, 2) * pulse->tvarscale_width;
        psum_new += pow(pulse->width - proposal, 2) * pulse->tvarscale_width;

        // Normalizing constants
        oldint += Rf_pnorm5(curr_width * sqrt(pulse->tvarscale_width) / curr_width_sd,
                            0, 1, 1.0, 1.0); // second 1.0 does the log xform for us
        newint += Rf_pnorm5(proposal * sqrt(pulse->tvarscale_width) / curr_width_sd,
                            0, 1, 1.0, 1.0); // first 1.0 says to use lower tail
        ++pulse;
      }

      prop_ratio = 0.5 / pow(curr_width_sd, 2) * (psum_old - psum_new);
      normalizing_ratio = oldint - newint;

      return prior_ratio + prop_ratio + normalizing_ratio;

    }

};

#endif

