#ifndef GUARD_ss_drawfixedeffects_h
#define GUARD_ss_drawfixedeffects_h

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
      double prior_mass_mean   = patient->priors->mass_mean;
      double prior_mass_var    = patient->priors->mass_variance;
      double curr_mass         = patient->estimates->mass_mean;
      double curr_mass_sd      = patient->estimates->mass_sd;
      std::list<PulseEstimate>::const_iterator pulse     = patient->pulses.begin();
      std::list<PulseEstimate>::const_iterator pulse_end = patient->pulses.end();

      // Prior Ratio
      prior_ratio = (pow(curr_mass - prior_mass_mean, 2) -
                     pow(proposal - prior_mass_mean, 2)) /
                    (2 * prior_mass_var);

      // 'likelihood' ratio -- Ratio of p(alpha|mu, nu, kappa)
      while (pulse != pulse_end) {
        psum_old += pow(pulse->mass - curr_mass, 2) * pulse->tvarscale_mass;
        psum_new += pow(pulse->mass - proposal, 2) * pulse->tvarscale_mass;

        // Normalizing constants
        oldint += Rf_pnorm5(curr_mass * sqrt(pulse->tvarscale_mass) / curr_mass_sd,
                            0, 1, 1.0, 1.0); // second 1.0 does the log xform for us
        newint += Rf_pnorm5(proposal * sqrt(pulse->tvarscale_mass) / curr_mass_sd,
                            0, 1, 1.0, 1.0); // first 1.0 says to use lower tail
        ++pulse;
      }

      prop_ratio = 0.5 / pow(curr_mass_sd, 2) * (psum_old - psum_new);
      normalizing_ratio = oldint - newint;

      return prior_ratio + prop_ratio + normalizing_ratio;

    }

};

#endif

