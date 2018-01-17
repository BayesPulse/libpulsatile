#ifndef GUARD_ss_drawfixedeffects_h
#define GUARD_ss_drawfixedeffects_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"


// draw_fixed_effects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the mean mass & width
//

class SS_DrawFixedEffects : public ModifiedMetropolisHastings<Patient, double, ProposalVariance>
{

  public:
    ProposalVariance pv;

    SS_DrawFixedEffects(double in_pv, int in_adjust_iter, int in_max_iter, double in_target_ratio) :
      ModifiedMetropolisHastings<Patient, double, ProposalVariance>::ModifiedMetropolisHastings() { 
        std::cout << "inpv is = " << in_pv << std::endl;
        ProposalVariance pv(in_pv, in_adjust_iter, in_max_iter, in_target_ratio);
      };

  private:
    bool const parameter_support(double val) { return (val > 0.0); }

    double const posterior_function(Patient *patient, double proposal) {

      //  TODO: Double check whether mass_variance is an SD or actually variance and make patientpriors match.
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

      //std::cout << "current mass is = " << curr_mass << std::endl;
      //std::cout << "proposal mass is = " << proposal << std::endl;

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

