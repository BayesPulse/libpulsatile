#ifndef GUARD_bpmod_singlesubject_draw_fixedeffects_h
#define GUARD_bpmod_singlesubject_draw_fixedeffects_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>


//
// SS_DrawFixedEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the mean mass & width
//

class SS_DrawFixedEffects :
  public ModifiedMetropolisHastings<Patient, bool, double, ProposalVariance>
{

  public:

    // Constructor
    //   lognormal: when true the pulse random effects are modeled on the LOG
    //   scale (papers' parameterization). The MMH container type here is `bool`
    //   (not `Patient`), so parameter_support cannot reach patient->lognormal_pulses;
    //   the flag is therefore threaded in through the constructor and stored as a
    //   member. Defaults to false to preserve the natural-scale behavior at call
    //   sites that do not pass it.
    SS_DrawFixedEffects(double in_pv,
                        int in_adjust_iter,
                        int in_max_iter,
                        double in_target_ratio,
                        bool for_width,
                        bool verbose,
                        int verbose_iter,
                        bool lognormal = false) :
      ModifiedMetropolisHastings <Patient, bool, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {

        lognormal_ = lognormal;

        // Choose which set of parameters to use: width or mass
        if (for_width) {
          prior_mean_     = &PatientPriors::width_mean;
          prior_variance_ = &PatientPriors::width_variance;
          est_mean_       = &PatientEstimates::width_mean;
          est_sd_         = &PatientEstimates::width_sd;
          tvarscale_      = &PulseEstimates::tvarscale_width;
          randomeffect_   = &PulseEstimates::width;
          parameter_name  = "Mean width";
        } else {
          prior_mean_     = &PatientPriors::mass_mean;
          prior_variance_ = &PatientPriors::mass_variance;
          est_mean_       = &PatientEstimates::mass_mean;
          est_sd_         = &PatientEstimates::mass_sd;
          tvarscale_      = &PulseEstimates::tvarscale_mass;
          randomeffect_   = &PulseEstimates::mass;
          parameter_name  = "Mean mass";
        }

      };

    // Exposed publicly so the closed-form log-ratio can be unit tested; the MH
    // base dispatches through its own private virtual, so this does not change
    // normal sampling behavior.
    double posterior_function(Patient *patient, double proposal, bool *notused) {

      double prior_ratio       = 0.0 ;
      double psum_old          = 0.0 ;
      double psum_new          = 0.0 ;
      double newint            = 0.0 ;
      double oldint            = 0.0 ;
      double normalizing_ratio = 0.0 ;
      double prop_ratio        = 0.0 ;
      PatientPriors *priors    = &patient->priors;
      PatientEstimates *est    = &patient->estimates;
      double prior_mass_mean   = (*priors).*prior_mean_;
      double prior_mass_var    = (*priors).*prior_variance_;
      double current           = (*est).*est_mean_;
      double stddev            = (*est).*est_sd_;
      double scale             = 0.0;
      double randomeffect      = 0.0;


      // Prior Ratio -- normal hyperprior on mu itself. Under the log-normal
      // parameterization mu is the mean of the LOG random effect but is still
      // given a normal hyperprior, so this term is identical in both branches.
      prior_ratio = (pow(current - prior_mass_mean, 2) -
                     pow(proposal - prior_mass_mean, 2)) /
                    (2 * prior_mass_var);

      // 'likelihood' ratio -- Ratio of p(theta|mu, sigma, kappa) for current and
      // proposed mu, where the per-pulse t-scale kappa weights each quadratic.
      for (auto &pulse : patient->pulses) {  // uses range based loop instead of iterators

        scale        = pulse.*tvarscale_;
        randomeffect = pulse.*randomeffect_;

        if (lognormal_) {
          // Log-normal parameterization (papers): log(theta) ~ N(mu, sigma^2/kappa).
          // Residuals are formed on the log scale, (log(theta) - mu), and there is
          // no truncation at 0, so the truncated-normal normalizing constant is
          // dropped. mu is not transformed, so no per-proposal Jacobian appears
          // (the log random effects are constants w.r.t. the proposal on mu).
          double logre = log(randomeffect);
          psum_old += pow(logre - current, 2) * scale;
          psum_new += pow(logre - proposal, 2) * scale;
        } else {
          // Natural-scale truncated-normal parameterization (research option).
          psum_old += pow(randomeffect - current, 2) * scale;
          psum_new += pow(randomeffect - proposal, 2) * scale;

          // Normalizing constants (truncation at 0)
          oldint += Rf_pnorm5(current * sqrt(scale) / stddev,
                              0, 1, 1.0, 1.0); // second 1.0 does the log xform for us
          newint += Rf_pnorm5(proposal * sqrt(scale) / stddev,
                              0, 1, 1.0, 1.0); // first 1.0 says to use lower tail
        }

      }

      prop_ratio = 0.5 / pow(stddev, 2) * (psum_old - psum_new);
      // Under log-normal there is no truncation, so oldint/newint stay 0 and the
      // normalizing ratio vanishes.
      normalizing_ratio = oldint - newint;

      return prior_ratio + prop_ratio + normalizing_ratio;

    }

  private:

    double PatientPriors::*prior_mean_;
    double PatientPriors::*prior_variance_;
    double PatientEstimates::*est_mean_;
    double PatientEstimates::*est_sd_;
    double PulseEstimates::*tvarscale_;
    double PulseEstimates::*randomeffect_; //pulse specific mass or width
    bool lognormal_;                       //log-scale pulse random effects?

    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; };

    // Under log-normal, mu is the mean of the LOG random effect and may be any
    // real; under natural-scale truncated-normal it must be positive.
    bool parameter_support(double val, bool *notused) {
      return lognormal_ ? true : (val > 0.0);
    }

};

#endif

