#ifndef GUARD_bp_population_draw_populsemeans_h
#define GUARD_bp_population_draw_populsemeans_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>  //check this.
#include <bp_datastructures/population.h>


//
// SS_DrawFixedEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the mean mass & width
//

class SS_DrawPopPulseMeans :
  public ModifiedMetropolisHastings<Patient, bool, double, ProposalVariance>
{

  public:

    // Constructor
    SS_DrawPopPulseMeans(double in_pv,
                        int in_adjust_iter,
                        int in_max_iter,
                        double in_target_ratio,
                        bool for_width,
                        bool verbose,
                        int verbose_iter) :
      ModifiedMetropolisHastings <Patient, bool, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) { 

        // Choose which set of parameters to use: width or mass
        if (for_width) {
          prior_mean_     = &PopulationPriors::width_mean;
          prior_variance_ = &PopulationPriors::width_variance;
          est_mean_       = &PatientPriors::width_mean;
          est_sd_         = &PatientPriors::width_sd;
          parameter_name  = "Pop Mean width";
        } else {
          prior_mean_     = &PopulationPriors::mass_mean;
          prior_variance_ = &PopulationPriors::mass_variance;
          est_mean_       = &PatientPriors::mass_mean;
          est_sd_         = &PatientPriors::mass_sd;
          parameter_name  = "Pop Mean mass";
        }

      };

  private:

    double PopulationPriors::*prior_mean_;
    double PopulationPriors::*prior_variance_;
    double PatientPriors::*est_mean_;
    double PatientPriors::*est_sd_;
   
    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; };

    bool parameter_support(double val, bool *notused) { return (val > 0.0); }

    double posterior_function(Patient *patient, double proposal, bool *notused) {

      double prior_ratio       = 0.0 ;
      double psum_old          = 0.0 ;
      double psum_new          = 0.0 ;
      double newint            = 0.0 ;
      double oldint            = 0.0 ;
      double normalizing_ratio = 0.0 ;
      double prop_ratio        = 0.0 ;
      PopulationPriors *priors    = &patient->priors;  //need help here?  This isn't in patient.
      PatientPriors *est    = &patient->priors;
      double prior_mass_mean   = (*priors).*prior_mean_;  //why mass, can we remove
      double prior_mass_var    = (*priors).*prior_variance_;  //why mass, can we remove.
      double current           = (*est).*est_mean_;
      double stddev            = (*est).*est_sd_;


      // Prior Ratio: This is ratio of truncated Normals
      prior_ratio = (pow(current - prior_mass_mean, 2) -
                     pow(proposal - prior_mass_mean, 2)) /
                    (2 * prior_mass_var);

      // 'likelihood' ratio -- Ratio of p(mu_alpha|m, v) for current and
      // proposed alpha
      for (auto &patient : patient->estimates) {  // uses range based loop instead of iterators

/**STOPPED HERE*/
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

