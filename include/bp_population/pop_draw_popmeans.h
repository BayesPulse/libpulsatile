#ifndef GUARD_bp_population_draw_popmeans_h
#define GUARD_bp_population_draw_popmeans_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>  //check this.
#include <bp_datastructures/population.h>


//
// Pop_DrawPopMeans
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the population mean of the patient mean mass, width, baseline & halflife
//

class Pop_DrawPopMeans :
  public ModifiedMetropolisHastings<Population, bool, double, ProposalVariance>
  {

  public:

    // Constructor
    Pop_DrawPopMeans(double in_pv,
                     int in_adjust_iter,
                     int in_max_iter,
                     double in_target_ratio,
                     bool for_width,
                     bool for_mass,
                     bool for_baseline,
                     bool verbose,
                     int verbose_iter) :
      ModifiedMetropolisHastings <Population, bool, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) { 

        // Choose which set of parameters to use: width, mass, baseline, or halflife
        if (for_width) {
          prior_mean_     = &PopulationPriors::width_mean;
          prior_variance_ = &PopulationPriors::width_variance;
          est_mean_       = &PatientPriors::width_mean;
          est_sd_         = &PatientPriors::width_s2s_sd;
          randomeffect_   = &PatientEstimates::width_mean;
          parameter_name  = "Pop Mean width";
        } else if (for_mass) {
          prior_mean_     = &PopulationPriors::mass_mean;
          prior_variance_ = &PopulationPriors::mass_variance;
          est_mean_       = &PatientPriors::mass_mean;
          est_sd_         = &PatientPriors::mass_s2s_sd;
          randomeffect_   = &PatientEstimates::mass_mean;
          parameter_name  = "Pop Mean mass";
        } else if (for_baseline) {
          prior_mean_     = &PopulationPriors::baseline_mean;
          prior_variance_ = &PopulationPriors::baseline_variance;
          est_mean_       = &PatientPriors::baseline_mean;
          est_sd_         = &PatientPriors::baseline_sd;
          randomeffect_    = &PatientEstimates::baseline;
          //randomeffect_   = &PatientEstimates::baseline_halflife.colptr(0);  //need help to check that this pulls the first value of this vector
          parameter_name  = "Pop Mean baseline";
        } else {
          prior_mean_     = &PopulationPriors::halflife_mean;
          prior_variance_ = &PopulationPriors::halflife_variance;
          est_mean_       = &PatientPriors::halflife_mean;
          est_sd_         = &PatientPriors::halflife_sd;
          randomeffect_    = &PatientEstimates::halflife;
          //randomeffect_   = &PatientEstimates::baseline_halflife.memptr(1);  //need help to check that this pulls the 2nd value of this vector
          parameter_name  = "Pop Mean halflife";
        }
      };

  private:

    double PopulationPriors::*prior_mean_;
    double PopulationPriors::*prior_variance_;
    double PatientPriors::*est_mean_;
    double PatientPriors::*est_sd_;
    double PatientEstimates::*randomeffect_; //patient specific mean mass, mean width, baseline, or halflife
   
    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; };

    bool parameter_support(double val, bool *notused) { return (val > 0.0); }

    double posterior_function(Population *population, double proposal, bool *notused) {

      double prior_ratio          = 0.0 ;
      double psum_old             = 0.0 ;
      double psum_new             = 0.0 ;
      double newint               = 0.0 ;
      double oldint               = 0.0 ;
      double normalizing_ratio    = 0.0 ;
      double prop_ratio           = 0.0 ;
      PopulationPriors *priors    = &population->popPriors;  //need help here?  This isn't in patient so I put in population which would have all patients and fixed prior info.
      PatientPriors *est          = &population->patPriors;
      double prior_mean           = (*priors).*prior_mean_;  //N: why mass, can we remove, M: removed
      double prior_var            = (*priors).*prior_variance_;  
      double current              = (*est).*est_mean_;
      double stddev               = (*est).*est_sd_;
      double randomeffect         = 0.0;
      
      //Rcpp::Rcout << "Prior mean: " << prior_mean << " Prior var: " << prior_var
      //            << "\nEst mean: " << current << " Est sd: " << stddev << "\n";


      // Prior Ratio: This is ratio of truncated Normals
      prior_ratio = (pow(current - prior_mean, 2) -
                     pow(proposal - prior_mean, 2)) /
                    (2 * prior_var);
      
      //Rcpp::Rcout << "Prior ratio: " << prior_ratio << "\n";

      // 'likelihood' ratio -- Ratio of p(mu_alpha|m, v) for current and
      // proposed alpha
      for (auto &patient : population->patients) {  // uses range based loop instead of iterators. Max: Need check here. Want to loop through patients in a population structure containing all patients

        randomeffect = patient.estimates.*randomeffect_;
        psum_old += pow(randomeffect - current, 2);
        psum_new += pow(randomeffect - proposal, 2);

        // Normalizing constants
        oldint += Rf_pnorm5(current / stddev,
                            0, 1, 1.0, 1.0); // second 1.0 does the log xform for us
        newint += Rf_pnorm5(proposal / stddev,
                            0, 1, 1.0, 1.0); // first 1.0 says to use lower tail

      }

      prop_ratio = 0.5 / pow(stddev, 2) * (psum_old - psum_new);
      normalizing_ratio = oldint - newint;
      
      //Rcpp::Rcout << "Prop ratio: " << prop_ratio << " Norm ratio: " << normalizing_ratio
      //            << "\n";

      return prior_ratio + prop_ratio + normalizing_ratio;

    }

};

#endif

