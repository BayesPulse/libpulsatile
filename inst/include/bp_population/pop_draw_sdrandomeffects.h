#ifndef GUARD_bpmod_population_draw_sd_randomeffects_h
#define GUARD_bpmod_population_draw_sd_randomeffects_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>
#include <bp_datastructures/population.h>


// NOTE: I separated out function definitions here 

//
// Pop_DrawSDRandomEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the sd of the pulse masses & widths
//   This version has the half-Cauchy prior for the SD of the masses and widths

class Pop_DrawSDRandomEffects :
  public ModifiedMetropolisHastings<Population, Population, double, ProposalVariance>
{

  public:

    // Constructor
    Pop_DrawSDRandomEffects(double in_pv,
                           int in_adjust_iter,
                           int in_max_iter,
                           double in_target_ratio,
                           bool for_width,
                           bool verbose,
                           int verbose_iter) :
      ModifiedMetropolisHastings <Population, Population, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {

        // Choose which set of parameters to use: width or mass
        if (for_width) {
          est_mean_       = &PatientEstimates::width_mean;
          est_sd_         = &PopulationEstimates::width_p2p_sd;
          tvarscale_      = &PulseEstimates::tvarscale_width;
          randomeffect_   = &PulseEstimates::width;
          sd_param_       = &PopulationPriors::width_p2p_sd_param;
          sd_param_rate_ = &PopulationPriors::width_p2p_sd_rate;
          parameter_name  = "SD of pulse widths";
        } else {
          est_mean_       = &PatientEstimates::mass_mean;
          est_sd_         = &PopulationEstimates::mass_p2p_sd;
          tvarscale_      = &PulseEstimates::tvarscale_mass;
          randomeffect_   = &PulseEstimates::mass;
          sd_param_       = &PopulationPriors::mass_p2p_sd_param;
          sd_param_rate_ = &PopulationPriors::mass_p2p_sd_rate;
          parameter_name  = "SD of pulse masses";
        }

      };

  private:

    double PatientEstimates::*est_mean_;
    double PopulationEstimates::*est_sd_;
    double PulseEstimates::*tvarscale_;
    double PulseEstimates::*randomeffect_; //pulse specific mass or width
    
    double PopulationPriors::*sd_param_; //pulse specific mass or width shape in gamma prior
    double PopulationPriors::*sd_param_rate_; //scale in gamma prior for mass or width

    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; };

    bool parameter_support(double val, Population *population);
    double posterior_function(Population *population, double proposal, Population *notused);
       //What should be plugged into posterior function, population vs. patient
};




//------------------------------------------------------------
// Defined functions for SD random effects MMH class
//------------------------------------------------------------

// parameter_support()
//   Defines whether the proposal value is within the parameter support
//   For the Cauchy this is just positive
//    TO DO: can we remove the Patient part of the function?
bool Pop_DrawSDRandomEffects::parameter_support(double val,  Population *population) {

 // PatientPriors *priors = &patient->priors;
  //double patient_sd_param = (*priors).*sd_param_;

  return (val > 0.0);

}


// posterior_function()
//   Calculates the acceptance ratio for use in modified metropolis hastings
//   sampler (inherited SS_DrawSDRandomEffects::sample() function)
//
double Pop_DrawSDRandomEffects::posterior_function(Population *population,
                                                   double proposal, 
                                                   Population *notused) {

  // Need to loop through all the patients.
  //Not sure how to do the temp allocations to loop.
    
  double stdx_old    = 0.0;
  double stdx_new    = 0.0;
  double old_int     = 0.0;
  double new_int     = 0.0;
  double first_part  = 0.0;
  double second_part = 0.0;
  double third_part  = 0.0;
  double fourth_part = 0.0;
  //PatientEstimates *est  = &patient->estimates;
  PopulationEstimates *popEst = &population->estimates;
  PopulationPriors *popPriors = &population->priors;
  //double patient_mean    = (*est).*est_mean_;
  double patient_sd      = (*popEst).*est_sd_;
  double patient_sd_param = (*popPriors).*sd_param_;
  double patient_sd_param_rate = (*popPriors).*sd_param_rate_;

  // Calculate pulse-specific portion of acceptance ratio
    for (auto &patient : population->patients) {

      PatientEstimates *est = &patient.estimates;
      double patient_mean   = (*est).*est_mean_;

      for (auto &pulse : patient.pulses) {

          // Normalizing constants for ratio of log likelihoods. They are truncated t-distributions
          stdx_old   = patient_mean / ( patient_sd  / sqrt(pulse.*tvarscale_) );
          stdx_new   = patient_mean / ( proposal / sqrt(pulse.*tvarscale_) );
          new_int   += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
          old_int   += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);

          // 3rd part of acceptance ratio: This is for ratio of log likelihoods, which are truncated t-distribuitons
          third_part += pulse.*tvarscale_ *
                        (pulse.*randomeffect_ - patient_mean) *
                        (pulse.*randomeffect_ - patient_mean);

      }

        // 1st and 2nd 'parts' of acceptance ratio: This is for the ratio of log likelihoods, which are truncated t-distn.
        first_part  += patient.get_pulsecount() * (log(patient_sd) - log(proposal));
        second_part += 0.5 * ((1 / (patient_sd * patient_sd)) - (1 / (proposal * proposal)));
    
        // 4th part of acceptance ratio: Ratio of priors (Gamma distribution)
        fourth_part = (patient_sd_param - 1) * (log(proposal * proposal) - log (patient_sd * patient_sd)) + patient_sd_param_rate * (1/(proposal * proposal) - 1/(patient_sd * patient_sd));
    }

  // Compute and return log rho
  return old_int - new_int + first_part + second_part * third_part + fourth_part;

};

#endif

