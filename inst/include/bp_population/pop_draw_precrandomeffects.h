#ifndef GUARD_bpmod_population_draw_prec_randomeffects_h
#define GUARD_bpmod_population_draw_prec_randomeffects_h

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

class Pop_DrawPrecRandomEffects :
  public ModifiedMetropolisHastings<Population, Population, double, ProposalVariance>
{

  public:

    // Constructor
    Pop_DrawPrecRandomEffects(double in_pv,
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
          est_prec_         = &PopulationEstimates::width_p2p_prec;
          tvarscale_      = &PulseEstimates::tvarscale_width;
          randomeffect_   = &PulseEstimates::width;
          prec_param_       = &PopulationPriors::width_p2p_prec_param;
          prec_param_rate_ = &PopulationPriors::width_p2p_prec_param_rate;
          parameter_name  = "Prec of pulse widths";
        } else {
          est_mean_       = &PatientEstimates::mass_mean;
          est_prec_         = &PopulationEstimates::mass_p2p_prec;
          tvarscale_      = &PulseEstimates::tvarscale_mass;
          randomeffect_   = &PulseEstimates::mass;
          prec_param_       = &PopulationPriors::mass_p2p_prec_param;
          prec_param_rate_ = &PopulationPriors::mass_p2p_prec_param_rate;
          parameter_name  = "Prec of pulse masses";
        }

      };

  private:

    double PatientEstimates::*est_mean_;
    double PopulationEstimates::*est_prec_;
    double PulseEstimates::*tvarscale_;
    double PulseEstimates::*randomeffect_; //pulse specific mass or width
    
    double PopulationPriors::*prec_param_; //pulse specific mass or width shape in gamma prior
    double PopulationPriors::*prec_param_rate_; //scale in gamma prior for mass or width

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
bool Pop_DrawPrecRandomEffects::parameter_support(double val,  Population *population) {

 // PatientPriors *priors = &patient->priors;
  //double patient_sd_param = (*priors).*sd_param_;

  return (val > 0.0);

}


// posterior_function()
//   Calculates the acceptance ratio for use in modified metropolis hastings
//   sampler (inherited SS_DrawSDRandomEffects::sample() function)
//
double Pop_DrawPrecRandomEffects::posterior_function(Population *population,
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
  double patient_prec      = (*popEst).*est_prec_;
  double patient_prec_param = (*popPriors).*prec_param_;
  double patient_prec_param_rate = (*popPriors).*prec_param_rate_;

  // Calculate pulse-specific portion of acceptance ratio
    for (auto &patient : population->patients) {

      PatientEstimates *est = &patient.estimates;
      double patient_mean   = (*est).*est_mean_;

      for (auto &pulse : patient.pulses) {

          // Normalizing constants for ratio of log likelihoods. They are truncated t-distributions
          stdx_old   = patient_mean * sqrt(patient_prec) * sqrt(pulse.*tvarscale_);
          stdx_new   = patient_mean * sqrt(proposal) * sqrt(pulse.*tvarscale_);
          new_int   += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
          old_int   += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);

          // 3rd part of acceptance ratio: This is for ratio of log likelihoods, which are truncated t-distribuitons
          third_part += pulse.*tvarscale_ *
                        (pulse.*randomeffect_ - patient_mean) *
                        (pulse.*randomeffect_ - patient_mean);

      }

        // 1st and 2nd 'parts' of acceptance ratio: This is for the ratio of log likelihoods, which are truncated t-distn.
        first_part  += patient.get_pulsecount() * (log(sqrt(proposal)) - log(sqrt(patient_prec)));
        second_part += 0.5 * (patient_prec - proposal);
    
        // 4th part of acceptance ratio: Ratio of priors (Gamma distribution)
        fourth_part = (patient_prec_param - 1) * (log (patient_prec) - log(proposal)) + patient_prec_param_rate * (proposal - patient_prec);
    }

  // Compute and return log rho
  return old_int - new_int + first_part + second_part * third_part + fourth_part;

};

#endif

