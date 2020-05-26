#ifndef GUARD_bpmod_population_draw_s2s_sd_h
#define GUARD_bpmod_population_draw_s2s_sd_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>  //check this
#include <bp_datastructures/population.h>


//
// Pop_DrawS2S_SD
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the sd of the individual pulse masses & widths
//   This version has the half-Cauchy prior for the SD of the masses and widths

class Pop_DrawS2S_SD :
  public ModifiedMetropolisHastings<Population, Population, double, ProposalVariance>
{

  public:

    // Constructor
    Pop_DrawS2S_SD(double in_pv,
                           int in_adjust_iter,
                           int in_max_iter,
                           double in_target_ratio,
                           bool for_width,
                           bool verbose,
                           int verbose_iter) :
      ModifiedMetropolisHastings <Ppopulation, Population, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {

        // Choose which set of parameters to use: width or mass
        if (for_width) {
          est_mean_       = &PatientPriors::width_mean;
          est_sd_         = &PatientPriors::width_s2s_sd;
          randomeffect_   = &PatientEstimates::width_mean;
          sd_param_       = &PopulationPriors::width_s2s_sd_param;
          parameter_name = "SD of mean pulse widths";
        } else {
            if (for_mass) {
                est_mean_       = &PatientPriors::mass_mean;
                est_sd_         = &PatientPriors::mass_s2s_sd;
                randomeffect_   = &PatientEstimates::mass_mean;
                sd_param_       = &PopulationPriors::mass_s2s_sd_param;
                parameter_name = "SD of mean pulse masses";
            } else {
                if (for_baseline) {
                    est_mean_       = &PatientPriors::baseline_mean;
                    est_sd_         = &PatientPriors::baseline_sd;
                    randomeffect_   = &PatientEstimates::baseline_halflife(0);
                    sd_param_       = &PopulationPriors::baseline_sd_param;
                    parameter_name = "SD of baselines";
                } else {
                    est_mean_       = &PatientPriors::halflife_mean;
                    est_sd_         = &PatientPriors::halflife_sd;
                    randomeffect_   = &PatientEstimates::baseline_halflife(1);
                    sd_param_       = &PopulationPriors::halflife_sd_param;
                    parameter_name = "SD of halflives";
                }
            }
        }
      };

  private:

    double PatientPriors::*est_mean_;
    double PatientPriors::*est_sd_;
    double PatientEstimates::*randomeffect_; //pulse specific mass or width
    
    double PopulationPriors::*sd_param_; //pulse specific mass or width

    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; };

    bool parameter_support(double val, Patient *patient);  //Need help here: Patient or population here
    double posterior_function(Population *Population, double proposal, Patient *notused);  //Need help here.  Patient at the end population?

};




//------------------------------------------------------------
// Defined functions for SD random effects MMH class
//------------------------------------------------------------

// parameter_support()
//   Defines whether the proposal value is within the parameter support
//   For the Cauchy this is just positive
//    TO DO: can we remove the Patient part of the function?
bool Pop_DrawS2S_SD::parameter_support(double val, Patient *patient) {

 // PatientPriors *priors = &patient->priors;
  //double patient_sd_param = (*priors).*sd_param_;

  return (val > 0.0);

}


// posterior_function()
//   Calculates the acceptance ratio for use in modified metropolis hastings
//   sampler (inherited SS_DrawSDRandomEffects::sample() function)
//
double Pop_DrawS2S_SD::posterior_function(Population *Population,
                                                  double proposal, 
                                                  Patient *notused) {

  double stdx_old    = 0.0;
  double stdx_new    = 0.0;
  double old_int     = 0.0;
  double new_int     = 0.0;
  double first_part  = 0.0;
  double second_part = 0.0;
  double third_part  = 0.0;
  double fourth_part = 0.0;
  PatientPriors *est  = &Population->estimates;  //Assumes there is a population structure with an estimates component, I think.
  PopulationPriors *priors = &Population->priors;
  double pop_mean    = (*est).*est_mean_;
  double pop_sd      = (*est).*est_sd_;
  double pop_sd_param = (*priors).*sd_param_;

  // Calculate pulse-specific portion of acceptance ratio
  for (auto &patient : Population->patients) {

    // Normalizing constants for ratio of log likelihoods. They are truncated t-distributions
    stdx_old   = pop_mean / pop_sd;
    stdx_new   = pop_mean / proposal;
    new_int   += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
    old_int   += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);

    // 3rd part of acceptance ratio: This is for ratio of log likelihoods, which are truncated t-distribuitons
    third_part += (patient.*randomeffect_ - pop_mean) *
                  (patient.*randomeffect_ - pop_mean);

  }

  // 1st and 2nd 'parts' of acceptance ratio: This is for the ratio of log likelihoods, which are truncated t-distn.
  first_part  = population->get_pulsecount() * (log(pop_sd) - log(proposal));
  second_part = 0.5 * ((1 / (pop_sd * pop_sd)) - (1 / (proposal * proposal)));
    
  // 4th part of acceptance ratio: Ratio of priors
    fourth_part = log(pop_sd_param + pop_sd * pop_sd) - log(pop_sd_param + proposal * proposal);

  // Compute and return log rho
  return old_int - new_int + first_part + second_part * third_part + fourth_part;

};

#endif

