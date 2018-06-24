#ifndef GUARD_bpmod_singlesubject_draw_sd_randomeffects_h
#define GUARD_bpmod_singlesubject_draw_sd_randomeffects_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>


// NOTE: I separated out function definitions here 

//
// SS_DrawFixedEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the mean mass & width
//

class SS_DrawSDRandomEffects : public ModifiedMetropolisHastings<Patient, Patient, double, ProposalVariance>
{

  public:
    // Constructor
    SS_DrawSDRandomEffects(double in_pv, 
                           int in_adjust_iter, 
                           int in_max_iter,
                           double in_target_ratio,
                           bool for_width) :
      ModifiedMetropolisHastings <Patient, Patient, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv,
                                 in_adjust_iter,
                                 in_max_iter,
                                 in_target_ratio) {

         // Choose which set of parameters to use: width or mass
          if (for_width) {
            est_mean_       = &PatientEstimates::width_mean;
            est_sd_         = &PatientEstimates::width_sd;
            tvarscale_      = &PulseEstimates::tvarscale_width;
            randomeffect_   = &PulseEstimates::width;
            sd_max_         = &PatientPriors::width_sd_max;
            // priors->mass_sd_max
          } else {
            est_mean_       = &PatientEstimates::mass_mean;
            est_sd_         = &PatientEstimates::mass_sd;
            tvarscale_      = &PulseEstimates::tvarscale_mass;
            randomeffect_   = &PulseEstimates::mass;
            sd_max_         = &PatientPriors::mass_sd_max;
          }

       };

  private:
    bool parameter_support(double val, Patient *patient);
    double posterior_function(Patient *patient, double proposal, Patient *notused);

    double PatientEstimates::*est_mean_;
    double PatientEstimates::*est_sd_;
    double PulseEstimates::*tvarscale_;
    double PulseEstimates::*randomeffect_; //pulse specific mass or width
    double PatientPriors::*sd_max_; //pulse specific mass or width

};




//------------------------------------------------------------
// Defined functions for SD random effects MMH class
//------------------------------------------------------------

// parameter_support()
//   Defines whether the proposal value is within the parameter support
bool SS_DrawSDRandomEffects::parameter_support(double val, Patient *patient) {

  PatientPriors *priors = patient->priors;
  double patient_sd_max = (*priors).*sd_max_;

  return (val > 0.0 && val < patient_sd_max);

}


// posterior_function()
//   Calculates the acceptance ratio for use in modified metropolis hastings
//   sampler (inherited SS_DrawSDRandomEffects::sample() function)
double SS_DrawSDRandomEffects::posterior_function(Patient *patient, 
                                                  double proposal, 
                                                  Patient *notused) {

  double stdx_old    = 0.0;
  double stdx_new    = 0.0;
  double old_int     = 0.0;
  double new_int     = 0.0;
  double first_part  = 0.0;
  double second_part = 0.0;
  double third_part  = 0.0;
  PatientEstimates *est  = patient->estimates;
  double patient_mean    = (*est).*est_mean_;
  double patient_sd      = (*est).*est_sd_;
  std::list<PulseEstimates>::const_iterator pulse     = patient->pulses.begin();
  std::list<PulseEstimates>::const_iterator pulse_end = patient->pulses.end();

  // Calculate pulse-specific portion of acceptance ratio
  while (pulse != pulse_end) {
    //std::cout << "Hi, I'm in SS_DrawSDRandomEffects::posterior_function()" << std::endl;

    // Normalizing constants
    stdx_old   = patient_mean / ( patient_sd  / sqrt((*pulse).*tvarscale_) );
    stdx_new   = patient_mean / ( proposal / sqrt((*pulse).*tvarscale_) );
    new_int   += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
    old_int   += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);

    // 3rd part of acceptance ratio
    third_part += (*pulse).*tvarscale_ * 
                  ((*pulse).*randomeffect_ - patient_mean) * 
                  ((*pulse).*randomeffect_ - patient_mean);

    // Next pulse
    ++pulse;
  }

  // 1st and 2nd 'parts' of acceptance ratio
  first_part  = patient->get_pulsecount() * (log(patient_sd) - log(proposal));
  second_part = 0.5 * ((1 / (patient_sd * patient_sd)) - (1 / (proposal * proposal)));

  // Compute and return log rho
  return old_int - new_int + first_part + second_part * third_part;

};

#endif

