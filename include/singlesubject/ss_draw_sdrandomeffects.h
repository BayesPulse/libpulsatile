#ifndef GUARD_ss_draw_sd_randomeffects_h
#define GUARD_ss_draw_sd_randomeffects_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"


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
    SS_DrawSDRandomEffects(double in_pv, int in_adjust_iter, int in_max_iter,
                           double in_target_ratio) :
      ModifiedMetropolisHastings
      <Patient, Patient, double,
       ProposalVariance>::ModifiedMetropolisHastings(in_pv,
                                                     in_adjust_iter,
                                                     in_max_iter,
                                                     in_target_ratio) { };

  private:
    bool parameter_support(double val, Patient *patient);
    double posterior_function(Patient *patient, double proposal, Patient *notused);

};




// 
// Defined functions for SD random effects MMH class
//

// parameter_support()
//   Defines whether the proposal value is within the parameter support
bool SS_DrawSDRandomEffects::parameter_support(double val, Patient *patient) {
  return (val > 0.0 && val < patient->priors->mass_sd_max);
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
  double sd_mass     = patient->estimates->mass_sd;
  double mean_mass   = patient->estimates->mass_mean;
  std::list<PulseEstimate>::const_iterator pulse     = patient->pulses.begin();
  std::list<PulseEstimate>::const_iterator pulse_end = patient->pulses.end();

  // Calculate pulse-specific portion of acceptance ratio
  while (pulse != pulse_end) {

    // Normalizing constants
    stdx_old   = mean_mass / ( sd_mass  / sqrt(pulse->tvarscale_mass) );
    stdx_new   = mean_mass / ( proposal / sqrt(pulse->tvarscale_mass) );
    new_int   += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
    old_int   += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);

    // 3rd part of acceptance ratio
    third_part +=
      pulse->tvarscale_mass * (pulse->mass - mean_mass) * (pulse->mass - mean_mass);

    // Next pulse
    ++pulse;
  }

  // 1st and 2nd 'parts' of acceptance ratio
  first_part  = patient->get_pulsecount() * (log(sd_mass) - log(proposal));
  second_part = 0.5 * ((1 / (sd_mass * sd_mass)) - (1 / (proposal * proposal)));

  // Compute and return log rho
  return old_int - new_int + first_part + second_part * third_part;

};

#endif

