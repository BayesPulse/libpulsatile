#ifndef GUARD_bpmod_singlesubject_draw_sd_randomeffectswidths_h
#define GUARD_bpmod_singlesubject_draw_sd_randomeffectswidths_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>


// NOTE: I separated out function definitions here 

//
// SS_DrawFixedEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the mean width & width
//

class SS_DrawSDRandomEffectsWidths : public ModifiedMetropolisHastings<Patient, Patient, double, ProposalVariance>
{

  public:
    // Constructor
    SS_DrawSDRandomEffectsWidths(double in_pv, int in_adjust_iter, int in_max_iter,
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
bool SS_DrawSDRandomEffectsWidths::parameter_support(double val, Patient *patient) {
  return (val > 0.0 && val < patient->priors->width_sd_max);
}


// posterior_function()
//   Calculates the acceptance ratio for use in modified metropolis hastings
//   sampler (inherited SS_DrawSDRandomEffects::sample() function)
double SS_DrawSDRandomEffectsWidths::posterior_function(Patient *patient, 
                                                  double proposal, 
                                                  Patient *notused) {

  double stdx_old    = 0.0;
  double stdx_new    = 0.0;
  double old_int     = 0.0;
  double new_int     = 0.0;
  double first_part  = 0.0;
  double second_part = 0.0;
  double third_part  = 0.0;
  double sd_width     = patient->estimates->width_sd;
  double mean_width   = patient->estimates->width_mean;
  std::list<PulseEstimates>::const_iterator pulse     = patient->pulses.begin();
  std::list<PulseEstimates>::const_iterator pulse_end = patient->pulses.end();

  // Calculate pulse-specific portion of acceptance ratio
  while (pulse != pulse_end) {

    // Normalizing constants
    stdx_old   = mean_width / ( sd_width  / sqrt(pulse->tvarscale_width) );
    stdx_new   = mean_width / ( proposal / sqrt(pulse->tvarscale_width) );
    new_int   += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
    old_int   += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);

    // 3rd part of acceptance ratio
    third_part +=
      pulse->tvarscale_width * (pulse->width - mean_width) * (pulse->width - mean_width);

    // Next pulse
    ++pulse;
  }

  // 1st and 2nd 'parts' of acceptance ratio
  first_part  = patient->get_pulsecount() * (log(sd_width) - log(proposal));
  second_part = 0.5 * ((1 / (sd_width * sd_width)) - (1 / (proposal * proposal)));

  // Compute and return log rho
  return old_int - new_int + first_part + second_part * third_part;

};

#endif

