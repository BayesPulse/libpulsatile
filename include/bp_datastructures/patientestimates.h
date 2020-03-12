#ifndef GUARD_patientestimates_h
#define GUARD_patientestimates_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/utils.h>

//
// patientestimates.h
//   defining the estimates object used at the patient level
//
// Author: Matt Mulvahill
// Notes:
//

using namespace Rcpp;




//
// PatientEstimates structure
//
//  TODO: abandoned inheritance approach -- clean up comments
//   1. PatientEstimates_Pop is the version of the 'common parms' mcmc chain
//   (patient-level estimates) that is used in the population model.
//
//   2. PatientEstimates_Single is the single-subject version, which adds two
//   parameters. These two parameters are estimated at the population level in
//   the pop model.
//

struct PatientEstimates {

  // Used in all models
  arma::vec baseline_halflife;
  //double baseline;
  //double halflife;
  double errorsq;    // model error (variance)
  double mass_mean;
  double width_mean;
  //int    pulse_count; // function of linked list instead?;
  // Always use these functions to get these values.  removed them as separate
  // member variables to ensure the result is always up-to-date
  double get_decay() { return log(2) / baseline_halflife(1); }
  double get_logerrorsq() { return log(errorsq); }

  // Population model only:
  double mass_p2p_sd_var;
  double mass_s2s_sd_var;
  double width_p2p_sd_var;
  double width_s2s_sd_var;
  double baseline_mean;
  double baseline_s2s_sd_var;
  double halflife_mean;
  double halflife_s2s_sd_var;
  double error_var;

  //
  // For population model
  //
  PatientEstimates(double sv_mass_mean,
                   double sv_mass_p2p_sd_var,
                   double sv_mass_s2s_sd_var,
                   double sv_width_mean,
                   double sv_width_p2p_sd_var,
                   double sv_width_s2s_sd_var,
                   double sv_baseline_mean,
                   double sv_baseline_s2s_sd_var,
                   double sv_halflife_mean,
                   double sv_halflife_s2s_sd_var,
                   double sv_error_var) {

    
    mass_mean           = sv_mass_mean;
    mass_p2p_sd_var     = sv_mass_p2p_sd_var;
    mass_s2s_sd_var     = sv_mass_s2s_sd_var;
    width_mean          = sv_width_mean;
    width_p2p_sd_var    = sv_width_p2p_sd_var;
    width_s2s_sd_var    = sv_width_s2s_sd_var;
    baseline_mean       = sv_baseline_mean;
    baseline_s2s_sd_var = sv_baseline_s2s_sd_var;
    halflife_mean       = sv_halflife_mean;
    halflife_s2s_sd_var = sv_halflife_s2s_sd_var;
    error_var           = sv_error_var;

    //pulse_count = 1;

  }


  //
  // For single-subject model only
  //
  double mass_sd;
  double width_sd;

  // Single-subject constructor:
  PatientEstimates(double sv_baseline,
                   double sv_halflife,
                   double sv_errorsq,
                   double sv_mass_mean,
                   double sv_width_mean,
                   double sv_mass_sd,
                   double sv_width_sd){
    baseline_halflife = { sv_baseline, sv_halflife };
    errorsq     = sv_errorsq;
    mass_mean   = sv_mass_mean;
    width_mean  = sv_width_mean;
    mass_sd  = sv_mass_sd;
    width_sd = sv_width_sd;
  }

};



#endif
