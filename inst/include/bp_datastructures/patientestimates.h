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
  double baseline;
  double halflife;
  double errorsq;    // model error (variance)
  double mass_mean;
  double width_mean;
  //int    pulse_count; // function of linked list instead?;
  // Always use these functions to get these values.  removed them as separate
  // member variables to ensure the result is always up-to-date
  double get_decay() { 
    //if (baseline_halflife.n_elem == 2) {
      return log(2) / baseline_halflife(1);
    //} else {
    //  return log(2) / halflife;
    //}
  }
  double get_logerrorsq() { return log(errorsq); }

  
  //
  // For population model
  //
  PatientEstimates(double sv_baseline,
                   double sv_halflife,
                   double sv_error_var,
                   double sv_mass_mean,
                   double sv_width_mean,
                   double sv_mass_prec,
                   double sv_width_prec,
                   bool notused){
    baseline_halflife = {sv_baseline, sv_halflife};
    baseline = sv_baseline;
    halflife = sv_halflife;
    errorsq = sv_error_var;
    mass_mean   = sv_mass_mean;
    width_mean  = sv_width_mean;
    mass_prec = sv_mass_prec;
    width_prec = sv_width_prec;
    
  }

  // Coordinates vector and double BL/HL
  //   Vector used for subject level draw
  //   doubles used for population draws
  void matchBLHL() {
    baseline = baseline_halflife(0);
    halflife = baseline_halflife(1);
  }

  //
  // For single-subject model only
  // We define the precision of the pulse-to-pulse variation in this version
  //
  double mass_prec;
  double width_prec;

  // Single-subject constructor:
  PatientEstimates(double sv_baseline,
                   double sv_halflife,
                   double sv_errorsq,
                   double sv_mass_mean,
                   double sv_width_mean,
                   double sv_mass_prec,
                   double sv_width_prec){
    baseline_halflife = { sv_baseline, sv_halflife };
    errorsq     = sv_errorsq;
    mass_mean   = sv_mass_mean;
    width_mean  = sv_width_mean;
    mass_prec = sv_mass_prec;
    width_prec = sv_width_prec;
  }

};



#endif
