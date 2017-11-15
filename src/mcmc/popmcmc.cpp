#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "mh.h"
#include "proposalvariance.h"
#include "proposalvariance.h"
#include "proposalvariance.h"
#include "proposalvariance.h"

using namespace Rcpp;

// popmcmc.cpp
//   MCMC loop for population model
//
// Author: Matt Mulvahill
// Created: 11/14/17
// Notes:
//

void run_mcmc_pop(PatientList patients, int iters) {


  //
  // Create mcmc objects
  //

  //   mass-related parameters
  drawPopMassMean   pop_mu_mean;
  drawPopMassSD     pop_mu_sd;
  drawSubjMassMean  subj_mu_mean;
  drawSubjMassSD    subj_mu_sd;
  drawPulseMassMean pulse_mu;

  //   width-related parameters
  drawPopWidthMean   pop_omega_mean;
  drawPopWidthSD     pop_omega_sd;
  drawSubjWidthMean  subj_omega_mean;
  drawSubjWidthSD    subj_omega_sd;
  drawPulseWidthMean pulse_omega;

  //   baseline and halflife related parameters
  drawPopBhMean pop_theta_mean;
  drawPopBhSd   pop_theta_sd;
  drawSubjBh    subj_theta;

  //   pulse location related parameters 
  drawPulseTimes pulse_tau;

  //   Birth-death object
  drawBirthDeath bd(&patient, pulse_mu, pulse_omega, pulse_tau);
                    

  for (i = 0; i < iters; i++) {
    for (j = 0; j < num_patients; j++) {

      bd.sample              (&patients[j]) ;
      pop_mu_mean.sample     (&patients[j]) ;
      pop_mu_sd.sample       (&patients[j]) ;
      subj_mu_mean.sample    (&patients[j]) ;
      subj_mu_sd.sample      (&patients[j]) ;
      pulse_mu.sample        (&patients[j]) ;
      pop_omega_mean.sample  (&patients[j]) ;
      pop_omega_sd.sample    (&patients[j]) ;
      subj_omega_mean.sample (&patients[j]) ;
      subj_omega_sd.sample   (&patients[j]) ;
      pulse_omega.sample     (&patients[j]) ;
      pop_theta_mean.sample  (&patients[j]) ;
      pop_theta_sd.sample    (&patients[j]) ;
      subj_theta.sample      (&patients[j]) ;
      pulse_tau.sample       (&patients[j]) ;
      error.sample           (&patients[j]) ;

    }

    print_object.print();

  }

  // Patient is updated internally, so return nothing

}
