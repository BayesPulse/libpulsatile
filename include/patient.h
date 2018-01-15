#ifndef GUARD_patient_h
#define GUARD_patient_h

#include <RcppArmadillo.h>
#include <RInside.h>                    // for the embedded R via RInside
#include "datastructures.h"
//#include <list>

//
// patient.h
//   defining the patient class and subclasses
//
// Author: Matt Mulvahill
// Created: 11/14/17
// Notes:
//


//
// Parent class of Patients and Populations, object that MMH and Gibbs
// recognize
//
//struct MCMCSamplingUnit { };


using namespace Rcpp;
//using namespace RcppArmadillo;



//
// Patient struct
//
struct Patient {

  //
  // For single-subject model
  //

  PatientData *data;
  PatientPriors *priors;
  PatientEstimates *estimates;
  std::list<PulseEstimate> pulses;
  std::list<PulseEstimate>::iterator piter = pulses.begin();
  std::list<PulseEstimate> responses;
  std::list<PulseEstimate>::iterator riter = responses.begin();

  // Single-subject constructor
  Patient(PatientData *in_data,
          PatientPriors *in_priors,
          PatientEstimates *in_parms) {

    data = in_data;
    priors = in_priors;
    estimates = in_parms;


    PulseEstimate firstpulse(in_data->fitstart, 1, 1, 1, 1, in_parms->get_decay(), in_data->concentration);
    pulses.push_back(firstpulse);
    //++piter;

  }

  //
  // For population models
  //

  Patient(PatientData *in_data,
          PatientEstimates *in_parms) {

    data = in_data;
    estimates = in_parms;
    // priors // uninitialized

  }


  //
  // Methods for both types/all-models
  //

  // Get current number of pulses
  int get_pulsecount() { return pulses.size(); };


  // likelihood()
  //   computes the current likelihood using the observed log-concentrations and
  //   mean concentration
  double likelihood(bool response_hormone) {
    std::list<PulseEstimate>::const_iterator emptyiter;
    return likelihood(response_hormone, emptyiter);
  }

  double likelihood(bool response_hormone,
                    std::list<PulseEstimate>::const_iterator pulse_excluded) {

    double like = 0;
    arma::vec mean_conc;
    arma::vec *conc;

    if (response_hormone) {
      conc = &data->response_concentration;
    } else {
      conc = &data->concentration;
    }

    // Sum across mean_contribs
    mean_conc = mean_concentration(response_hormone, pulse_excluded);

    arma::vec tmp = (*conc) - mean_conc;
    //std::cout << "conc = " << (*conc) << std::endl;
    //std::cout << "mean_conc = " << mean_conc << std::endl;
    //std::cout << "tmp = " << tmp << std::endl;
    like  = arma::sum(tmp);
    //std::cout << "baseline = " << estimates->baseline << std::endl;
    //std::cout << "like = " << like << std::endl;
    //std::cout << "errorsq = " << estimates->errorsq << std::endl;
    //std::cout << "logerrorsq = " << estimates->get_logerrorsq() << std::endl;
    like  = pow(like, 2);
    //std::cout << "like^2 = " << like << std::endl;
    like /= (-2.0 * estimates->errorsq);
    //std::cout << "like/(-2*errsq) = " << like << std::endl;
    like += -0.5 * conc->n_elem * (1.8378771 + estimates->get_logerrorsq());
    //std::cout << "like + -.5*144*... = " << like << std::endl;

    std::cout << "Final likelihood = " << like << std::endl;

    return like;

  }


  //
  // mean_concentration()
  //   this takes each pulse's mean_contrib vector and sums across them
  //
  arma::vec mean_concentration(bool response_hormone) {

    std::list<PulseEstimate>::const_iterator emptyiter;
    return mean_concentration(response_hormone, emptyiter);

  }

  arma::vec mean_concentration(bool response_hormone,
                               std::list<PulseEstimate>::const_iterator pulse_excluded) {

    arma::vec mean_conc(data->concentration.n_elem);
    mean_conc.fill(0);
    std::list<PulseEstimate>::iterator pulse_iter;
    std::list<PulseEstimate>::const_iterator pulselist_end;

    if (response_hormone) {
      pulse_iter    = responses.begin();
      pulselist_end = responses.end();
    } else {
      pulse_iter    = pulses.begin();
      pulselist_end = pulses.end();
    }

    // Add the contribution to the mean from each pulse
    arma::vec mctrb(data->concentration.n_elem);
    int i = 1;
    while (pulse_iter != pulselist_end) {
      if (pulse_iter != pulse_excluded) {
        mctrb = pulse_iter->get_mean_contribution(data->concentration,
                                                  estimates->get_decay());
        mean_conc += mctrb;
        //std::cout << " Pulse num " << i << " and its mean contrib: " << mctrb << std::endl;
        //std::cout << " Pulse num " << i << " and cumsum mean concentration: " << mean_conc << std::endl;
      }
      ++pulse_iter;
      ++i;
    }

    //std::cout << "Baseline = " << estimates->baseline << std::endl;
    // Add the baseline contribution and log
    mean_conc += estimates->baseline;
    mean_conc = log(mean_conc);

    return mean_conc;

  }




};






////////////////////////////////////////////////////////////
// Chains -- use R data types? YES
////////////////////////////////////////////////////////////
//struct PulseChain { NumericMatrix pulse };
//struct PatientChain { NumericMatrix patient };
//struct PopulationChain { NumericMatrix population };
//struct AssociationChain { NumericMatrix association };

#endif
