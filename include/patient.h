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
  double likelihood(bool response_hormone,
                    //Patient *patient,
                    Pulse *pulse_excluded) {

    int i;
    int N = get_pulsecount();  // function not yet defined **
    double like = 0; // likelihood to be calculated
    arma::vec mean_conc;
    arma::vec *data;

    if (response_hormone) {
      data = data.response_concentration;
    } else {
      data = data.concentration;
    }

    // Sum across mean_contribs
    mean_conc = mean_concentration(patient, pulse_excluded);
    for (i = 0; i < N; i++) {
      like += pow(data(i) - mean_conc(i), 2);  // should be able to get rid of the loop by element diffing these vectors and squaring
    }

    like /= (-2.0 * patient->estimates.errorsq);
    like += -0.5 * N * (1.8378771 + patient->estimates.logerrorsq);

    return like;

  }

  // mean_concentration()
  //   this takes each pulse's mean_contrib vector and sums across them
  arma::vec mean_concentration(bool response_pulses) {

    std::list<PulseEstimate>::const_iterator emptyiter;
    return mean_concentration(response_pulses, emptyiter);

  }
  arma::vec mean_concentration(bool response_pulses,
                               std::list<PulseEstimate>::const_iterator pulse_excluded) {

    arma::vec mean_conc(data->concentration.n_elem);
    std::list<PulseEstimate>::iterator pulse_iter;
    std::list<PulseEstimate>::const_iterator pulselist_end;

    if (response_pulses) {
      pulse_iter    = responses.begin();
      pulselist_end = responses.end();
    } else {
      pulse_iter    = pulses.begin();
      pulselist_end = pulses.end();
    }

    // Add the contribution to the mean from each pulse
    //++pulse_iter; // move to first pulse, not sure if works**
    int i = 0; // for testing only
    std::cout << "Patient's halflife = " << estimates->halflife << " and decay rate = " << estimates->get_decay() << std::endl;
    while (pulse_iter != pulselist_end) {
      i++; 
      if (pulse_iter != pulse_excluded) {
        std::cout << "Pulse number " << i << "'s mean concentration" << std::endl;
        mean_conc += pulse_iter->get_mean_contribution(data->concentration, estimates->get_decay());
      }
      ++pulse_iter;
    }

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
