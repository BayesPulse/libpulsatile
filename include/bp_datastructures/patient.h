#ifndef GUARD_patient_h
#define GUARD_patient_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include "bp_datastructures/patientdata.h"
#include "bp_datastructures/patientestimates.h"
#include "bp_datastructures/patientpriors.h"
#include "bp_datastructures/pulseestimates.h"



//
// patient.h
//   defining the patient class and subclasses
//
//   NOTE: getting too big -- clean up/reorganize
//


//
// Parent class of Patients and Populations, object that MMH and Gibbs
// recognize
//
//struct MCMCSamplingUnit { };


using namespace Rcpp;

typedef std::list<PulseEstimates> PulseList;
typedef PulseList::iterator PulseIter;
typedef PulseList::const_iterator PulseConstIter;


//
// Patient struct
//
struct Patient {

  // Member objects for all models
  PatientData *data;
  PatientEstimates *estimates;
  PulseList pulses;
  PulseIter piter = pulses.begin();
  PulseList responses;
  PulseIter riter = responses.begin();

  //
  // For single-subject model
  //

  // Member objects
  PatientPriors *priors;

  // Constructor
  Patient(PatientData *in_data,
          PatientPriors *in_priors,
          PatientEstimates *in_parms) {

    data      = in_data;
    priors    = in_priors;
    estimates = in_parms;

    PulseEstimates firstpulse(in_data->fitstart,
                             1, 1, 1, 1,
                             in_parms->get_decay(), 
                             in_data->concentration);
    pulses.push_back(firstpulse);

  }

  //
  // For population models
  //

  // Constructor
  Patient(PatientData *in_data,
          PatientEstimates *in_parms) {

    data = in_data;
    estimates = in_parms;
    // priors member obj not used

  }


  //
  // Methods for both types/all-models
  //

  // get_pulsecount()
  //   Get current number of pulses
  int get_pulsecount() { return pulses.size(); };

  // get_sumerrorsquared()
  //   Sums of squared error for error gibbs
  double get_sumerrorsquared(bool response_hormone) {

    arma::vec mean = mean_concentration(response_hormone);
    arma::vec *conc;

    if (response_hormone) {
      conc = &data->response_concentration;
    } else {
      conc = &data->concentration;
    }

    return arma::accu((*conc - mean) % (*conc - mean)); // % Schur product (elementwise multiplication)

  }

  // likelihood()
  //   computes the current likelihood using the observed log-concentrations and
  //   mean concentration as requested. Not stored.
  //   there is a version for a) using all pulses and b) one for excluding one
  //   pulse.
  double likelihood(bool response_hormone) {
    PulseIter emptyiter;
    return likelihood(response_hormone, emptyiter);
  }

  double likelihood(bool response_hormone, PulseIter pulse_excluded) {

    double like = 0;
    arma::vec *conc;

    if (response_hormone) {
      conc = &data->response_concentration;
    } else {
      conc = &data->concentration;
    }

    // Calculate likelihood
    like  = arma::accu(*conc - mean_concentration(response_hormone, pulse_excluded));
    like  = like * like;
    like /= (-2.0 * estimates->errorsq);
    like += -0.5 * conc->n_elem * (1.8378771 + estimates->get_logerrorsq());

    return like;

  }

  // Calculate all partial likelihoods (excluding each pulse(i))
  arma::vec get_partial_likelihood(bool response_hormone) {

    PulseIter exclude_pulse = pulses.begin();;
    PulseIter pulse_end     = pulses.end();;
    arma::vec partials(get_pulsecount());

    int i = 0;
    while(exclude_pulse != pulse_end) {
      partials(i) = likelihood(response_hormone, exclude_pulse);
      //std::cout << "partials(" << i << ") = " << partials(i) << std::endl;
      exclude_pulse++;
      i++;
    }

    return partials;

  }

  // mean_concentration()
  //   this takes each pulse's mean_contrib vector and sums across them
  //   there is a version for a) using all pulses and b) one for excluding one
  //   pulse.
  arma::vec mean_concentration(bool response_hormone) {

    PulseIter emptyiter;
    return mean_concentration(response_hormone, emptyiter);

  }

  arma::vec mean_concentration(bool response_hormone, PulseIter pulse_excluded)
  {

    arma::vec mean_conc(data->concentration.n_elem);
    mean_conc.fill(0);
    PulseIter pulse_iter;
    PulseConstIter pulselist_end;

    if (response_hormone) {
      pulse_iter    = responses.begin();
      pulselist_end = responses.end();
    } else {
      pulse_iter    = pulses.begin();
      pulselist_end = pulses.end();
    }

    // Add the contribution to the mean from each pulse
    arma::vec mctrb(data->concentration.n_elem);
    int i = 1; // i think this is extraneous
    while (pulse_iter != pulselist_end) {
      if (pulse_iter != pulse_excluded) {
        mctrb = pulse_iter->get_mean_contribution(data->concentration,
                                                  estimates->get_decay());
        mean_conc += mctrb;
      }
      ++pulse_iter;
      ++i;
    }

    // Add the baseline contribution and log
    mean_conc += estimates->baseline_halflife(0);
    mean_conc = log(mean_conc);

    return mean_conc;

  }

  // calc_sr_strauss()
  //   Calculates sum(S(R)), the exponent on the gamma parameter in the Strauss
  //   process/prior for pulse location. Used for Strauss prior in birth_death
  //   and mh_time_strauss.
  //   location is time/loc to test against other pulses
  int calc_sr_strauss(double location) {
    PulseIter emptyiter;
    return calc_sr_strauss(location, &(*emptyiter));
  }

  int calc_sr_strauss(double location, PulseEstimates * pulse_excluded) {

    int s_r = 0;       // Sum of indicators where diff < 20
    double difference; // Time difference
    PulseIter pulse = pulses.begin();
    PulseConstIter pulse_end = pulses.end();

    while (pulse != pulse_end) {
      if (&(*pulse) != pulse_excluded) { // TODO: Test that pulse is actually excluded!
        // skip if node is same that location is from;
        difference = fabs(location - pulse->time);
        // increment by 1 if diff<R
        s_r = (difference < priors->strauss_repulsion_range) ? s_r + 1 : s_r; 
      }
      pulse++;
    }

    // sum(S(R)) - scalar value for sum of # pulses too close to each other
    return(s_r); 

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
