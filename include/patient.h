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

using namespace Rcpp;
using namespace RcppArmadillo;



//
// Patient class
//
class Patient {

  public:

    Patient(NumericMatrix data,
            NumericVector priors,
            NumericVector starting_values) {

      read_data(data);
      PatientPriors priors();
      PatientData data();
      PatientEstimates estimates();

    }

    // Read in data
    // TODO: need to use PatientData constructors instead
    void read_data(NumericMatrix data) {

      // col 0 = time; col 1 = (driver) conc.; col 2 = (response conc.)
      int ncol = data.ncol();
      if (ncol == 2) {
        read_data_in(data(_, 0), data(_, 1));
      } else if (ncol == 3) {
        read_data_in(data(_, 0),  data(_, 1),  data(_, 2));
      } else {
        // **error handling**
      }

    }

    int get_pulsecount() { return estimates.pulse_count };


  private:

    // need access to pulses, pulse_count, priors, data.. basically all of this
    //double likelihood; // calculate as-needed?
    PatientPriors priors;
    PatientData data;
    PatientEstimates estimates;
    AssocEstimates association; // Note: should move to Population class, check how its estimated.
    std::list<PulseEstimate> pulses;


};




////////////////////////////////////////////////////////////
// Chains -- use R data types? YES
////////////////////////////////////////////////////////////
//struct PulseChain { NumericMatrix pulse };
//struct PatientChain { NumericMatrix patient };
//struct PopulationChain { NumericMatrix population };
//struct AssociationChain { NumericMatrix association };

#endif
