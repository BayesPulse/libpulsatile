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
// OnlyPatient struct
//   For single-subject model
//
struct OnlyPatient {
//: MCMCSamplingUnit

  // need access to pulses, pulse_count, priors, data.. basically all of this
  PatientData *data;
  PatientPriors_Single *priors;
  PatientEstimates_Single *estimates;
  std::list<PulseEstimate> pulses;
  std::list<PulseEstimate>::iterator piter = pulses.begin();
  std::list<PulseEstimate> responses;
  std::list<PulseEstimate>::iterator riter = responses.begin();

  // Constructor
  OnlyPatient(PatientData *in_data,
              PatientPriors_Single *in_priors,
              PatientEstimates_Single *in_parms) {

    data = in_data;
    priors = in_priors;
    estimates = in_parms;

  }

  // Get current number of pulses
  int get_pulsecount() { return pulses.size(); };


  // likelihood()
  //   computes the current likelihood using the observed log-concentrations and
  //   mean concentration
  double likelihood(Patient *patient,
                    Pulse *pulse_excluded,
                    bool response_hormone = FALSE){

    int i;
    int N = get_pulsecount();  // function not yet defined **
    double like = 0; // likelihood to be calculated
    arma::vec mean_conc;
    arma::vec *data;

    if (response_hormone) {
      data = patient->data.response_concentration;
    } else {
      data = patient->data.concentration;
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
  arma::vec mean_concentration(Pulse **patient, // list of pulses
                               Pulse *pulse_excluded,
                               int pulse_count)
  {

    int i;
    Pulse *pulse;
    int N = pulse_count; //patient.get_pulsecount();  // function not yet defined **
    arma::vec mean_conc(N);

    // Add the contribution to the mean from each pulse
    pulse = Patient->pulses->next; // move to first pulse, not sure if works**
    while (pulse != NULL) {
      if (pulse != pulse_excluded) {
        for (i = 0; i < N; i++)
          mean_conc(i) += pulse->mean_contrib(i); // mean contrib not yet fully defined **
      }
      node = pulse->next;
    }

    // Add the baseline contribution and log
    mean_conc += baseline;
    mean_conc = log(mean_conc);

    return mean_conc;
  }




};




//
// OnePatient struct
//  For pop model (OnePatient of many) 
//
struct OnePatient {
//: MCMCSamplingUnit

  OnePatient(PatientData *in_data,
             PatientEstimates_Pop *in_parms) {

    data = in_data;
    estimates = in_parms;

  }

  int get_pulsecount() { return pulses.size(); };

  PatientData *data;
  PatientEstimates_Pop *estimates;
  std::list<PulseEstimate> pulses;
  std::list<PulseEstimate>::iterator piter = pulses.begin();
  std::list<PulseEstimate> responses;
  std::list<PulseEstimate>::iterator riter = responses.begin();

};



// Example handler class
//class Student_info { 
//
//  public:
//    // constructors and copy control
//    Student_info(): cp(0) { } 
//    Student_info(std::istream& is): cp(0) { read(is); } 
//    Student_info(const Student_info&);
//    Student_info& operator=(const Student_info&); 
//    ~Student_info() { delete cp; }
//
//    // operations
//    std::istream& read(std::istream&);
//    std::string name() const {
//      if (cp) return cp->name();
//      else throw std::runtime_error("uninitialized Student");
//    }
//    double grade() const {
//      if (cp) return cp->grade();
//      else throw std::runtime_error("uninitialized Student");
//    }
//
//    static bool compare(const Student_info& s1,
//                        const Student_info& s2) { return s1.name() < s2.name(); }
//
//  private:
//    Core* cp;
//
//}


////////////////////////////////////////////////////////////
// Chains -- use R data types? YES
////////////////////////////////////////////////////////////
//struct PulseChain { NumericMatrix pulse };
//struct PatientChain { NumericMatrix patient };
//struct PopulationChain { NumericMatrix population };
//struct AssociationChain { NumericMatrix association };

#endif
