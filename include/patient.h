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

  OnlyPatient(PatientData *in_data,
              PatientPriors_Single *in_priors,
              PatientEstimates_Single *in_parms) {

    data = in_data;
    priors = in_priors;
    estimates = in_parms;

  }

  //int get_pulsecount() { return estimates.pulse_count };

  // need access to pulses, pulse_count, priors, data.. basically all of this
  PatientData *data;
  PatientPriors_Single *priors;
  PatientEstimates_Single *estimates;
  //AssocEstimates association; // Note: should move to Population class, check how its estimated.
  std::list<PulseEstimate> pulses;
  std::list<PulseEstimate>::iterator piter = pulses.begin();
  std::list<PulseEstimate> responses;

};


//
// aPatient struct
//  For pop model (OnePatient of many) 
//
struct OnePatient {
//: MCMCSamplingUnit

  OnePatient(PatientData *in_data,
             //PatientPriors_Pop *in_priors,
             PatientEstimates_Pop *in_parms) {

    data = in_data;
    //priors = in_priors;
    estimates = in_parms;

  }

  //int get_pulsecount() { return estimates.pulse_count };

  // need access to pulses, pulse_count, priors, data.. basically all of this
  PatientData *data;
  PatientEstimates_Pop *estimates;
  std::list<PulseEstimate> pulses;
  std::list<PulseEstimate> responses;

  //PopulationEstimates *popest; // Note: should move to population class
  //AssocEstimates association; // Note: should move to Population class, check how its estimated.

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
