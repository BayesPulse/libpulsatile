#ifndef GUARD_chains_h
#define GUARD_chains_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include "bp_datastructures/patient.h"
#include "bp_datastructures/patientdata.h"
#include "bp_datastructures/patientestimates.h"
#include "bp_datastructures/patientpriors.h"
#include "bp_datastructures/pulseestimates.h"



//
// chains.h
//   defining the chains class and subclasses
//
//   NOTE: getting too big -- clean up/reorganize
//


//
// Parent class of chainss and Populations, object that MMH and Gibbs
// recognize
//

// TODO: split out patient chain with constructor bool args for population and
// assocation, then pop creates X number of these + the pop chain

using namespace Rcpp;


struct Chains {


  // 
  // Constructors
  //

  // Base constructor (single subject, single hormone)
  Chains(int iterations, int thin, int burnin, bool response_hormone) {

    int num_outputs = (iterations - burnin) / thin;

    // Create list of length num_outputs for pulse chains
    pulse_chains(num_outputs);

    // Create matrix of nrow num_outputs for patient chain
    // NOTE: How does multiple patient chains fit in here? 
    patient_chain(num_outputs, 9);
    association_chain(num_outputs, 4);

  }

  // Population constructor
  // Chains(int iterations, int thin, int burnin, bool response_hormone, int num_patients) :
  //   Chains(iterations, thin, burnin, response_hormone) {

  //   }


  //
  // Member object definitions
  //
  List pulse_chains;
  NumericMatrix one_set_of_pulses; // objects added to pulse_chains list, ncol=8
  NumericMatrix patient_chain;
  NumericMatrix association_chain;
  NumericMatrix population_chain;

  //
  // Primary user-facing function -- adds 1 iteration to Chains object
  //
  void save_sample(Patient * pat) { };


  //
  // Supporting/Internal functions
  //

  // Function for adding attributes to patient_chain
  NumericMatrix addattribs_patient_chain(NumericMatrix out) {

    out.names() = CharacterVector::create("iteration", "num_pulses", "baseline",
                                          "mean_pulse_mass", "mean_pulse_width",
                                          "halflife", "model_error", "sd_mass",
                                          "sd_widths");
    out.attr("class") = "patient_chain";

    return out;

  }

  // Function for adding attributes to association_chain
  NumericMatrix addattribs_association_chain(NumericMatrix out) {

    // TODO: clarify names -- check Karen's code and AssocEstimates in population.h
    out.names() = CharacterVector::create("iteration", "cluster_size",
                                          "cluster_width");
    out.attr("class") = "association_chain";

    return out;

  }

  // Function for adding attributes to one_set_of_pulses
  NumericMatrix addattribs_set_of_pulses(NumericMatrix out) {

    out.names() = CharacterVector::create("iteration", "total_num_pulses",
                                         "pulse_num", "location",  "mass",
                                         "width", "eta_mass", "eta_width");
    out.attr("class") = "one_set_of_pulses";

    return out;

  }

};

//struct PulseChain { NumericMatrix pulse; };
//struct PatientChain { NumericMatrix patient; };
//struct PopulationChain { NumericMatrix population; };
//struct AssociationChain { NumericMatrix association; };

#endif
