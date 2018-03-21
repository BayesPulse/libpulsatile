#ifndef GUARD_chains_h
#define GUARD_chains_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
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


using namespace Rcpp;


struct Chains {


  Chains(int iterations, int thin, int burnin) {

    int num_outputs = (iterations - burnin) / thin;
    // Create list of length num_outputs for pulse chains
    pulse_chains(num_outputs);
    // Create matrix of nrow num_outputs for patient chain
    // NOTE: How does multiple patient chains fit in here? 
    patient_chain(num_outputs, 9);

  }

  List pulse_chains;
  NumericMatrix patient_chain;

  void save_sample(Patient * pat) { };


};

struct PulseChain { NumericMatrix pulse };
struct PatientChain { NumericMatrix patient };
struct PopulationChain { NumericMatrix population };
struct AssociationChain { NumericMatrix association };

#endif
