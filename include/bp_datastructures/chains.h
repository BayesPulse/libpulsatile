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

//struct PulseChain { NumericMatrix pulse; };
//struct PatientChain { NumericMatrix patient; };
//struct PopulationChain { NumericMatrix population; };
//struct AssociationChain { NumericMatrix association; };


//
// Parent class of chains and Populations, object that MMH and Gibbs
// recognize
//

// TODO: split out patient chain with constructor bool args for population and
// assocation, then pop creates X number of these + the pop chain

using namespace Rcpp;





class Chains {

  public: 

    //
    // Constructors
    //

    // Base constructor (single subject, single hormone)
    Chains(int in_iterations, int in_thin, int in_burnin, bool in_response_hormone) :
      num_outputs((in_iterations - in_burnin) / in_thin),
      patient_chain(num_outputs, 9, arma::fill::zeros) {

      //patient_chain.fill(0.0);
      iterations       = in_iterations;
      thin             = in_thin;
      burnin           = in_burnin;
      response_hormone = in_response_hormone;

      //std::fill(patient_chain.begin(), patient_chain.end(), 0.0);

      // Model type
      model_type = "single_subject";

    }

    // Population constructor
    // Chains(int iterations, int thin, int burnin, bool response_hormone, int num_patients) :
    //   Chains(iterations, thin, burnin, response_hormone) {
    //   association_chain(num_outputs, 4);
    //  model_type = "population";
    //   }


    //
    // Member scalar definitions
    //
    String model_type;
    int iterations;
    int thin;
    int burnin;
    bool response_hormone;
    int num_outputs;

    //
    // Member object definitions
    //
    std::vector<arma::mat> pulse_chains;
    arma::mat one_set_of_pulses; // objects added to pulse_chains list, ncol=8
    arma::mat patient_chain;
    arma::mat association_chain;
    arma::mat population_chain;

  //
  // Primary user-facing functions -- adds 1 iteration to Chains object
  //

  // Record this iteration 
  void save_sample(Patient * pat, int iter);
  //void save_sample(Population * pop);

  // Return chains function
  List output();



  //
  // Supporting/Internal functions
  //

  // Functions for adding attributes to chains
  NumericMatrix addattribs_patient_chain(NumericMatrix out);
  //NumericMatrix addattribs_association_chain(NumericMatrix out);
  NumericMatrix addattribs_set_of_pulses(NumericMatrix out);
  //// Function for adding attributes to population_chain
  //// TODO: Straighten out variance vs SD terms and why is there a variance AND
  //// SD term for mass/width
  //NumericMatrix addattribs_population_chain(NumericMatrix out);

};








// Record this iteration 
void Chains::save_sample(Patient * pat, int iter) {

  if (iter > burnin & (iter % thin) == 0) {

    // Calculations used repeatedly in this function
    int output_num = (iter - burnin) / thin;
    double pulsecount = (double)pat->get_pulsecount();

    // Fill patient chain w/ current patient-level estimates
    patient_chain(output_num, 0) = (double)iter;
    patient_chain(output_num, 1) = pulsecount;
    patient_chain(output_num, 2) = pat->estimates->baseline_halflife(0);
    patient_chain(output_num, 3) = pat->estimates->mass_mean;
    patient_chain(output_num, 4) = pat->estimates->width_mean;
    patient_chain(output_num, 5) = pat->estimates->baseline_halflife(1);
    patient_chain(output_num, 6) = pat->estimates->errorsq;
    patient_chain(output_num, 7) = pat->estimates->mass_sd;
    patient_chain(output_num, 8) = pat->estimates->width_sd;

    // Create a matrix of current pulse-level estimates and add matrix to the
    //   vector chain

    // Construct matrix of parameters that need constructiong (current iter,
    // number of pulses, pulse number id) and add to a matrix
    arma::vec itervec(pulsecount); itervec.fill((double)iter);
    arma::vec pcvec(pulsecount); pcvec.fill(pulsecount);
    arma::mat constructed_parms(pulsecount, 3, arma::fill::zeros);
    constructed_parms.col(0) = itervec;
    constructed_parms.col(1) = pcvec;
    constructed_parms.col(2) = arma::linspace<arma::vec>(1, pulsecount, pulsecount);

    // Create matrix of current pulse-level estimates
    arma::mat pulseparms(pulsecount, 6, arma::fill::zeros);
    int i = 0;
    for (auto &pulse : pat->pulses) {
      pulseparms.row(i) = pulse.get_vector_of_values();
      i++;
    }

    // Combine constructed (ids) and pulse-level estimates and add to vector of
    // matrices
    arma::mat new_mat = arma::join_rows(constructed_parms, pulseparms);
    pulse_chains.push_back(new_mat);

  }

};

//// Return chains function
//List Chains::output() {
//
//  // Add names to each output chain
//  NumericMatrix patient_chain_r = as<NumericMatrix>(wrap(patient_chain));
//  patient_chain_r = addattribs_patient_chain(patient_chain_r);
//  // would like to do this, but requires patient_chain class:
//  // patient_chain.add_attribs(); 
//
//  //pulse_chains = addattribs_pulse_chain(pulse_chains);
//
//  // Create list object combining all chains & other output
//  List out = List::create(Named("patient_chain") = patient_chain);//,
//                        //Named("pulse_chains") = pulse_chains);
//
//  return out;
//
//}
//
////void save_sample(Population * pop) {
//
////};
//
//
//
// Supporting/Internal functions
//

// Function for adding attributes to patient_chain
NumericMatrix Chains::addattribs_patient_chain(NumericMatrix out) {

  colnames(out) = CharacterVector::create("iteration",
                                          "num_pulses",
                                          "baseline",
                                          "mass_mean",
                                          "width_mean",
                                          "halflife",
                                          "model_error",
                                          "mass_sd",
                                          "width_sd");
  out.attr("class") = "patient_chain";

  return out;

}

//// Function for adding attributes to association_chain
//NumericMatrix Chains::addattribs_association_chain(NumericMatrix out) {
//
//  // TODO: clarify names -- check Karen's code and AssocEstimates in population.h
//  colnames(out) = CharacterVector::create("iteration", 
//                                          "cluster_size",
//                                          "cluster_width");
//  out.attr("class") = "association_chain";
//
//  return out;
//
//}

// Function for adding attributes to one_set_of_pulses
NumericMatrix Chains::addattribs_set_of_pulses(NumericMatrix out) {

  colnames(out) = CharacterVector::create("iteration",
                                          "total_num_pulses",
                                          "pulse_num",
                                          "location",
                                          "mass",
                                          "width",
                                          "eta_mass",
                                          "eta_width");
  out.attr("class") = "one_set_of_pulses";

  return out;

}

// Function for adding attributes to population_chain
// TODO: Straighten out variance vs SD terms and why is there a variance AND
// SD term for mass/width
//NumericMatrix Chains::addattribs_population_chain(NumericMatrix out) {
//  colnames(out) = CharacterVector::create("iteration",
//                                          "baseline_mean",
//                                          "baseline_variance",
//                                          "halflife_mean",
//                                          "halflife_variance",
//                                          "mass_mean",
//                                          "mass_variance",
//                                          "mass_mean_sd",
//                                          "width_mean",
//                                          "width_variance",
//                                          "width_mean_sd");
//  out.attr("class") = "population_chain";
//
//  return out;
//
//}

#endif
