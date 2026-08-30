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

typedef std::vector<arma::mat>  MatrixVector;

using namespace Rcpp;



//
// chains.h
//   defining the chains class and subclasses
//

//struct PulseChain { NumericMatrix pulse; };
//struct PatientChain { NumericMatrix patient; };


//
// Parent class of chains and Populations, object that MMH and Gibbs
// recognize
//

// TODO: split out patient chain with constructor bool args for population and
// association, then pop creates X number of these + the pop chain

// Chains class
//   currently only single-subject (pulse and common)
//
// Description:
//   Key components:
//     - constructor (response y/n)
//     - add R attributes
//     - save_sample() function
//     - output() function
class Chains {

  public:

    // Base constructor (single subject, single hormone)
    Chains(int in_iterations, int in_thin, int in_burnin, bool
           in_response_hormone, bool in_verbose, int in_verbose_iter)
      : num_outputs((in_iterations - in_burnin) / in_thin)
      , patient_chain(num_outputs, 9, arma::fill::zeros) {

        //patient_chain.fill(0.0);
        iterations       = in_iterations;
        thin             = in_thin;
        burnin           = in_burnin;
        response_hormone = in_response_hormone;
        verbose          = in_verbose;
        verbose_iter     = in_verbose_iter;

      }

    //----------------------------------------
    // Member definitions
    //----------------------------------------
    // Member scalar definitions
    int iterations;
    int thin;
    int burnin;
    bool response_hormone;
    int num_outputs;
    bool verbose; 
    int verbose_iter;
    MatrixVector pulse_chains;
    arma::mat patient_chain;

    //----------------------------------------
    // Primary user-facing functions
    //----------------------------------------
    // Record this iteration -- adds 1 iteration to Chains object
    void save_sample(Patient * pat, int iter);
    // Print pulse/patient diagnostic info
    void print_diagnostic_output(Patient * patient, int iter);
    // Return chains function
    List output(Patient * pat);

  private: 
    //----------------------------------------
    // Supporting/Internal functions
    //----------------------------------------
    // Functions for adding attributes to chains
    NumericMatrix addattribs_patient_chain(arma::mat out);
    List addattribs_pulse_chain(MatrixVector out);
    NumericMatrix addattribs_set_of_pulses(NumericMatrix out);

};


//------------------------------------------------------------
// User-facing functions
//------------------------------------------------------------

// Member Function: Record this iteration
inline void Chains::save_sample(Patient * pat, int iter) {

  int r_iter = iter + 1;

  if ((r_iter > burnin) & ((r_iter % thin) == 0)) {

    // Calculations used repeatedly in this function
    int output_index = ((r_iter - burnin) / thin) - 1;
    double pulsecount = (double)pat->get_pulsecount();

    // Fill patient chain w/ current patient-level estimates
    patient_chain(output_index, 0) = (double)r_iter;
    patient_chain(output_index, 1) = pulsecount;
    patient_chain(output_index, 2) = pat->estimates.baseline_halflife(0);
    patient_chain(output_index, 3) = pat->estimates.mass_mean;
    patient_chain(output_index, 4) = pat->estimates.width_mean;
    patient_chain(output_index, 5) = pat->estimates.baseline_halflife(1);
    patient_chain(output_index, 6) = pat->estimates.errorsq;
    patient_chain(output_index, 7) = pat->estimates.mass_sd;
    patient_chain(output_index, 8) = pat->estimates.width_sd;

    // Create a matrix of current pulse-level estimates and add matrix to the
    //   vector chain

    // Construct matrix of parameters that need constructiong (current r_iter,
    // number of pulses, pulse number id) and add to a matrix
    arma::vec itervec(pulsecount); itervec.fill((double)r_iter);
    arma::vec pcvec(pulsecount); pcvec.fill(pulsecount);
    arma::mat constructed_parms(pulsecount, 3, arma::fill::zeros);
    constructed_parms.col(0) = itervec;
    constructed_parms.col(1) = pcvec;
    constructed_parms.col(2) = arma::linspace<arma::vec>(1, pulsecount, pulsecount);

    // Create matrix of current pulse-level estimates
    arma::mat pulseparms(pulsecount, 5, arma::fill::zeros);
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

// Member Function: Return chains object (R list)
inline List Chains::output(Patient * pat) {

  // Add names to each output chain
  NumericMatrix patient_chain_r = addattribs_patient_chain(patient_chain);
  List pulse_chain_r            = addattribs_pulse_chain(pulse_chains);
  NumericVector fitrange(2); 
  fitrange(0) = pat->data.fitstart;
  fitrange(1) = pat->data.fitend;

  // Create list object combining all chains & other output
  List out = List::create(Named("patient_chain") = patient_chain_r,
                          Named("pulse_chains")  = pulse_chain_r,
                          Named("time_range")     = fitrange);

  return out;

}

// Member Function: Print chain/patient diagnostic info
//   (proposal variance diagnostics handled by MH class)
inline void Chains::print_diagnostic_output(Patient * patient, int iter) {

  if (verbose == 1 && (iter % verbose_iter == 0)) {

    Rcpp::Rcout            << "\n\n\n"                                <<
      "Iteration = "       << iter                                    <<
      " Likelihood = "     << patient->likelihood(false)              << "\n" <<
      "Pulse count = "     << patient->get_pulsecount()               <<
      " Baseline = "       << patient->estimates.baseline_halflife(0) <<
      " Half-life = "      << patient->estimates.baseline_halflife(1) <<
      " Error variance = " << patient->estimates.errorsq              << "\n" <<
      "Mass mean = "       << patient->estimates.mass_mean            <<
      " Mass SD = "        << patient->estimates.mass_sd              <<
      " Width mean = "     << patient->estimates.width_mean           <<
      " Width SD = "       << patient->estimates.width_sd             << "\n" <<

      //" Current pulse-specific parms: " << "\n" << 
      //"Pulse No. Time  Mass  Width\n" << pulse_chains.back() <<
      std::endl;

  }

};


//------------------------------------------------------------
// Supporting/Internal functions
//------------------------------------------------------------

// Member Function: Prep patient_chain for exporting
inline NumericMatrix Chains::addattribs_patient_chain(arma::mat in) {

  // Convert arma obj to Rcpp
  NumericMatrix out = as<NumericMatrix>(wrap(in));

  //Rcout << "Hello world, I'm addattribs_patient_chain()" << std::endl;
  // Add R attributes
  //CharacterVector classes = CharacterVector::create("patient_chain");
  //Rcout << "classes = " << classes << std::endl;
  //out.attr("class") = classes;

  colnames(out) = CharacterVector::create("iteration",
                                          "num_pulses",
                                          "baseline",
                                          "mass_mean",
                                          "width_mean",
                                          "halflife",
                                          "model_error",
                                          "mass_sd",
                                          "width_sd");

  return out;

}

// Member Function: Function for adding attributes to one_set_of_pulses
inline NumericMatrix Chains::addattribs_set_of_pulses(NumericMatrix out) {

  colnames(out) = CharacterVector::create("iteration",
                                          "total_num_pulses",
                                          "pulse_num",
                                          "location",
                                          "mass",
                                          "width",
                                          "eta_mass",
                                          "eta_width"//,
                                          //"lambda"
                                          );
  //out.attr("class") = "one_set_of_pulses";

  return out;

}

//
inline List Chains::addattribs_pulse_chain(std::vector<arma::mat>  in) {

  //MatrixVector::iterator iter = pulse_chains.begin() ;
  //while( iter != iter.end() ){

  List out(in.size());
  //NumericMatrix out;
  int i = 0;
  for (auto &one_iter : in) {  // uses range based loop instead of iterators
    NumericMatrix this_iter = as<NumericMatrix>(wrap(one_iter));
    this_iter = addattribs_set_of_pulses(this_iter);
    out[i] = this_iter;
    //Rcout << "Class of one set of pulses: " << std::endl;
    i++;
  }

  return out;

}

#endif
