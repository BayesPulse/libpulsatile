#ifndef GUARD_chains_h
#define GUARD_chains_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include "bp_datastructures/patient.h"
#include "bp_datastructures/population.h"
#include "bp_datastructures/patientdata.h"
#include "bp_datastructures/patientestimates.h"
#include "bp_datastructures/patientpriors.h"
#include "bp_datastructures/pulseestimates.h"

typedef std::vector<std::vector<arma::mat>> MatrixVectorVector;
typedef std::vector<arma::mat> MatrixVector;

using namespace Rcpp;

class PopChains {

  public:
    // Base constructor (population, single hormone)
    PopChains(int in_iterations, int in_thin, int in_burnin,
           bool in_response_hormone, bool in_verbose, 
           int in_verbose_iter, int in_num_patients)
      : num_outputs((in_iterations - in_burnin) / in_thin)
      , pop_chain(num_outputs, 11, arma::fill::zeros) {

      iterations = in_iterations;
      thin = in_thin;
      burnin = in_burnin;
      response_hormone = in_response_hormone;
      verbose = in_verbose;
      verbose_iter = in_verbose_iter;
      num_patients = in_num_patients;

      for (int i = 0; i < num_patients; i++) {
        patient_chains.push_back(arma::mat(num_outputs, 6,
                                           arma::fill::zeros));
      }

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
    int num_patients;
    MatrixVectorVector pulse_chains;
    MatrixVector patient_chains;
    arma::mat pop_chain;

    //----------------------------------------
    // Primary user-facing functions
    //----------------------------------------
    // Record this iteration -- adds 1 iteration to Chains object
    void save_sample(Population * pop, int iter);
    // Print pulse/patient diagnostic info
    void print_diagnostic_output(Population * pop, int iter);
    // Return chains function
    List output(Population * pop);

  private: 
    //----------------------------------------
    // Supporting/Internal functions
    //----------------------------------------
    // Functions for adding attributes to chains
    NumericMatrix addattribs_patient_chain(arma::mat out);
    NumericMatrix addattribs_pop_chain(arma::mat out);
    List addattribs_pulse_chain(MatrixVector out);
    NumericMatrix addattribs_set_of_pulses(NumericMatrix out);

};

//------------------------------------------------------------
// User-facing functions
//------------------------------------------------------------

// Member Function: Record this iteration
void PopChains::save_sample(Population * pop, int iter) {

  int r_iter = iter + 1;


  if ((r_iter > burnin) & ((r_iter % thin) == 0)) {

    // Calculations used repeatedly in this function
    int output_index = ((r_iter - burnin) / thin) - 1;
    int pat_iter = 0;

    // Fill population chain w/ current population-level estimates
    pop_chain(output_index, 0) = (double)r_iter;
    pop_chain(output_index, 1) = pop->patPriors.mass_p2p_sd;
    pop_chain(output_index, 2) = pop->patPriors.width_p2p_sd;
    pop_chain(output_index, 3) = pop->patPriors.width_mean;
    pop_chain(output_index, 4) = pop->patPriors.mass_mean;
    pop_chain(output_index, 5) = pop->patPriors.baseline_mean;
    pop_chain(output_index, 6) = pop->patPriors.halflife_mean;
    pop_chain(output_index, 7) = pop->patPriors.width_s2s_sd;
    pop_chain(output_index, 8) = pop->patPriors.mass_s2s_sd;
    pop_chain(output_index, 9) = pop->patPriors.baseline_sd;
    pop_chain(output_index, 10) = pop->patPriors.halflife_sd;

    for(auto pat : pop->patients) {

      double pulsecount = (double)pat.get_pulsecount();

      // Fill patient chain w/ current patient-level estimates
      patient_chains[pat_iter](output_index, 0) = (double)r_iter;
      patient_chains[pat_iter](output_index, 1) = pulsecount;
      patient_chains[pat_iter](output_index, 2) = pat.estimates.mass_mean;
      patient_chains[pat_iter](output_index, 3) = pat.estimates.width_mean;
      patient_chains[pat_iter](output_index, 4) = pat.estimates.errorsq;
      patient_chains[pat_iter](output_index, 5) = pat.likelihood(false);

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
      for (auto pulse : pat.pulses) {
        pulseparms.row(i) = pulse.get_vector_of_values();
        i++;
      }

      // Combine constructed (ids) and pulse-level estimates and add to vector of
      // matrices
      arma::mat new_mat = arma::join_rows(constructed_parms, pulseparms);
      std::vector<arma::mat> temp;
      pulse_chains.push_back(temp);
      pulse_chains[pat_iter].push_back(new_mat);

      pat_iter++;

    }

  }

};

// Member Function: Return chains object (R list)
List PopChains::output(Population * pop) {

  int patientCount = pop->get_patientcount();

  // Add names to each output chain
  NumericMatrix pop_chain_r     = addattribs_pop_chain(pop_chain);
  List pulse_chains_r           = List(patientCount);
  List patient_chains_r         = List(patientCount);

  for(int i = 0; i < patientCount; i++) {

    NumericMatrix patient_chain_r = addattribs_patient_chain(patient_chains[i]);
    List pulse_chain_r            = addattribs_pulse_chain(pulse_chains[i]);

    String patientString("Patient ");
    patientString += std::to_string(i+1);    

    patient_chains_r[i] = patient_chain_r;
    pulse_chains_r[i]   = pulse_chain_r;

  }

  //NumericMatrix patient_chain_r = addattribs_patient_chain(patient_chains[2]);
  //List pulse_chain_r            = addattribs_pulse_chain(pulse_chains[1]);

  NumericVector fitrange(2); 
  fitrange(0) = pop->patients[1].data.fitstart;
  fitrange(1) = pop->patients[1].data.fitend;

  // Create list object combining all chains & other output
  List out = List::create(Named("pop_chain") = pop_chain_r,
                          Named("patient_chains") = patient_chains_r,
                          Named("pulse_chains")  = pulse_chains_r,
                          Named("time_range")     = fitrange);

  return out;

}

// Member Function: Print chain/patient/population diagnostic info 
//   (proposal variance diagnostics handled by MH class)
void PopChains::print_diagnostic_output(Population * pop, int iter) {

  if (verbose == 1 && (iter % verbose_iter == 0)) {

  int pat_iter = 0;

  Rcpp::Rcout              << "\n\n\n"                         <<
    "Population Level------------------------------------"     << "\n" <<
    "Iteration = "         << iter                             <<
    " Mass P2P sd = "      << pop->patPriors.mass_p2p_sd       <<
    " Width P2P sd = "     << pop->patPriors.width_p2p_sd      <<
    " Width mean = "       << pop->patPriors.width_mean        << "\n" <<
    "Mass mean = "         << pop->patPriors.mass_mean         <<
    " Baseline mean = "    << pop->patPriors.baseline_mean     <<
    " Halflife mean = "    << pop->patPriors.halflife_mean     <<
    " Width S2S sd = "     << pop->patPriors.width_s2s_sd      << "\n" <<
    "Mass S2S sd = "       << pop->patPriors.mass_s2s_sd       <<
    " Baseline S2S sd = "  << pop->patPriors.baseline_sd       <<
    " Halflife S2S sd = "  << pop->patPriors.halflife_sd       << "\n\n";

  for(auto pat : pop->patients) {

      Rcpp::Rcout            <<
        "Patient = "         << pat_iter                       <<
        " Likelihood = "     << pat.likelihood(false)          << "\n" <<
        "Pulse count = "     << pat.get_pulsecount()           <<
        " Error variance = " << pat.estimates.errorsq          << "\n" <<
        "Mass mean = "       << pat.estimates.mass_mean        <<
        " Width mean = "     << pat.estimates.width_mean       << "\n\n";

        //" Current pulse-specific parms: " << "\n" << 
        //"Pulse No. Time  Mass  Width\n" << pulse_chains.back() <<

        pat_iter++;
    }

  }

};


//------------------------------------------------------------
// Supporting/Internal functions
//------------------------------------------------------------

// Member Function: Prep population chain for exporting
NumericMatrix PopChains::addattribs_pop_chain(arma::mat in) {

  // Convert arma obj to Rcpp
  NumericMatrix out = as<NumericMatrix>(wrap(in));

  //Rcout << "Hello world, I'm addattribs_pop_chain()" << std::endl;

  colnames(out) = CharacterVector::create("iteration",
                                          "mass_p2p_sd",
                                          "width_p2p_sd",
                                          "width_mean",
                                          "mass_mean",
                                          "baseline_mean",
                                          "halflife_mean",
                                          "width_s2s_sd",
                                          "mass_s2s_sd",
                                          "baseline_sd",
                                          "halflife_sd"

  );

  return out;

}

// Member Function: Prep patient_chain for exporting
NumericMatrix PopChains::addattribs_patient_chain(arma::mat in) {

  // Convert arma obj to Rcpp
  NumericMatrix out = as<NumericMatrix>(wrap(in));

  //Rcout << "Hello world, I'm addattribs_patient_chain()" << std::endl;
  // Add R attributes
  //CharacterVector classes = CharacterVector::create("patient_chain");
  //Rcout << "classes = " << classes << std::endl;
  //out.attr("class") = classes;

  colnames(out) = CharacterVector::create("iteration",
                                          "num_pulses",
                                          "mass_mean",
                                          "width_mean",
                                          "model_error",
                                          "likelihood");

  return out;

}

// Member Function: Function for adding attributes to one_set_of_pulses
NumericMatrix PopChains::addattribs_set_of_pulses(NumericMatrix out) {

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
List PopChains::addattribs_pulse_chain(MatrixVector in) {

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

//struct PopulationChain { NumericMatrix population; };


  //void save_sample(Population * pop);
    //arma::mat population_chain;



    // Population constructor
    // Chains(int iterations, int thin, int burnin, bool response_hormone, int num_patients) :
    //   Chains(iterations, thin, burnin, response_hormone) {
    //   association_chain(num_outputs, 4);
    //  model_type = "population";
    //   }



  //// Function for adding attributes to population_chain
  //// TODO: Straighten out variance vs SD terms and why is there a variance AND
  //// SD term for mass/width
  //NumericMatrix addattribs_population_chain(NumericMatrix out);




//void save_sample(Population * pop) {

//};

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
