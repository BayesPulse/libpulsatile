#ifndef GUARD_populationchains_h
#define GUARD_populationchains_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/population.h>
#include <bp_datastructures/patient.h>
#include <bp_datastructures/chains.h>

typedef std::vector<arma::mat> MatrixVector;

using namespace Rcpp;

//
// populationchains.h
//   Chains class for population model MCMC output
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// Stores MCMC output for population model with three levels:
//   1. Population-level parameters (μ_α, μ_ω, θ_b, θ_h, υ_α, υ_ω, σ_b, σ_h, σ_α, σ_ω, σ²_e)
//   2. Subject-level parameters (baseline, halflife, mass_mean, width_mean for each subject)
//   3. Pulse-level parameters (location, mass, width, eta_mass, eta_width for each pulse in each subject)
//

class PopulationChains {

    // Validate the MCMC configuration BEFORE any member is initialized. The
    // arma::mat members are sized from num_outputs in the initializer list, so
    // an invalid thin (e.g. 0 -> integer divide-by-zero) or burnin >= iterations
    // (negative size) must be rejected here rather than in the constructor body,
    // which runs only after the members are already (mis-)constructed.
    static int validate_and_count(int iterations, int thin, int burnin) {
      if (thin <= 0) {
        Rcpp::stop("thin must be positive (got %d)", thin);
      }
      if (burnin >= iterations) {
        Rcpp::stop("burnin must be less than iterations (burnin=%d, iterations=%d)",
                   burnin, iterations);
      }
      int n = (iterations - burnin) / thin;
      if (n <= 0) {
        Rcpp::stop("No outputs would be saved with these settings (burnin=%d, iterations=%d, thin=%d)",
                   burnin, iterations, thin);
      }
      return n;
    }

  public:

    // Constructor
    PopulationChains(int in_iterations,
                      int in_thin,
                      int in_burnin,
                      int in_num_subjects,
                      bool in_verbose,
                      int in_verbose_iter)
      : num_outputs(validate_and_count(in_iterations, in_thin, in_burnin))
      , num_subjects(in_num_subjects)
      , population_chain(num_outputs, 12, arma::fill::zeros)  // 11 pop params + iteration
      , subject_chains(in_num_subjects, arma::mat(num_outputs, 6, arma::fill::zeros)) {  // 5 params + iteration per subject

        // MCMC configuration already validated by validate_and_count() above,
        // which runs before the arma::mat members are sized.

        iterations   = in_iterations;
        thin         = in_thin;
        burnin       = in_burnin;
        verbose      = in_verbose;
        verbose_iter = in_verbose_iter;

        // Initialize pulse chains vector for each subject
        for (int s = 0; s < num_subjects; s++) {
          MatrixVector subj_pulse_chains;
          pulse_chains.push_back(subj_pulse_chains);
        }
      }

    //----------------------------------------
    // Member variables
    //----------------------------------------
    int iterations;
    int thin;
    int burnin;
    int num_outputs;
    int num_subjects;
    bool verbose;
    int verbose_iter;

    // Storage for chains
    arma::mat population_chain;                      // Population-level parameters
    std::vector<arma::mat> subject_chains;           // Subject-level parameters (one matrix per subject)
    std::vector<MatrixVector> pulse_chains;          // Pulse-level parameters (vector of matrices per subject)

    //----------------------------------------
    // Primary user-facing functions
    //----------------------------------------
    void save_sample(Population *pop, int iter);
    void print_diagnostic_output(Population *pop, int iter);
    List output(Population *pop);

  private:
    //----------------------------------------
    // Supporting functions
    //----------------------------------------
    NumericMatrix addattribs_population_chain(arma::mat out);
    NumericMatrix addattribs_subject_chain(arma::mat out, int subject_id);
    List addattribs_pulse_chain(MatrixVector out, int subject_id);
    NumericMatrix addattribs_set_of_pulses(NumericMatrix out);
};


//------------------------------------------------------------
// Implementation
//------------------------------------------------------------

// Save one MCMC sample
inline void PopulationChains::save_sample(Population *pop, int iter) {

  int r_iter = iter + 1;

  if ((r_iter > burnin) && ((r_iter % thin) == 0)) {

    int output_index = ((r_iter - burnin) / thin) - 1;

    //
    // 1. Save population-level parameters
    //
    population_chain(output_index, 0)  = (double)r_iter;
    population_chain(output_index, 1)  = pop->estimates.mass_mean;
    population_chain(output_index, 2)  = pop->estimates.width_mean;
    population_chain(output_index, 3)  = pop->estimates.baseline_mean;
    population_chain(output_index, 4)  = pop->estimates.halflife_mean;
    population_chain(output_index, 5)  = pop->estimates.mass_mean_sd;
    population_chain(output_index, 6)  = pop->estimates.width_mean_sd;
    population_chain(output_index, 7)  = pop->estimates.baseline_sd;
    population_chain(output_index, 8)  = pop->estimates.halflife_sd;
    population_chain(output_index, 9)  = pop->estimates.mass_sd;
    population_chain(output_index, 10) = pop->estimates.width_sd;
    population_chain(output_index, 11) = pop->estimates.errorsq;

    //
    // 2. Save subject-level parameters for each subject
    //
    for (int s = 0; s < num_subjects; s++) {
      Patient &subject = pop->subjects[s];

      subject_chains[s](output_index, 0) = (double)r_iter;
      subject_chains[s](output_index, 1) = subject.estimates.baseline_halflife(0);
      subject_chains[s](output_index, 2) = subject.estimates.baseline_halflife(1);
      subject_chains[s](output_index, 3) = subject.estimates.mass_mean;
      subject_chains[s](output_index, 4) = subject.estimates.width_mean;
      subject_chains[s](output_index, 5) = (double)subject.get_pulsecount();

      //
      // 3. Save pulse-level parameters for each pulse in this subject
      //
      int pulsecount = subject.get_pulsecount();

      // Create matrix for this iteration's pulses
      arma::vec itervec(pulsecount); itervec.fill((double)r_iter);
      arma::vec pcvec(pulsecount); pcvec.fill(pulsecount);
      arma::mat constructed_parms(pulsecount, 3, arma::fill::zeros);
      constructed_parms.col(0) = itervec;
      constructed_parms.col(1) = pcvec;
      constructed_parms.col(2) = arma::linspace<arma::vec>(1, pulsecount, pulsecount);

      // Get pulse parameters
      arma::mat pulseparms(pulsecount, 5, arma::fill::zeros);
      int i = 0;
      for (auto &pulse : subject.pulses) {
        pulseparms.row(i) = pulse.get_vector_of_values();
        i++;
      }

      // Combine and save
      arma::mat pulse_mat = arma::join_rows(constructed_parms, pulseparms);
      pulse_chains[s].push_back(pulse_mat);
    }
  }
}


// Print diagnostic output
inline void PopulationChains::print_diagnostic_output(Population *pop, int iter) {

  if (verbose && (iter % verbose_iter == 0) && iter > 0) {

    Rcpp::Rcout << "\n========================================" << std::endl;
    Rcpp::Rcout << "Iteration = " << iter << std::endl;
    Rcpp::Rcout << "========================================" << std::endl;
    Rcpp::Rcout << "Total pulses: " << pop->get_total_pulses() << std::endl;
    Rcpp::Rcout << "\nPopulation parameters:" << std::endl;
    Rcpp::Rcout << "  μ_α (mass mean)      = " << pop->estimates.mass_mean << std::endl;
    Rcpp::Rcout << "  μ_ω (width mean)     = " << pop->estimates.width_mean << std::endl;
    Rcpp::Rcout << "  θ_b (baseline mean)  = " << pop->estimates.baseline_mean << std::endl;
    Rcpp::Rcout << "  θ_h (halflife mean)  = " << pop->estimates.halflife_mean << std::endl;
    Rcpp::Rcout << "  υ_α (mass mean SD)   = " << pop->estimates.mass_mean_sd << std::endl;
    Rcpp::Rcout << "  υ_ω (width mean SD)  = " << pop->estimates.width_mean_sd << std::endl;
    Rcpp::Rcout << "  σ_b (baseline SD)    = " << pop->estimates.baseline_sd << std::endl;
    Rcpp::Rcout << "  σ_h (halflife SD)    = " << pop->estimates.halflife_sd << std::endl;
    Rcpp::Rcout << "  σ_α (pulse mass SD)  = " << pop->estimates.mass_sd << std::endl;
    Rcpp::Rcout << "  σ_ω (pulse width SD) = " << pop->estimates.width_sd << std::endl;
    Rcpp::Rcout << "  σ²_e (error var)     = " << pop->estimates.errorsq << std::endl;
    Rcpp::Rcout << "========================================\n" << std::endl;
  }
}


// Return output list
inline List PopulationChains::output(Population *pop) {

  // Add attributes to population chain
  NumericMatrix pop_chain_r = addattribs_population_chain(population_chain);

  // Add attributes to subject chains
  List subject_chains_r(num_subjects);
  for (int s = 0; s < num_subjects; s++) {
    subject_chains_r[s] = addattribs_subject_chain(subject_chains[s], s + 1);
  }

  // Add attributes to pulse chains
  List pulse_chains_r(num_subjects);
  for (int s = 0; s < num_subjects; s++) {
    pulse_chains_r[s] = addattribs_pulse_chain(pulse_chains[s], s + 1);
  }

  // Create output list
  List out = List::create(
    Named("population_chain") = pop_chain_r,
    Named("subject_chains")   = subject_chains_r,
    Named("pulse_chains")     = pulse_chains_r,
    Named("num_subjects")     = num_subjects
  );

  out.attr("class") = CharacterVector::create("bp_population_fit");

  return out;
}


// Add attributes to population chain
inline NumericMatrix PopulationChains::addattribs_population_chain(arma::mat in) {

  NumericMatrix out = as<NumericMatrix>(wrap(in));

  colnames(out) = CharacterVector::create(
    "iteration",
    "mass_mean",
    "width_mean",
    "baseline_mean",
    "halflife_mean",
    "mass_mean_sd",
    "width_mean_sd",
    "baseline_sd",
    "halflife_sd",
    "mass_sd",
    "width_sd",
    "errorsq"
  );

  return out;
}


// Add attributes to subject chain
inline NumericMatrix PopulationChains::addattribs_subject_chain(arma::mat in, int subject_id) {

  NumericMatrix out = as<NumericMatrix>(wrap(in));

  colnames(out) = CharacterVector::create(
    "iteration",
    "baseline",
    "halflife",
    "mass_mean",
    "width_mean",
    "num_pulses"
  );

  out.attr("subject_id") = subject_id;

  return out;
}


// Add attributes to pulse chain
inline List PopulationChains::addattribs_pulse_chain(MatrixVector in, int subject_id) {

  List out(in.size());

  int i = 0;
  for (auto &one_iter : in) {
    NumericMatrix this_iter = as<NumericMatrix>(wrap(one_iter));
    this_iter = addattribs_set_of_pulses(this_iter);
    out[i] = this_iter;
    i++;
  }

  out.attr("subject_id") = subject_id;

  return out;
}


// Add attributes to one set of pulses
inline NumericMatrix PopulationChains::addattribs_set_of_pulses(NumericMatrix out) {

  colnames(out) = CharacterVector::create(
    "iteration",
    "total_num_pulses",
    "pulse_num",
    "location",
    "mass",
    "width",
    "eta_mass",
    "eta_width"
  );

  return out;
}


#endif
