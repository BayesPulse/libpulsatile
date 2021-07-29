#ifndef GUARD_bpmod_singlesubject_draw_baselinehalflife_h
#define GUARD_bpmod_singlesubject_draw_baselinehalflife_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>


//
// SS_DrawBaselineHalflife
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the baseline and halflife
//

class SS_DrawBaselineHalflife : 
  public ModifiedMetropolisHastings<Patient, bool, arma::vec, ProposalVariance2p>
{

  public:

    //
    // Constructor
    //   pass the proposal variance parameters to the constructor
    //
    SS_DrawBaselineHalflife(arma::vec in_pv, // double or arma::vec
                            int in_adjust_iter,
                            int in_max_iter,
                            double in_target_ratio,
                            bool verbose,
                            int verbose_iter) :
      ModifiedMetropolisHastings<Patient, bool, arma::vec, ProposalVariance2p>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) { };
    ~SS_DrawBaselineHalflife() {};

  private:

    std::string parameter_name = "Baseline/Half-life";
    std::string get_parameter_name() { return parameter_name; };

    bool parameter_support(arma::vec val, bool *notused) { 
      return ( val(0) > 0.0 && val(1) > 0.0 ); 
    }

    double posterior_function(Patient *patient,
                              arma::vec proposal,
                              bool *notused) {

      double priorb_old       = 0.0;
      double priorh_old       = 0.0;
      double priorb_new       = 0.0;
      double priorh_new       = 0.0;
      double prior_ratio      = 0.0;
      double acceptance_ratio = 0.0;
      double prior_baseline_mean     = patient->priors.baseline_mean;
      double prior_baseline_variance = patient->priors.baseline_variance;
      double prior_halflife_mean     = patient->priors.halflife_mean;
      double prior_halflife_variance = patient->priors.halflife_variance;
      arma::vec current      = patient->estimates.baseline_halflife;
      double curr_likelihood = patient->likelihood(false);

      // Compute ratio of prior densities 
      //  NOTE: this could be matrix algebra -- best if priors were stored
      //  as vectors -- come back to this.
      priorb_old  = current(0) - prior_baseline_mean;
      priorb_old *= 0.5 * priorb_old / prior_baseline_variance;
      priorh_old  = current(1) - prior_halflife_mean;
      priorh_old *= 0.5 * priorh_old / prior_halflife_variance;

      priorb_new  = proposal(0) - prior_baseline_mean;
      priorb_new *= 0.5 * priorb_new / prior_baseline_variance;
      priorh_new  = proposal(1) - prior_halflife_mean;
      priorh_new *= 0.5 * priorh_new / prior_halflife_variance;

      prior_ratio = priorb_old + priorh_old - priorb_new - priorh_new;

      // Update stored estimates to proposal in order to calculate likelihood
      // under proposal
      patient->estimates.baseline_halflife = proposal;

      // Calculate proposed likelihood and then calculate likelihood ratio 
      acceptance_ratio = prior_ratio + (patient->likelihood(false) - curr_likelihood);

      // Reset estimates to current value before exiting (decision made in sample())
      patient->estimates.baseline_halflife = current;

      //Debug output
      //Rcpp::Rcout << "Prior BL mean: " << prior_baseline_mean
      //            << " Prior BL var: " << prior_baseline_variance << "\n"
      //            << "Prior HL mean: " << prior_halflife_mean
      //            << " Prior HL var: " << prior_halflife_variance << "\n"
      //            << "Current: " << current(0) << " " << current(1) << "\n"
      //            << "Acc ratio: " << acceptance_ratio << "\n\n";

      return acceptance_ratio;

    }

};

#endif


