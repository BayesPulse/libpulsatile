#ifndef GUARD_ss_draw_baselinehalflife_h
#define GUARD_ss_draw_baselinehalflife_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"


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
                            double in_target_ratio) :
      ModifiedMetropolisHastings<Patient, bool, arma::vec, ProposalVariance2p>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio) { };
    ~SS_DrawBaselineHalflife() {};

//  private:
//
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
      double prior_baseline_mean     = patient->priors->baseline_mean;
      double prior_baseline_variance = patient->priors->baseline_variance;
      double prior_halflife_mean     = patient->priors->halflife_mean;
      double prior_halflife_variance = patient->priors->halflife_variance;
      arma::vec current      = patient->estimates->baseline_halflife;
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
      patient->estimates->baseline_halflife = proposal;

      // Calculate proposed likelihood and then calculate likelihood ratio 
      acceptance_ratio = prior_ratio + (patient->likelihood(false) - curr_likelihood);

      // Reset estimates to current value before exiting (decision made in sample())
      patient->estimates->baseline_halflife = current;

      return acceptance_ratio;

    }

};

#endif


