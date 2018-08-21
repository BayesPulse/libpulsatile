#ifndef GUARD_bp_mcmc_metropolishastings_h
#define GUARD_bp_mcmc_metropolishastings_h

#include <Rcpp.h>
#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/proposalvariance.h>
//#include <bp_datastructures/patient.h>
//#include <bp_datastructures/datastructures.h>
#include <bp_mcmc/utils.h>

// metropolishastings.h
//   Abstract template class for defining Modified Metropolis Hastings samplers
//
// Author: Matt Mulvahill
// Created: 10/13/17
// Notes:
//   Outstanding questions:
//    - Where does this class get implemented? -- implementation handled here,
//      instantiation handled in either mcmc function or MMH classes?
//    - draw_proposal() takes SD's - straighten this out.
//    - How do I get priors, parms, data into MMH methods?
//    - Pointer *is --may be-- best.  Is there a meaningful/simple way to take
//      simple args and return a value? Or is it best to update internally
//

// T - sampling unit (e.g.1. pulse iterator)  (e.g.2 patient)
// U - Container of sampling unit  (e.g.1. patient iterator) (e.g.2 population)
// SAMPLETYPE - type of object being sampled (double/int/arma::vec)
// PV - proposal variance type corresponding to SAMPLETYPE (double/int/arma::mat)
template <typename T, typename U, typename SAMPLETYPE, typename PV>
class ModifiedMetropolisHastings
{

  public:
    //T * sampling_unit; // either patient or population class
    //SAMPLETYPE current_val;   // current sample value (double or arma::vec)
    PV pv;                   // needs to be a ProposalVariance object

    // sample from posterior
    //   - runs 1 iteration
    //   - changes internally (pass-by-reference)
    //   - pass Patient as pointer
    //  SAMPLETYPE will be double (or int) or arma::vec depending on single or two
    //    parameter MMH
    //    - pass container (or one level higher in hierarchy) as container
    void sample(T *sampling_unit, SAMPLETYPE *current_val, U *container, int iter) {

      double accept_prob, alpha;

      // Draw new proposal
      SAMPLETYPE proposal = draw_proposal(*current_val, pv.getpsd());
      bool supported      = parameter_support(proposal, container);

      if (!supported) {

        pv.addreject(iter);

      } else {

        accept_prob = posterior_function(sampling_unit, proposal, container);
        alpha = (0 < accept_prob) ? 0 : accept_prob;

        if (log(R::runif(0, 1)) < alpha) {

          pv.addaccept(iter);
          *current_val = proposal;

        } else {

          pv.addreject(iter);

        }
      }
    }

    void sample(T *sampling_unit, SAMPLETYPE *current_val, int iter) {
      U empty_container;
      U * empty = &empty_container;
      sample(sampling_unit, current_val, empty, iter);
    }

  protected:

    //
    // Constructors
    //   - Don't need multiple constructors for now.
    //ModifiedMetropolisHastings() { ProposalVariance pv; }
    //ModifiedMetropolisHastings(PV proposal_variance) {
    //  pv = proposal_variance;
    //}
    ModifiedMetropolisHastings(SAMPLETYPE in_pv, // double or arma::vec
                               int in_adjust_iter,
                               int in_max_iter,
                               double in_target_ratio) :
      pv(in_pv, in_adjust_iter, in_max_iter, in_target_ratio) { }

  private:

    PulseUtils pu;

    double draw_proposal(double current, double proposal_sd) {

      Rcpp::RNGScope rng_scope;
      double new_value = 0;
      new_value = Rf_rnorm(current, proposal_sd);

      return new_value;

    };

    arma::vec draw_proposal(arma::vec current, arma::mat proposal_sd){
      return pu.rmvnorm(current, proposal_sd);
    };

    // Parameter support/Truncation logic
    virtual bool parameter_support(SAMPLETYPE val, U *container) = 0;
    // Posterior function/log(rho) calculation
    virtual double posterior_function(T *sampling_unit, SAMPLETYPE proposal, U *container) = 0;

    // Some estimates can only have one MMH object for multiple objects being
    // sampled -- since pulses are born/die, we can't have an MMH for each pulse
    // -- all pulses need one MMH for each location, mass, and widths.
    // so we need a function for iterating over them, that calls sample() from
    // within it.
    //void sample_pulses(U *container, SAMPLETYPE *current_val); 

};






#endif
