#ifndef GUARD_metropolishastings_h
#define GUARD_metropolishastings_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include "proposalvariance.h"
#include "utils.h"

// metropolishastings.h
//   Abstract class for defining Modified Metropolis Hastings samplers
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

// T = sampling unit = Patient or Population (or pulse?) 
// S = current val type
// PV = proposal variance object type
template <typename T, typename S, typename PV>
class ModifiedMetropolisHastings
{

  public:
    //T * sampling_unit; // either patient or population class
    //S current_val;   // current sample value (double or arma::vec)
    PV pv;                   // needs to be a ProposalVariance object

    // sample from posterior
    //   - runs 1 iteration
    //   - inputs/outputs or changes internally?
    //   - pass Patient as pointer?
    //  S will be double (or int) or arma::vec depending on single or two
    //    parameter MMH
    void sample(T *sampling_unit, S *current_val) {

      double accept_prob, alpha;
      std::cout << "current_val is = " << *current_val << std::endl;

      // Draw new proposal
      std::cout << "Proposal variance and SD is = " << pv.getpv() << " & " << pv.getpsd() << 
        std::endl;
      S proposal     = draw_proposal(current_val, pv.getpsd());
      bool supported = parameter_support(proposal);
      std::cout << "Proposal is = " << proposal << std::endl;

      if (!supported) {

        pv.addreject();
        std::cout << "reject (not supported) = " << pv.getratio() << std::endl;
        //return current_val;

      } else {

        accept_prob = posterior_function(sampling_unit, proposal);
        alpha = (0 < accept_prob) ? 0 : accept_prob;

        if (log(R::runif(0, 1)) < alpha) {

          pv.addaccept();
          *current_val = proposal;
          std::cout << "accepted! = " << pv.getratio() << std::endl;
          //return proposal;

        } else {

          pv.addreject();
          std::cout << "reject = " << pv.getratio() << std::endl;
          //return current_val;

        }
      }
    }


  protected:
    // Constructors
    ModifiedMetropolisHastings() {}
    ModifiedMetropolisHastings(PV proposal_variance);
    ModifiedMetropolisHastings(double in_pv,
                               int in_adjust_iter,
                               int in_max_iters,
                               double in_target_ratio);

  private:
    PulseUtils pu;
    double draw_proposal(const double * current, double proposal_sd) {
      double p = Rf_rnorm(*current, proposal_sd);
      std::cout << "(inside draw_proposal)  proposal_sd is = " << proposal_sd << std::endl;
      std::cout << "(inside draw_proposal)  current is = " << *current << std::endl;
      std::cout << "(inside draw_proposal)  proposal is = " << p << std::endl;
      return p; 
    };
    arma::vec draw_proposal(const arma::vec * current, arma::mat proposal_sd){
      return pu.rmvnorm(current, proposal_sd);
    };
    virtual bool const parameter_support(S val); // i.e. truncation logic
    virtual double const posterior_function(T * sampling_unit, S proposal);   // logrho_calculation

};






#endif
