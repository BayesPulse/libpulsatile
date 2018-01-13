//#ifndef GUARD_metropolishastings_h
//#define GUARD_metropolishastings_h
//
//#include "proposalvariance.h"
//#include "utils.h"
//#include <RcppArmadillo.h>
//#include <RInside.h>
//
//// metropolishastings.h
////   Abstract class for defining Modified Metropolis Hastings samplers
////
//// Author: Matt Mulvahill
//// Created: 10/13/17
//// Notes:
////   Outstanding questions:
////    - Where does this class get implemented? -- implementation handled here,
////      instantiation handled in either mcmc function or MMH classes?
////    - draw_proposal() takes SD's - straighten this out.
////    - How do I get priors, parms, data into MMH methods?
////    - Pointer *is --may be-- best.  Is there a meaningful/simple way to take
////      simple args and return a value? Or is it best to update internally
////
//
//template <typename T>
//class ModifiedMetropolisHastings
//{
//
//  public:
//    // sample from posterior
//    //   - runs 1 iteration
//    //   - inputs/outputs or changes internally?
//    //   - pass Patient as pointer?
//    //  S will be double or arma::vec depending on single or two parameter MMH
//    template <typename S> S sample()
//    {
//
//      double accept_prob, alpha;
//
//      // Draw new proposal
//      S proposal     = draw_proposal(current, pv.getpv());
//      bool supported = parameter_support(proposal);
//
//      if (!supported) {
//
//        pv.addreject();
//        return current;
//
//      } else {
//
//        accept_prob = posterior_function();
//        alpha = (0 < accept_prob) ? 0 : accept_prob;
//
//        if (log(Rf_runif()) < alpha) {
//
//          pv.addaccept();
//          return proposal;
//
//        } else {
//
//          pv.addreject();
//          return current;
//
//        }
//      }
//    }
//
//
//  protected:
//    ModifiedMetropolisHastings(T proposal_variance);
//    ModifiedMetropolisHastings(T proposal_variance,
//                               int adjust_at_iter,
//                               int max_iters);
//
//  private:
//    double draw_proposal(double current, double proposal_sd) {
//      return Rf_rnorm(current, proposal_sd);
//    };
//    arma::vec draw_proposal(arma::vec current, arma::mat proposal_sd){
//      return pulseutils::rmvnorm(current, proposal_sd);
//    };
//    virtual bool parameter_support();  // i.e. truncation logic
//    virtual double posterior_function(); // logrho_calculation
//    ProposalVariance pv;                   // needs to be a ProposalVariance object
//
//};
//
//
//
//
//
//
//#endif
