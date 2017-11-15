#ifndef GUARD_metropolishastings_h
#define GUARD_metropolishastings_h

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "proposalvariance.h"
using namespace Rcpp;

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
//    - Pointer may be best.  Is there a meaningful/simple way to take simple
//      args and return a value? Or is it best to update internally

class ModifiedMetropolisHastings
{

  public:
    // sample from posterior
    //   - runs 1 iteration
    //   - inputs/outputs or changes internally?
    //   - pass Patient as pointer?
    template<typename T> T sample();

  private:
    double draw_proposal(double current, double pv.getpv()) {
      return Rf_rnorm(current, proposal_sd);
    };
    arma::vec draw_proposal(arma::vec current, arma::mat pv.getpv()){
      return Rf_rmultinorm(current, proposal_sd);
    };
    virtual parameter_support();  // i.e. truncation logic
    virtual posterior_function(); // logrho_calculation
    virtual pv;                   // needs to be a ProposalVariance object

};


// MCMC sampling function
//  T will be double or arma::vec depending on single or two parameter MMH
template <typename T>
T ModifiedMetropolisHastings::sample()
{

  double accept_prob, alpha;

  // Draw new proposal
  T proposal     = draw_proposal(current, pv.getpv());
  bool supported = parameter_support(proposal);

  if (!supported) {

    pv.addreject();
    return current;

  } else {

    accept_prob = posterior_function();
    alpha = (0 < accept_prob) ? 0 : accept_prob;

    if (log(Rf_runif()) < alpha) {

      pv.addaccept();
      return proposal;

    } else {

      pv.addreject();
      return current;

    }
  }
}




////////////////////////////////////////////////////////////
// Change these over to classes inheriting abstract class
////////////////////////////////////////////////////////////


class ModifiedMetropolisHastings : public sd_widths {

};


class ModifiedMetropolisHastings : public sd_mass {

};

class ModifiedMetropolisHastings : public baseline {

};

class ModifiedMetropolisHastings : public halflife {

};


#endif
