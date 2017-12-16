// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "counter.h"

#ifndef GUARD_proposalvariance_h
#define GUARD_proposalvariance_h

//
// proposalvariance.h
//   Definitions for the ProposalVariance class
//
// Author: Matt Mulvahill
// Created: 10/13/17
//
// NOTE: Sec9.3.1 in Acc C++ define here to tell compiler to avoid function-call
// overhead, so all simple fn's are defined here -- not sure if this works for
// non-const functions or virtual functions.
//
// NOTE: ProposalVariance has a constructor only to document the requirements
// for the constructor.  It can't actually be used since its an abstract class.
//

class ProposalVariance {

  public:
    ProposalVariance(double inpv);
    ~ProposalVariance();
    double getpv() const { return pv; };

  private:
    double pv;
    Counter counter;
};

ProposalVariance::ProposalVariance(double inpv)
  : pv(inpv)
{
  counter();
}

#endif
