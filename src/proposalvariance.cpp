// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "proposalvariance.h"

// -------------------------------------
// initialize_proposals() for each type
// -------------------------------------
// One for one-variable and one for two-variable MMH/pv.
// create/set pv and psd objects (proposal var-covar and correlation matrices)
void ProposalVariance::initialize_proposals(double this_pv)
{

  set_proposals(this_pv);

}

// -------------------------------------
// set_proposals() for each type
// -------------------------------------
void ProposalVariance::set_proposals(double this_pv)
{

  pv = this_pv;
  psd = sqrt(pv);

}


// -------------------------------------
// adjustpv() for each type
// -------------------------------------
// Acceptance adjustment routine for one-parameter modified MH
void ProposalVariance::adjustpv()
{

  double y = 1.0 + 1000.0 * pow(getratio() - target_ratio, 3);
  if (y < 0.9)      set_proposals(pv * 0.9);
  else if (y > 1.1) set_proposals(pv * 1.1);

}
