#include <cmath>

// Constructor
void ProposalVariance::ProposalVariance(double initial_pv,
                                        int adjust_at_iter,
                                        int max_iters = 25000); {

  pv             = initial_pv;
  adjust_at_iter = adjust_at_iter;
  max_iters      = max_iters;
  accept_ct      = 0;
  iter_ct        = 0;
  ratio          = 0;

}

// destructor -- need to define here?
ProposalVariance::~ProposalVariance() { }

void ProposalVariance::adjustpv()
{

}




// Acceptance adjustment routine for one-parameter modified MH
// args: x   = current acceptance rate
//       *X  = current proposal variance
// returns: none -- updates internally
void ProposalVariance::adjustpv(double x, double *X, double target_pv = 0.35)
{

  double y = 1.0 + 1000.0 * pow(target_pv, 3);
  if (y < 0.9) {
    y = 0.9;
    *X *= y;
  } else if (y > 1.1) {
    y = 1.1;
    *X *= y;
  }

}


// Acceptance adjustment routine for two-parameter modified MH
//   used previously for baseline and halflife
// args: x    = current acceptance rate
//       **X  = current proposal var-covar matrix
//       corr = correlation between the two proposal variances
// returns: none -- updates internally
void ProposalVariance::adjustpv(double x, double **X, double corr, double target_pv = 0.25)
{

  // y - new diagonal elements of proposal variance-covariance matrix based on
  // inputs
  double y = 1.0 + 1000.0 * pow(target_pv, 3);
  if (y < .90) {
    y = .90;
    X[0][0] *= y;
    X[1][1] *= y;
    X[0][1]  = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
  } else if (y > 1.1) {
    y = 1.1;
    X[0][0] *= y;
    X[1][1] *= y;
    X[0][1]  = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
  }

}

