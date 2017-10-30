
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





/*----------------------------------------------------------------------------**
/  Acceptance adjustment routines from poppulsepaper/Pop-UnifLogNormPrior
/-----------------------------------------------------------------------------*/

/*********************************************************************/
              /*START OF adjust_acceptance SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*adjust_acceptance: this adjusts the proposal variances based on the inputted
                    acceptance rate and proposal variance. If the acceptance rate is
                    too high or too low, proposal variance will be adjusted.
    ARGUMENTS: double x; the inputted acceptance rate; usually inputted as the
                  acceptance counter divided by the attempt counter
               double *X; the current proposal variance
    RETURNS: None; update to the proposal variance is made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 y: new proposal variance based on inputs

SUBROUTINES USED
 None  */
 /**************************************************************************/

void adjust_acceptance(double x, double *X)
{

  double y = 1.0 + 1000.0 * (x - 0.35) * (x - 0.35) * (x - 0.35);
  if (y < 0.9) y = 0.9;
  if (y > 1.1) y = 1.1;

  *X *= y;

}

/*********************************************************************/
              /*START OF adjust2_acceptance SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*adjust2_acceptance: this adjusts the proposal variance-covriance matrix for
                     baseline and halflife based on the inputted acceptance rate
                     and proposal matrix. If the acceptance rate is too high or
                     too low, proposal variance will be adjusted.
    ARGUMENTS: double x; the inputted acceptance rate; usually inputted as the
                  acceptance counter divided by the attempt counter
               double **X; the current proposal variance-covariance matrix
               double corr; the correlation between the two proposal variances
    RETURNS: None; update to the proposal variance is made internally*/
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 y: new diagonal elements of proposal variance-covariance matrix based on inputs

SUBROUTINES USED
 None  */
 /**************************************************************************/

void adjust_acceptance(double x, double **X, double corr)
{

  double y = 1. + 1000.*(x-.25)*(x-.25)*(x-.25);

  if (y < .90) {
    y = .90;
    X[0][0] *= y;
    X[1][1] *= y;
    X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
  }

  if (y > 1.1) {
    y = 1.1;
    X[0][0] *= y;
    X[1][1] *= y;
    X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
  }

}

