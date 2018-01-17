//-----------------------------------------------------------------------------
// START OF draw_fixed_effects SUBROUTINE*/
//
//
// draw_fixed_effects: this runs the Gibbs sampler draw for subject-specific
// mean masses and widths;
// ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
// Priors *priors; the current values of the prior parameters;
// Common_parms *parms; the current values of the common parameters;
// unsigned long *seed; seed values needed for the randon number generator;
// RETURNS: None; all updates are made internally
//
//********************************************************************/
//
// VARIABLE DEFINITIONS
// j: generic counter
// numnode: counter of pulses
//
// gvar: variance of distribution used in Gibbs sampler
// *subject: current list of subjects and their characteristics
// *subnode: for a given subject, list of pulses and their characteristics
//
// SUBROUTINES USED
// kiss: found in randgen.c; draws from U(0,1) distribution
// rnorm: found in randgen.c; draws from the normal distribution
//------------------------------------------------------------------------------
void draw_fixed_effects(Node_type *list, 
                        Priors *priors, 
                        Common_parms *parms, 
                        double sdfem, 
                        double sdfew,
                        long *afem, long *nfem, long *afew, long *nfew ) {  
  int j, newint, oldint, numnode;
  long *accept_counter;
  double normalizing_ratio, acceptance_ratio, theta[2], old_prior, new_prior,
         prior_ratio, alpha, psum_old,
         psum_new, prop_ratio, stdxold, stdxnew;
  Node_type *node;

  accept_counter = (long *)calloc(2, sizeof(long));

  // Acceptance counters
  accept_counter[0]  = *afem;
  accept_counter[1]  = *afew;
  // Iteration counters
  (*nfem)++;
  (*nfew)++;

  // Proposed values
  theta[0] = Rf_rnorm(parms->theta[0], sdfem);
  theta[1] = Rf_rnorm(parms->theta[1], sdfew);

  // Draw a new pair of pulse mass and width on the subject level
  for (j = 0; j < 2; j++) {

    if (theta[j]>0) {

      // Prior Ratio  -- DONE
      old_prior    = (parms->theta[j] - priors->fe_mean[j]) * 
                     (parms->theta[j] - priors->fe_mean[j]); 
      new_prior    = (theta[j] - priors->fe_mean[j]) * 
                     (theta[j] - priors->fe_mean[j]); 
      prior_ratio  = (old_prior - new_prior) / (2 * priors->fe_variance[j]);

      // likelihood ratio (Proposal Ratio)
      psum_old = 0;
      psum_new = 0;
      newint   = 0;
      oldint   = 0;

      numnode = 0;
      node = list->succ;

      // 'likelihood' ratio -- Ratio of p(alpha|mu, nu, kappa)
      while (node != NULL) {
        psum_old += (node->theta[j] - parms->theta[j]) * 
                    (node->theta[j] - parms->theta[j]) * node->eta[j];
        psum_new += (node->theta[j] - theta[j]) * 
                    (node->theta[j] - theta[j]) * node->eta[j];

        // Normalizing constants
        stdxnew   = theta[j]        * sqrt(node->eta[j]) / parms->re_sd[j];
        stdxold   = parms->theta[j] * sqrt(node->eta[j]) / parms->re_sd[j];
        oldint += Rf_pnorm5(stdxold, 0, 1, 1.0, 1.0); // second 1.0 does the log xform for us 
        newint += Rf_pnorm5(stdxnew, 0, 1, 1.0, 1.0); // first 1.0 says to use lower tail

        node = node->succ;
        numnode++;
      }

      prop_ratio = 0.5 / (parms->re_sd[j] * parms->re_sd[j]) * 
                   (psum_old - psum_new);
      normalizing_ratio = oldint - newint; 

      acceptance_ratio = prior_ratio + prop_ratio + normalizing_ratio;
      alpha = (0 < acceptance_ratio) ? 0 : acceptance_ratio;

      // If log(U) < log rho, accept the proposed value
      // Increase acceptance count by 1
      if (log(Rf_runif(0, 1)) < alpha) {
        accept_counter[j]++;
        parms->theta[j] = theta[j];
      } 

    } 

  } // end of loop through mass & width 

  *afem = accept_counter[0];
  *afew = accept_counter[1];
  free(accept_counter);

}
