
//-----------------------------------------------------------------------------
//
// draw_random_effects:
//   this runs the M-H draw for individual pulse masses and widths
//
//   NOTE: now uses a truncated t prior on alpha_i and omega_i
//
//   ARGUMENTS:
//     double **ts         - this is the matrix of observed data (a
//                           column of times and a column of log(concentration)
//     Node_type *list     - this is the current list of pulses that exist
//     Common_parms *parms - the current values of the common parameters
//     int N               - the number of observations in **ts
//     double *like        - the current value of the likelihood
//     double v1           - the proposal variance for individual pulse masses
//     double v2           - the proposal variance for individual pulse widths
//-----------------------------------------------------------------------------
void draw_random_effects(double **ts, 
                         Node_type *list, 
                         Common_parms *parms, 
                         int N, 
                         double *like,
                         double v1, 
                         double v2, 
                         long *arem, long *nrem, 
                         long *arew, long *nrew) {
  //double current_like;
  int i;                  // Generic counters
  int j;                  // Generic counters
  double logrho;          // Acceptance ratio prior to min(0,logrho)
  double prior_old;       // Prior portion of accept. ratio due to current value
  double prior_new;       // Prior portion of accept. ratio due to proposed value
  double old_val;         // For holding old value of mass/width
  long *accept_counter; // Internal array for counting acceptances
  double *pRE;            // Array of proposed RE mass/width value
  double prior_ratio;     // Portion of accept. ratio due to priors
  double like_ratio;      // Portion of accept. ratio due to likelihoods
  double alpha;           // Acceptance ratio after min(0,logrho)
  double plikelihood;     // Likelihood with proposed value
  double *old_contrib;    // For holding mean contrib using current value
  Node_type *node;        // Pointer to linklist of pulses


  // Go to start of pulses 
  node = list->succ;

  // Set acceptance counts equal to temporary vector
  accept_counter    = (long *)calloc(2, sizeof(long));
  accept_counter[0] = *arem;
  accept_counter[1] = *arew;

  // Go through each existing pulse 
  while (node != NULL) {

    // Allocate Memory
    old_contrib = (double *)calloc(N, sizeof(double));
    pRE         = (double *)calloc(2, sizeof(double));

    // Increase the denominators of the acceptance rates
    (*nrem)++;
    (*nrew)++;

    // Draw proposed values of current pulse's mass and width
    pRE[0] = Rf_rnorm(node->theta[0], v1);
    pRE[1] = Rf_rnorm(node->theta[1], v2);

    if (pRE[0] > 0.0 && pRE[1] > 0.01 && pRE[1] < 100){
      // Determine if we accept or reject proposed pulse mass then determine
      // if we accept or reject proposed pulse width
      for (j = 0; j < 2; j++){

        // Compute the log of the ratio of the priors
        prior_old = node->theta[j] - parms->theta[j];
        prior_old *= 0.5*prior_old;
        prior_new = pRE[j] - parms->theta[j];
        prior_new *= 0.5*prior_new;
        prior_ratio = prior_old - prior_new;
        prior_ratio /= parms->re_sd[j];
        prior_ratio /= parms->re_sd[j];

        // Save the current value of mass/width
        old_val = node->theta[j];

        // Set the pulse's mass/width equal to the proposed value
        node->theta[j] = pRE[j];

        // Save the mean_contrib for that pulse
        for (i = 0; i < N; i++) {
          old_contrib[i] = node->mean_contrib[i];
        }

        // Recalculate that pulse's mean_contrib assuming proposed mass/width 
        mean_contribution(node, ts, parms, N);

        // Calculate likelihood assuming proposed mass/width 
        plikelihood = likelihood(list, ts, parms, N, list);

        like_ratio = plikelihood - *like;

        // Compute the log of the ratio between the two likelihoods

        // Calculate log rho; set alpha equal to min(0,log rho) 
        alpha = (0 < (logrho = (prior_ratio + like_ratio))) ? 0:logrho;

        // If log U < log rho, accept the proposed value, increase acceptance
        // counter 
        if (log(Rf_runif(0, 1)) < alpha) {
          accept_counter[j]++;
        } else {
        // Otherwise, reject the proposed value, set pulse's mass/width back to 
        // saved current value, and set pulse's mean_contrib back to saved value

          node->theta[j] = old_val;
          for (i = 0; i < N; i++) {
            node->mean_contrib[i] = old_contrib[i];
          }

        }
      } // end of loop through mass & width 
    }

    // Advance to next pulse
    node = node->succ;
    free(pRE);
    free(old_contrib);

  } // end of loop through pulses

  // Set counter equal to temporary vector values
  *arem = accept_counter[0];
  *arew = accept_counter[1];
  // free memory
  free(accept_counter);

}
