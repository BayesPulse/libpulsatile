void draw_re_sd(Node_type *list,
                Priors *priors, 
                Common_parms *parms, 
                double v1, 
                double v2,
                long *arevm, long *nrevm, long *arevw, long *nrevw) {
  int j, 
      num_pulses;
  long *accept_counter;   
  double *new_sd, 
         stdx_old, 
         stdx_new, 
         new_int, 
         old_int, 
         first_part,
         second_part, 
         third_part, 
         log_rho;
  Node_type *node;

  // Allocate Memory
  new_sd         = (double *)calloc(2, sizeof(double));
  accept_counter = (long *)calloc(2, sizeof(long));

  // Add 1 to the counters for acceptance rates of sigma_a and sigma_w
  (*nrevm)++;
  (*nrevw)++;

  // Assign current acceptance counts to temporary vector
  accept_counter[0] = *arevm;
  accept_counter[1] = *arevw;

  // Draw proposed values for sigma_a and sigma_w
  new_sd[0] = Rf_rnorm(parms->re_sd[0], v1);
  new_sd[1] = Rf_rnorm(parms->re_sd[1], v2);

  // Accept or Reject sigma_a, then accept or reject for sigma_w
  for (j = 0; j < 2; j++) {

    // We only can accept the proposed value if it is positive
    if (new_sd[j] > 0 && new_sd[j] < priors->re_sdmax[j]) {

      node = list->succ;
      num_pulses = 0;
      third_part = 0;
      old_int = new_int = 0;
      while (node != NULL) {

        // Normalizing constants
        stdx_old   = parms->theta[j] / ( parms->re_sd[j] / sqrt(node->eta[j]) );
        stdx_new   = parms->theta[j] / ( new_sd[j]       / sqrt(node->eta[j]) );
        new_int   += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
        old_int   += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);

        // Count pulses
        num_pulses++;
        // 3rd 'part' of acceptance ratio
        third_part += node->eta[j] * 
                      (node->theta[j] - parms->theta[j]) *
                      (node->theta[j] - parms->theta[j]);

        // Next pulse
        node = node->succ;
      }

      // 1st and 2nd 'parts' of acceptance ratio
      first_part  = (num_pulses) * (log(parms->re_sd[j]) - log(new_sd[j]));
      second_part = 0.5 * 
        ((1 / (parms->re_sd[j] * parms->re_sd[j])) - (1 / (new_sd[j] * new_sd[j])));

      // Compute log rho, and set alpha equal to min(log rho,0)
      log_rho = old_int - new_int + first_part + second_part * third_part;
      log_rho = fmin(0, log_rho);

      // If log(U) < log rho, accept the proposed value
      if (log(Rf_runif(0, 1)) < log_rho) {
        accept_counter[j]++;
        parms->re_sd[j] = new_sd[j];
      }

    }

  } 

  // Set acceptance count equal to temp vector components
  *arevm = accept_counter[0];
  *arevw = accept_counter[1];

  free(new_sd);
  free(accept_counter);

}
