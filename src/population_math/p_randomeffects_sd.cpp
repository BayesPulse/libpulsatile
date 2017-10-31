
//void draw_re_sd(Node_type *list,
//                Priors *priors, 
//                Common_parms *parms)


#include "mh.h"

class mmh_sd_widths : public ModifiedMetropolisHastings<mmh_sd_widths> {
  double draw_proposal();
}






// Draw proposed values for sigma_a and sigma_w
double mmh_sd_widths::draw_proposal(double current_sd, double proposal_sd) {
  return Rf_rnorm(current_sd, proposal_sd);
}

// Is proposal within the parameter support?
bool parameter_support(double proposal, double min = 0, double max) {
  // re_sdmax should be a const from the priors class, how to get it in here?
  return (proposal > min && proposal < max);
}


double posterior_function(Patient patient, ProposalVariance pv) {

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

