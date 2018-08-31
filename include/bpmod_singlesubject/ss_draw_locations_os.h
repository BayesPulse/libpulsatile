#ifndef GUARD_bpmod_singlesubject_draw_locations_os_h
#define GUARD_bpmod_singlesubject_draw_locations_os_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_datastructures/patient.h>
//#include <bp_mcmc/mh.h>
#include <bpmod_singlesubject/ss_draw_locations.h>
#include <bp_mcmc/utils.h>


//
// SS_DrawLocationsOS
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sampling the pulse locations with the every-3rd-order-statistic prior
//

class SS_DrawLocationsOS : public SS_DrawLocations
{

  public:

    // Constructors
    SS_DrawLocationsOS(double in_pv, int in_adjust_iter, int in_max_iter,
                     double in_target_ratio, bool in_verbose, int in_verbose_iter) :
      SS_DrawLocations(in_pv, in_adjust_iter, in_max_iter, in_target_ratio,
                       in_verbose, in_verbose_iter) { };

    // Destructor
    ~SS_DrawLocationsOS() { }

  private:

    // Test whether drawn value is within the parameter's support
    bool parameter_support(double val, Patient *patient) {
      return ((val <= patient->data.fitend) &&
              (val > patient->data.fitstart));
    }

    //
    // posterior_function()
    //   for strauss location prior mmh
    //
    double posterior_function(PulseIter *pulse,
                              double proposed_time,
                              Patient *patient) {

      // Internal variables
      PulseIter pulseiter          = *pulse;
      double previous_time         = 0;
      double next_time             = 0;
      double timediff_previous     = 0;
      double timediff_previous_new = 0;
      double timediff_next         = 0;
      double timediff_next_new     = 0;
      double prior_ratio           = 0;
      double plikelihood           = 0;
      double acceptance_ratio      = 0;
      // Extracted variables
      double current_time         = pulseiter->time;
      double curr_likelihood      = patient->likelihood(false);
      arma::vec curr_mean_contrib = pulseiter->get_mean_contribution(patient->data.time,
                                                                     patient->estimates.get_decay());

      // Start calculations
      previous_time = (pulseiter == patient->pulses.begin()) ?  
        patient->data.fitstart : std::prev(pulseiter)->time;
      next_time     = (std::next(pulseiter) == patient->pulses.end()) ?  
        patient->data.fitend : std::next(pulseiter)->time; 
      //Rcpp::Rcout << "current_time  = " << current_time   << std::endl;
      //Rcpp::Rcout << "proposed_time = " << proposed_time  << std::endl;
      //Rcpp::Rcout << "previous_time = " << previous_time  << std::endl;
      //Rcpp::Rcout << "next_time     = " << next_time      << std::endl;

      //if (pulseiter != patient->pulses.begin()) {
      //  // if we're not at the first pulse, use previous pulse locaiton
      //  previous_time = pulseiter.prev()->time;
      //} else {
      //  // Otherwise, use fitstart
      //  previous_time = patient->data.fitstart;
      //}
      //if (pulseiter != patient->pulses.end()) {
      //  // If we're not at the last pulse, use the next pulse location
      //  next_time = pulseiter.next()->time;
      //} else {
      //  // Otherwise, use fitend
      //  next_time = patient->data.fitend;
      //}

      // Take time diffs
      timediff_previous_new  = proposed_time - previous_time;
      timediff_previous      = current_time - previous_time;
      timediff_next_new      = next_time - proposed_time;
      timediff_next          = next_time - current_time;

      // Combine it all for the prior ratio
      prior_ratio = (patient->priors.num_orderstat-1) *
        log((timediff_previous_new * timediff_next_new) /
            (timediff_previous * timediff_next));

      // Save current time and set its time to proposed value, then calculate
      // the likelihood under the proposed value returns log-likelihood  --
      // calculate with new pulse->time
      current_time = pulseiter->time;
      pulseiter->time  = proposed_time;
      plikelihood  = patient->likelihood(false); 

      // Calculate the likelihood ratio 
      acceptance_ratio = prior_ratio + (plikelihood - curr_likelihood);

      // Reset pulse->time to current (sample() chooses whether to keep) and
      // get_mean_contribution() will recalc that when requested.
      pulseiter->time = current_time;

      return acceptance_ratio;

    }

};


#endif

