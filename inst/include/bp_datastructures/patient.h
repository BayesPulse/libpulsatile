#ifndef GUARD_patient_h
#define GUARD_patient_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include "bp_datastructures/patientdata.h"
#include "bp_datastructures/patientestimates.h"
#include "bp_datastructures/patientpriors.h"
#include "bp_datastructures/pulseestimates.h"



//
// patient.h
//   defining the patient class and subclasses
//
//   NOTE: getting too big -- clean up/reorganize
//


//
// Parent class of Patients and Populations, object that MMH and Gibbs
// recognize
//
//struct MCMCSamplingUnit { };


using namespace Rcpp;

typedef std::list<PulseEstimates> PulseList;
typedef PulseList::iterator PulseIter;
typedef PulseList::const_iterator PulseConstIter;


//
// Patient struct
//
struct Patient {

  // Member objects for all models
  PatientData data;
  PatientEstimates estimates;
  PulseList pulses;
  PulseIter piter = pulses.begin();
  PulseList responses;
  PulseIter riter = responses.begin();

  //
  // For single-subject model
  //

  // Member objects
  PatientPriors priors;

  // Constructor
  Patient(PatientData in_data,
          PatientPriors in_priors,
          PatientEstimates in_parms) :
      data(in_data)
    , estimates(in_parms)
    , priors(in_priors) {

    PulseEstimates firstpulse(in_data.fitstart,
                              1, 1, 1, 1,
                              in_parms.get_decay(), 
                              in_data.concentration);
    pulses.push_back(firstpulse);

  };


  //
  // For population models
  //

  // Constructor
  Patient(PatientData in_data,
          PatientEstimates in_parms) :
      data(in_data)
    , estimates(in_parms) 
    , priors() {
      // priors member obj not used

    PulseEstimates firstpulse(in_data.fitstart,
                              1, 1, 1, 1,
                              in_parms.get_decay(), 
                              in_data.concentration);
    pulses.push_back(firstpulse);
  };


  //
  // Methods for both types/all-models
  //

  // get_pulsecount()
  //   Get current number of pulses
  int get_pulsecount() { return pulses.size(); };

  // get_sumerrorsquared()
  //   Sums of squared error for error gibbs
  double get_sumerrorsquared(bool response_hormone) {

    arma::vec mean = mean_concentration(response_hormone);
    arma::vec conc;

    if (response_hormone) {
      conc = data.response_concentration;
    } else {
      conc = data.concentration;
    }

    arma::vec diff = conc - mean;
    arma::vec squareddiff =  diff % diff; // % Schur product (elementwise multiplication)
    double ssq = arma::accu(squareddiff);

    //Rcpp::Rcout << "Conc vector:\n" << conc << std::endl;
    //Rcpp::Rcout << "Mean vector:\n" << mean << std::endl;
    //Rcpp::Rcout << "Diff vector:\n" << diff << std::endl;
    //Rcpp::Rcout << "Square diff vector:\n" << squareddiff << std::endl;
    //Rcpp::Rcout << "Sums of squares:" << ssq << std::endl;

    return ssq;

  };




  // likelihood()
  //   computes the current likelihood using the observed log-concentrations and
  //   mean concentration as requested. Not stored.
  //   there is a version for a) using all pulses and b) one for excluding one
  //   pulse.
  double likelihood(bool response_hormone) {
    PulseIter emptyiter;
    return likelihood(response_hormone, emptyiter);
  };

  double likelihood(bool response_hormone, PulseIter pulse_excluded) {

    double like = 0;
    arma::vec conc;
    arma::vec diffs;

    if (response_hormone) {
      conc = data.response_concentration;
    } else {
      conc = data.concentration;
    }

    // Calculate likelihood
    diffs = conc - mean_concentration(response_hormone, pulse_excluded);
    diffs = arma::square(diffs);

    //Rcpp::Rcout << "Diffs: ";
    //for(auto diff : diffs) { Rcpp::Rcout << diff << " "; }
    //Rcpp::Rcout << "\n";
    //Rcpp::Rcout << "Mean conc: " << mean_concentration(response_hormone, pulse_excluded);

    like  = arma::accu(diffs);
    //Rcpp::Rcout << "\nLike 1: " << like;

    like /= (-2.0 * estimates.errorsq);
    //Rcpp::Rcout << "\nLike 2: " << like
    //            << "\nErrorsq: " << estimates.errorsq;

    like += -0.5 * conc.n_elem * (1.8378771 + estimates.get_logerrorsq());
    //Rcpp::Rcout << "\nLike 3: " << like
    //            << "\nLogErrorSq: " << estimates.get_logerrorsq()
    //            << "\n\n";

    return like;

  };

  // Calculate all partial likelihoods (excluding each pulse(i))
  arma::vec get_partial_likelihood(bool response_hormone) {

    PulseIter exclude_pulse = pulses.begin();;
    PulseIter pulse_end     = pulses.end();;
    arma::vec partials(get_pulsecount());

    int i = 0;
    while(exclude_pulse != pulse_end) {
      partials(i) = likelihood(response_hormone, exclude_pulse);
      exclude_pulse++;
      i++;
    }

    return partials;

  };

  // mean_concentration()
  //   this takes each pulse's mean_contrib vector and sums across them
  //   there is a version for a) using all pulses and b) one for excluding one
  //   pulse.
  arma::vec mean_concentration(bool response_hormone) {

    PulseIter emptyiter;
    return mean_concentration(response_hormone, emptyiter);

  };

  arma::vec mean_concentration(bool response_hormone, PulseIter pulse_excluded)
  {

    arma::vec mean_conc(data.concentration.n_elem);
    mean_conc.fill(0);
    PulseIter pulse_iter;
    PulseConstIter pulselist_end;

    if (response_hormone) {
      pulse_iter    = responses.begin();
      pulselist_end = responses.end();
    } else {
      pulse_iter    = pulses.begin();
      pulselist_end = pulses.end();
    }

    // Add the contribution to the mean from each pulse
    arma::vec mctrb(data.concentration.n_elem);
    mctrb.fill(0);

    while (pulse_iter != pulselist_end) {
      if (pulse_iter != pulse_excluded) {
        mctrb      = pulse_iter->get_mean_contribution(data.time, estimates.get_decay());
        mean_conc += mctrb;
      }
      ++pulse_iter;
    }

    // Add the baseline contribution and log
    if (estimates.baseline_halflife.n_elem == 2) {
      mean_conc += estimates.baseline_halflife(0);
    } else {
      mean_conc += estimates.baseline;
    }
    //Rcpp::Rcout << "Sum(mean_conc) = " << sum(mean_conc) << "\n";
    mean_conc = log(mean_conc);
    //Rcpp::Rcout << "Sum(log(mean_conc)) = " << sum(mean_conc) << "\n\n";

    return mean_conc;

  };

  // calc_sr_strauss()
  //   Calculates sum(S(R)), the exponent on the gamma parameter in the Strauss
  //   process/prior for pulse location. Used for Strauss prior in birth_death
  //   and mh_time_strauss.
  //   location is time/loc to test against other pulses
  int calc_sr_strauss(double location) {
    PulseIter emptyiter;
    return calc_sr_strauss(location, &(*emptyiter));
  };

  int calc_sr_strauss(double location, PulseEstimates * pulse_excluded) {

    int s_r = 0;       // Sum of indicators where diff < 20
    double difference = 0.; // Time difference
    PulseIter pulse = pulses.begin();
    PulseConstIter pulse_end = pulses.end();

    while (pulse != pulse_end) {
      if (&(*pulse) != pulse_excluded) { // TODO: Test that pulse is actually excluded!
        // skip if node is same that location is from;
        difference = abs(location - pulse->time);
        // increment by 1 if diff<R
        s_r = (difference < priors.strauss_repulsion_range) ? s_r + 1 : s_r; 
      }
      pulse++;
    }

    // sum(S(R)) - scalar value for sum of # pulses too close to each other
    return(s_r); 

  };



  void fix_estimates(Rcpp::List fix_params,
                     Rcpp::NumericVector masses_vec,
                     Rcpp::NumericVector width_vec,
                     Rcpp::NumericVector mass_tvarscale_vec,
                     Rcpp::NumericVector width_tvarscale_vec,
                     Rcpp::NumericVector location_vec) {
  
    double position, new_mass, new_width, new_tvarscale_mass, new_tvarscale_width,
           new_t_sd_mass, new_t_sd_width;
    int l = 0;
    int i = 0;
    int j = 0;

    if(fix_params["pulse_count"]) {

        // Get number of pulses for patient
        j = priors.pulse_count;
        Rcpp::RNGScope rng_scope;

        // Generate first pulse (must be done outside loop, otherwise push_back()
        //  adds first pulse on top of the one generated by pulse constructor)
        position = (fix_params["pulse_location"]) 
          ? location_vec(l) : Rf_runif(data.fitstart, data.fitend);
        new_tvarscale_mass = (fix_params["pulse_mass_sdscale"]) 
          ? mass_tvarscale_vec(l) : Rf_rgamma(2, 0.5);
        new_tvarscale_width = (fix_params["pulse_width_sdscale"]) 
          ? width_tvarscale_vec(l) : Rf_rgamma(2, 0.5);

        // Can fix values to known value (rather than gammas) for testing purposes
        //new_tvarscale_mass = (fix_params["mass_tvarscale"]) ? testMKappaVec(l) : 10;
        //new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : 10;
  
        if(fix_params["pulse_mass"]) {
          new_mass = masses_vec(l);
        } else {
          new_t_sd_mass = (1/sqrt(estimates.mass_prec)) / sqrt(new_tvarscale_mass);
          new_mass = -1.0;
          while (new_mass < 0) {
            new_mass = Rf_rnorm(estimates.mass_mean, new_t_sd_mass);
            //new_mass = Rf_rnorm(10, .1);
          }
        }

        if(fix_params["pulse_width"]) {
          new_width = width_vec(l);
        } else {
          new_t_sd_width = (1/sqrt(estimates.width_prec)) / sqrt(new_tvarscale_width);
          new_width = -1.0;
          while (new_width < 0) {
            new_width = Rf_rnorm(estimates.width_mean, new_t_sd_width);
            //new_width = Rf_rnorm(70, 1);
          }
        }
      
        PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                                 new_tvarscale_width, estimates.get_decay(),
                                 data.time);

        pulses.front() = new_pulse;

        l++;

        // Fix remainder of pulses using loop
        for(int k = 0; k < j-1; k++) {
          position = (fix_params["pulse_location"]) 
             ? location_vec(l) : Rf_runif(data.fitstart, data.fitend);
           new_tvarscale_mass = (fix_params["pulse_mass_sdscale"]) 
             ? mass_tvarscale_vec(l) : Rf_rgamma(2, 0.5);
           new_tvarscale_width = (fix_params["pulse_width_sdscale"]) 
             ? width_tvarscale_vec(l) : Rf_rgamma(2, 0.5);

           // Can fix values to known value (rather than gammas) for testing purposes
           //new_tvarscale_mass = (fix_params["mass_tvarscale"]) ? testMKappaVec(l) : 10;
           //new_tvarscale_width = (!test_tvarscale_width) ? testWKappaVec(l) : 10;
        
           if(fix_params["pulse_mass"]) {
             new_mass = masses_vec(l);
           } else {
             new_t_sd_mass = (1/sqrt(estimates.mass_prec)) / sqrt(new_tvarscale_mass);
             new_mass = -1.0;
             while (new_mass < 0) {
               new_mass = Rf_rnorm(estimates.mass_mean, new_t_sd_mass);
               //new_mass = Rf_rnorm(10, .1);
             }
           }

           if(fix_params["pulse_width"]) {
             new_width = width_vec(l);
           } else {
             new_t_sd_width = (1/sqrt(estimates.width_prec)) / sqrt(new_tvarscale_width);
             new_width = -1.0;
             while (new_width < 0) {
               new_width = Rf_rnorm(estimates.width_mean, new_t_sd_width);
               //new_width = Rf_rnorm(70, 1);
             }
           }
          
           PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                                    new_tvarscale_width, estimates.get_decay(),
                                    data.time);

           pulses.push_back(new_pulse);

           l++;
        }

        i++;

      Rcpp::Rcout << "Pulses Added Manually\n";

      Rcpp::Rcout << "Details:\n";
      j = l = 1;
      for(auto pulse : pulses) {
        Rcpp::Rcout << std::setprecision(3) << j << " "
                    << l << " "
                    << pulse.mass << " "
                    << pulse.width << " "
                    << pulse.time << " "
                    << pulse.tvarscale_mass << " "
                    << pulse.tvarscale_width << " " << "\n";
        l++;
      }

      Rcpp::Rcout << "\n";

    }
  
  };


};

#endif
