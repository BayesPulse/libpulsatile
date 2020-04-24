#ifndef GUARD_bpmod_singlesubject_birthdeath_h
#define GUARD_bpmod_singlesubject_birthdeath_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/bp_datastructures.h>



//
// birthdeath.h
//   defining the birthdeath class 
//


class BirthDeathProcess
{

  public:
    void sample(Patient *patient, bool response_hormone, int iter);
  private:
    PulseUtils pu;
    void add_new_pulse(Patient *patient, double position);
    void remove_pulse(Patient *patient, arma::vec death_rates, int pulse_count);
    double calculate_total_deathrate(arma::vec death_rates, double pulse_count);
    double calculate_total_deathrate_original(arma::vec death_rates, double pulse_count);
    arma::vec calc_death_rate_strauss(Patient *patient,
                                      arma::vec partial_likelihood,
                                      int pulse_count,
                                      bool response_hormone);
    // arma::vec calc_death_rate_os(Patient *patient,
    //                              double *partial_likelihood,
    //                              double total_death_rate,
    //                              int pulse_count,
    //                              response_hormone);


};



//
// sample()
//
void BirthDeathProcess::sample(Patient *patient, bool response_hormone, int iter) {

  Rcpp::RNGScope rng_scope;

  int aaa                 = 0;   // Counter for # BD iterations
  int max_num_node        = 60;  // Max number of pulses allowed before forced death (hardcoded parm)
  int pulse_count         = 0;
  double S                = 0.0; // Current progress in virtual time
  double T                = 1.0; // Stop birthdeath loop when S exceeds T
  double total_birth_rate = 0.0; // A constant birth rate needed for the process
  double birth_rate       = 0.0; // 'Instantaneous' birth rate
  double fitstart = patient->data.fitstart;
  double fitend   = patient->data.fitend;

  // TODO: Update this ---!!!!--- Add loc_prior argument to sample();
  int strauss = 1;


  // Step 1. Calculate total birth rate (and birth rate/min for strauss)
  total_birth_rate = (double)patient->priors.pulse_count;
  // If Strauss, calculate instantaneous birth rate (which is the prior intensity)
  if (strauss == 1) {
    birth_rate = total_birth_rate / (fitend - fitstart);
  }


  // Start Birth-death loop Run until break reached
  do {

    aaa++; // iters counter

    // 2. Calculate death rate for each component conditional on Strauss or order
    //    statistic -- using pulse count and partial likelihoods
    pulse_count = patient->get_pulsecount();
    arma::vec partial_likelihood = patient->get_partial_likelihood(response_hormone);

   // std::cout << "Partial_likelihood length: " << partial_likelihood.n_elem << "\n";
   // std::cout << "First element: " << partial_likelihood[0] << "\n";

    // 3. Calculate individual death rate for each pulse
    arma::vec death_rates(pulse_count);
    if (strauss == 1) {
      death_rates = calc_death_rate_strauss(patient, partial_likelihood,
                                           pulse_count, response_hormone);
    } //else {
    //  death_rates = calc_death_rate_os(list, pulse_count, partial_likelihood,
    //                                  full_likelihood, total_birth_rate, 
    //                                  patient->priors->pulse_count, priors->orderstat);
    //}
    


    // 4. Calculate Total death rate and convert death_rates to probabilities
    // (or cumulative probabilities)
    // can't do in separate functions/steps due to precision issues (need
    // full precision of total death rate to get probbilits to sum to 1.
    //total_death_rate = calculate_total_deathrate(death_rates, pulse_count);
    //total_death_rate = calculate_total_deathrate_original(death_rates, pulse_count);
    double total_death_rate = 0.0;
    double max = 0.0;

    if (death_rates.size() != 0) {

      total_death_rate = death_rates(0);

      for (int i = 1; i < pulse_count; i++) {
        max        = (total_death_rate > death_rates(i)) ? total_death_rate : death_rates(i);
        total_death_rate = log(exp(total_death_rate - max) + exp(death_rates(i) - max)) + max;
      }

      for (int i = 0; i < pulse_count; i++) {
        death_rates(i) -= total_death_rate;
        death_rates(i)  = exp(death_rates(i));
      }

      // Convert probabilities to cumulative probabilities -- not necessary
      // anymore for one_rmultinom()
      //for (int i = 1; i < pulse_count; i++) {
      //  death_rates(i) += death_rates(i-1);
      //}

      if (total_death_rate > 500 ) {
        total_death_rate = 1e300;
      } else {
        total_death_rate = exp(total_death_rate); 
      }

    } else { 
      total_death_rate = 0; 
    }

    if (pulse_count <= 1) { 
      total_death_rate = 0; 
    }


    // 5. Update virtual time (how long to run BD step) - Draw from exp(B+D) and add to current S
    //    If S exceeds T or if we've run this too many times, break
    double lambda = 1 / ( total_birth_rate + total_death_rate );
    S += Rf_rexp(lambda);
    if (S > T)      break;
    if (aaa > 5000) break;


    // 6. Calculate probability of birth - force death/birth if too many/few
    double probability_of_birth = 0.0;
    if (pulse_count <= 1) { 
      probability_of_birth = 1.1;
    } else if (pulse_count >= max_num_node) { //TODO: may want to remove this upper limit for production..
      probability_of_birth = -0.1;
    } else { 
      probability_of_birth = total_birth_rate / (total_birth_rate + total_death_rate); 
    }


    // 7. Select birth if random uniform is small than the probability of birth
    // (and for strauss, if position drawn is valid), otherwise choose a pulse to kill.
    if (Rf_runif(0, 1) < probability_of_birth) { // If U < B/(B+D), a birth occurs 

      // Generate new position
      double position = Rf_runif(fitstart, fitend);
      int accept_pos  = 1;

      // If using Strauss prior, run accept/reject for new position by finding
      // the sum of the pulse locations too close to the new position,
      // calculating Papangelou's conditional intensity function and dividing by
      // the birth rate
      if (strauss == 1) {
        int sum_s_r      = patient->calc_sr_strauss(position);
        double papas_cif = birth_rate * pow(patient->priors.strauss_repulsion, sum_s_r);
        double b_ratio   = papas_cif / birth_rate;
        accept_pos       = (Rf_runif(0, 1) < b_ratio) ? 1 : 0;
      }

      // If it's a valid position, generate initial parms and insert node
      if (accept_pos == 1) add_new_pulse(patient, position);

    } else { 
      // Choose and remove a pulse
      remove_pulse(patient, death_rates, pulse_count);
    }
    
    //std::cout << "Birth Rate: " << birth_rate << "\n";
    //std::cout << "Death rate size: " << death_rates.size() << "\n";
    //std::cout << "Probability of Birth: " << probability_of_birth << "\n";
    //std::cout << "Pulse count: " << patient->get_pulsecount() << "\n\n";
    

  } while (true);


};




//-----------------------------------------------------------------------------
// add_new_pulse()
//-----------------------------------------------------------------------------
void BirthDeathProcess::add_new_pulse(Patient *patient, double position) {

  Rcpp::RNGScope rng_scope;
  double new_mass  = -1.;
  double new_width = -1.;
  double new_tvarscale_mass  = Rf_rgamma(2, 0.5);
  double new_tvarscale_width = Rf_rgamma(2, 0.5);
  double new_t_sd_mass  = patient->estimates.mass_sd / sqrt(new_tvarscale_mass);
  double new_t_sd_width = patient->estimates.width_sd / sqrt(new_tvarscale_width);

  while (new_mass < 0) {
    new_mass = Rf_rnorm(patient->estimates.mass_mean, new_t_sd_mass);
  }
  while (new_width < 0) {
    new_width = Rf_rnorm(patient->estimates.width_mean, new_t_sd_width);
  }

  // Create new pulse and insert
  PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                           new_tvarscale_width, patient->estimates.get_decay(),
                           patient->data.time);

  patient->pulses.push_back(new_pulse);


};




//-----------------------------------------------------------------------------
// remove_pulse()
//   Pick a node to remove, find and remove it
//-----------------------------------------------------------------------------
void BirthDeathProcess::remove_pulse(Patient *patient, 
                                     arma::vec death_rates, 
                                     int pulse_count) {

  PulseIter pulse = patient->pulses.begin();
  int remove      = pu.one_rmultinom(death_rates);
  for (int i = 0; i < remove; i++)  pulse++;
  pulse           = patient->pulses.erase(pulse);

};




//
// calculate_total_deathrate()
//
// new, simplified version (may have numerical problems)
//double BirthDeathProcess::calculate_total_deathrate(arma::vec death_rates,
//                                                    double pulse_count) {
//
//  double total_death_rate = 0.0;
//
//  if (death_rates.size() != 0) {
//    total_death_rate = arma::accu(death_rates);
//
//    if (total_death_rate > 500 ) {
//      total_death_rate = 1e300;
//    } else {
//      total_death_rate = exp(total_death_rate); 
//    }
//
//  } else { 
//    total_death_rate = 0; 
//  }
//
//  if (pulse_count <= 1) { 
//    total_death_rate = 0; 
//  }
//
//  return total_death_rate;
//};
//
//
//
//
//// original version
//double BirthDeathProcess::calculate_total_deathrate_original(arma::vec
//                                                             death_rates, double
//                                                             pulse_count) {
//
//  double total_death_rate = 0.0;
//  double max = 0.0;
//
//  if (death_rates.size() != 0) {
//
//    total_death_rate = death_rates(0);
//
//    for (int i = 1; i < pulse_count; i++) {
//      max        = (total_death_rate > death_rates(i)) ? total_death_rate : death_rates(i);
//      total_death_rate = log(exp(total_death_rate - max) + exp(death_rates(i) - max)) + max;
//    }
//
//    for (int i = 0; i < pulse_count; i++) {
//      death_rates(i) -= total_death_rate;
//      death_rates(i)  = exp(death_rates(i));
//    }
//
//    for (int i = 1; i < pulse_count; i++) {
//      death_rates(i) += death_rates(i-1);
//    }
//
//    if (total_death_rate > 500 ) {
//      total_death_rate = 1e300;
//    } else {
//      total_death_rate = exp(total_death_rate); 
//    }
//
//  } else { 
//    total_death_rate = 0; 
//  }
//
//  if (pulse_count <= 1) { 
//    total_death_rate = 0; 
//  }
//
//  return total_death_rate;
//
//};




//-----------------------------------------------------------------------------
//  calc_death_rates_strauss()
//    Calculates a vector of death rates, one for each existing pulse.
//-----------------------------------------------------------------------------
arma::vec BirthDeathProcess::calc_death_rate_strauss(Patient *patient,
                                                     arma::vec partial_likelihood,
                                                     int pulse_count,
                                                     bool response_hormone) {

  arma::vec death_rates(pulse_count); // Individual death rate array;

  // Begin death rate calculations
  if (pulse_count > 1) {

    // Calculate death rates (priors all cancel)
    death_rates = partial_likelihood - patient->likelihood(response_hormone);

  } else { 

    // if we have 0 or 1 pulses return NULL/0 (i.e. don't kill anything)
    if (pulse_count == 0) { 
      return NULL;
    } else {
      death_rates(0) = -1e300;
    }

  }

  return death_rates;

};




//
// calc_death_rate_os() 
//   Calculates a vector of death rates, one for each existing pulse
//
//arma::vec BirthDeathProces::calc_death_rate_os(Patient *patient,
//                             double *partial_likelihood, 
//                             double total_death_rate,
//                             int pulse_count, 
//                             response_hormone) {

//  int i;              // Generic counter
//  int j;              // Generic counter
//  double x;           // Individual death rate
//  double *death_rates; // Vector of death rates
//  double coef_denom;  // Part of death rate calculation
//  double coef_num;    // Part of death rate calculation
//  arma::vec death_rates(pulse_count);

//  if (pulse_count > 1) {
//    node = list->succ;
//    i = 0;

//    // Calculate the coefficient of the distribution of the taus conditional
//    // on number of pulses. In this portion, I have an extra num_node in the
//    // numerator 
//    coef_num = 1;
//    for (j = 1; j < mmm; j++) { coef_num *= j; }
//    coef_denom = mmm * num_node;
//    for (j = 1; j < mmm; j++) { coef_denom *= (mmm * num_node + j); }

//    while (node != NULL) {

//      if (mmm > 0) {
//        // This computes the portion of death rate due to the Poisson prior on
//        // N, the birth rate, and the likelihood */
//        // num_node = number of pulses
//        // r = prior on pulse count
//        x = log(num_node * total_birth_rate / r) + 
//          partial_likelihood[i] - full_likelihood;
//      } else {
//        x = log(total_birth_rate / num_node) + 
//          partial_likelihood[i] - full_likelihood;
//      }

//      // Now, compute the portion of the death rate due to distribution of the
//      // taus conditional on number of pulses.
//      if (mmm > 1) {

//        if (i ==  0) {

//          // If we are on the first pulse
//          x += (mmm-1) * log(fitend - fitstart) +
//            (mmm-1) * log((node->succ->time - fitstart) /
//                          ((node->succ->time - node->time) * 
//                           (node->time - fitstart))) +
//            log(coef_num / coef_denom);

//        } else if (i > 0) {

//          // If we are not on the first pulse
//          if (node->succ) {
//            x += (mmm-1) * log(fitend - fitstart) +
//              (mmm-1) * log((node->succ->time - node->pred->time) /
//                            ((node->succ->time - node->time) * 
//                             (node->time - node->pred->time))) +
//              log(coef_num / coef_denom);

//          } else { // If we are not on the first or last pulse 
//            x += (mmm-1) * log(fitend - fitstart) + 
//              (mmm-1) * log((fitend - node->pred->time) / 
//                            ((fitend - node->time) * 
//                             (node->time - node->pred->time))) +
//              log(coef_num / coef_denom);
//          }

//        }

//      }

//      // if pulse equal to NaN, then set to large value/guarantee death
//      // necessary when starting value is < fitstart (orderstat part of
//      // death calc causes this due to neg. value in log) 
//      if (isnan(x)) {
//        x = 1e300;
//      }
//      // Save to death rate vector 
//      death_rates[i] = x;

//      // Advance to next pulse
//      node = node->succ;

//      // Advance counter 1
//      i++;

//    }

//  } else {

//    // If we have 0/1 pulses return NULL/0
//    if (num_node == 0) {
//      return NULL;
//    } else {
//      death_rates = (double *)calloc(num_node,sizeof(double));
//      death_rates[0] = -1e300;
//    }

//  }

//  return death_rates;  

//};



#endif

