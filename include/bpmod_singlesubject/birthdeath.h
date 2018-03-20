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
    void sample(Patient *patient, bool response_hormone);
  private:
    PulseUtils pu;
    void add_new_pulse(Patient *patient, double position);
    void remove_pulse(Patient *patient, arma::vec death_rate, int pulse_count);
    double calculate_total_deathrate(arma::vec death_rate, double pulse_count);
    double calculate_total_deathrate_original(arma::vec death_rate, double pulse_count);
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

void BirthDeathProcess::sample(Patient *patient, bool response_hormone) {

  int aaa                 = 0;   // Counter for # BD iterations
  int max_num_node        = 60;  // Max number of pulses allowed before forced death (hardcoded parm)
  int pulse_count         = 0;
  double max              = 0.0; // used in two different places (double check this)
  double S                = 0.0; // Current progress in virtual time
  double T                = 1.0; // Stop birthdeath loop when S exceeds T
  double total_birth_rate = 0.0; // A constant birth rate needed for the process
  double birth_rate       = 0.0; // 'Instantaneous' birth rate
  double position         = 0.0; // Value of new pulse's location
  double total_death_rate = 0.0;
  arma::vec death_rate;       // Vector of death rates for each pulse
  double fitstart = patient->data->fitstart;
  double fitend   = patient->data->fitend;
  //// Spatial BD and Strauss declarations
  int sum_s_r = 0;            // Sum(S(R)) for birth of new pulse
  double papas_cif;           // Papangelou's cond'l intensity fn for birthed
  double b_ratio;             // Ratio of papas_cif/birth_rate for accept/rej

  // TEMP -- add as argument
  int strauss = 1;
  std::cout << "In the birth-death sample function" << std::endl;
  std::cout << "Current likelihood is " << patient->likelihood(false) << std::endl;

  //-----------------------------------------
  // Prepare for birth death loop. Save and
  // calculate values and allocate memory
  //-----------------------------------------
  total_birth_rate = (double)patient->priors->pulse_count;

  std::cout << "Total birth rate is " << total_birth_rate << std::endl;

  // If Strauss, calculate instantaneous birth rate (which is the prior intensity)
  if (strauss == 1) {
    birth_rate = total_birth_rate / (fitend - fitstart);
  }
  std::cout << "Instantaneous birth rate is " << birth_rate << std::endl;


  // Start Birth-death loop Run until break reached
  while (1) {

    aaa++; // iters counter

    // Calculate death rate for each component conditional on Strauss or order
    // statistic -- using pulse count and partial likelihoods
    pulse_count = patient->get_pulsecount();
    std::cout << "Pulse count is " << pulse_count << std::endl;

    arma::vec partial_likelihood = patient->get_partial_likelihood(response_hormone);;
    //std::cout << "Partial likelihoods are " << partial_likelihood << std::endl;

    if (strauss == 1) {
      death_rate = calc_death_rate_strauss(patient, partial_likelihood,
                                           pulse_count, response_hormone);
    } //else {
    //  death_rate = calc_death_rate_os(list, pulse_count, partial_likelihood,
    //                                  full_likelihood, total_birth_rate, 
    //                                  patient->priors->pulse_count, priors->orderstat);
    //}


    // Compute total death rate D  -- multi-step due to precision issues
    //total_death_rate = calculate_total_deathrate(death_rate, pulse_count);
    total_death_rate = calculate_total_deathrate_original(death_rate, pulse_count);

    // Calculate probability of birth
    // Set Pr(Birth), force death/birth if too many/few
    if (pulse_count <= 1) { 
      max = 1.1;
    } else if (pulse_count >= max_num_node) { //TODO: may want to remove this upper limit for production..
      max = -0.1;
    } else { 
      max = total_birth_rate / (total_birth_rate + total_death_rate); 
    }


    // Update virtual time (how long to run BD step) - Draw from exp(B+D) and add to current S 
    double lambda = 1 /( total_birth_rate + total_death_rate);
    std::cout << "total birth rate = " << total_birth_rate << " total_death_rate = " << total_death_rate << std::endl;
    std::cout << "lambda = " << lambda << std::endl;
    S += Rf_rexp(lambda);
    // If S exceeds T or if we've run this too many times, break
    if (S > T)      { 
      std::cout << "In break 1" << std::endl;
      std::cout << "S = " << S << " and T = " << T << std::endl;
      break; 
    }
    if (aaa > 5000) { 
      std::cout << "In break 2" << std::endl;
      break; 
    }


    // Select Birth or Death and proceed with either
    if (Rf_runif(0, 1) < max) { // If U < B/(B+D), a birth occurs 

      std::cout << "In possible birth step" << std::endl;

      // Generate new position
      position = Rf_runif(fitstart, fitend);
      int accept_pos = 1;

      // If using Strauss prior, run accept/reject for new position
      if (strauss == 1) {
        sum_s_r    = patient->calc_sr_strauss(position); 
        papas_cif  = birth_rate * pow(patient->priors->strauss_repulsion, sum_s_r);
        b_ratio    = papas_cif / birth_rate;

        accept_pos = (Rf_runif(0, 1) < b_ratio) ? 1 : 0;
      }

      // If it's a valid position, generate initial parms and insert node
      if (accept_pos == 1) {
        std::cout << "In for sure birth step" << std::endl;
        add_new_pulse(patient, position);
      }

    } else { // Otherwise, a death occurs 
        std::cout << "In for sure death step" << std::endl;
        remove_pulse(patient, death_rate, pulse_count);
    }

  }


};




//
// add_new_pulse()
//
void BirthDeathProcess::add_new_pulse(Patient *patient, double position) {

  double new_mass = 0.0;
  double new_width = 0.0;
  double new_tvarscale_mass, new_tvarscale_width;
  new_tvarscale_mass = new_tvarscale_width = 1.0;

  while (new_mass < 0) {
    new_mass = Rf_rnorm(patient->estimates->mass_mean, patient->estimates->mass_sd);
  }
  while (new_width < 0) {
    new_width = Rf_rnorm(patient->estimates->width_mean, patient->estimates->width_sd);
  }
  // In old version, eta not set. Trying with set to 1.0 
  //while (tvarscale_mass < 0) {
  //  new_scale_mass = Rf_rnorm(patient->estimates->width_mean, patient->estimates->width_sd);
  //}

  // Create new pulse and insert
  PulseEstimates new_pulse(position, new_mass, new_width, new_tvarscale_mass,
                          new_tvarscale_width, patient->estimates->get_decay(),
                          patient->data->time);

  patient->pulses.push_back(new_pulse);

};




//
// remove_pulse()
//
void BirthDeathProcess::remove_pulse(Patient *patient, 
                                     arma::vec death_rate, 
                                     int pulse_count) {

  PulseIter pulse = patient->pulses.begin();
  int remove = 0;

  // Pick a node to remove, find and remove it
  remove = pu.one_rmultinom(death_rate) + 1; 
  for (int i = 0; i < remove; i++)  pulse++;
  pulse = patient->pulses.erase(pulse);

};




//
// calculate_total_deathrate()
//
// new, simplified version (may have numerical problems)
double BirthDeathProcess::calculate_total_deathrate(arma::vec death_rate, double pulse_count) {

  double total_death_rate = 0.0;

  if (death_rate.size() != 0) {
    total_death_rate = arma::accu(death_rate);

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

  return total_death_rate;
};




// original version
double BirthDeathProcess::calculate_total_deathrate_original(arma::vec death_rate, double pulse_count) {

  double total_death_rate = 0.0;
  double max = 0.0;

  if (death_rate.size() != 0) {

    total_death_rate = death_rate(0);

    for (int i = 1; i < pulse_count; i++) {
      max        = (total_death_rate > death_rate(i)) ? total_death_rate : death_rate(i);
      total_death_rate = log(exp(total_death_rate - max) + exp(death_rate(i) - max)) + max;
    }

    for (int i = 0; i < pulse_count; i++) {
      death_rate(i) -= total_death_rate;
      death_rate(i)  = exp(death_rate(i));
    }

    for (int i = 1; i < pulse_count; i++) {
      death_rate(i) += death_rate(i-1);
    }

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

  return total_death_rate;

};




// 
//  calc_death_rate_strauss()
//    Calculates a vector of death rates, one for each existing pulse
//
// Notes: 
//   Two possible version:
//      1) if we can use distr of birth_rate, it's simple -- just partial -
//      full likelihoods.
//      2) if we have to use our set birth_rate (which I suspect based on
//      Stephens2000), we have to multiply by papas^-1 * birth_rate
//
//   Resolution:
//     1) is the correct approach. We aren't actually setting the birth rate.
//     Instead we are providing the intensity parameter as the strauss and as
//     long as death rate and birth rate intensity parms are the same, they
//     cancel. Option 2) would require not using an accept/reject in birth step
//     and instead generating from the actual Strauss density (i.e. really
//     difficult and computationally expensive.
// 
//-----------------------------------------------------------------------------
arma::vec BirthDeathProcess::calc_death_rate_strauss(Patient *patient,
                                                     arma::vec partial_likelihood,
                                                     int pulse_count,
                                                     bool response_hormone) {

  arma::vec death_rate(pulse_count); // Individual death rate array;

  // Begin death rate calculations
  if (pulse_count > 1) {

    // Calculate death rates (priors all cancel)
    death_rate = partial_likelihood - patient->likelihood(response_hormone);

  } else { 

    // if we have 0 or 1 pulses return NULL/0 (i.e. don't kill anything)
    if (pulse_count == 0) { 
      return NULL;
    } else {
      death_rate(0) = -1e300;
    }

  }

  return death_rate;

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
//  double *death_rate; // Vector of death rates
//  double coef_denom;  // Part of death rate calculation
//  double coef_num;    // Part of death rate calculation
//  arma::vec death_rate(pulse_count);

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
//      death_rate[i] = x;

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
//      death_rate = (double *)calloc(num_node,sizeof(double));
//      death_rate[0] = -1e300;
//    }

//  }

//  return death_rate;  

//};



#endif

