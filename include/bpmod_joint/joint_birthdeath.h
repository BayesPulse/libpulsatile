#ifndef GUARD_bpmod_joint_birthdeath_h
#define GUARD_bpmod_joint_birthdeath_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/bp_datastructures.h>
#include <bp_datastructures/associationparameters.h>
#include <bpmod_joint/joint_update_lambda.h>
#include <bpmod_singlesubject/birthdeath.h>

//
// joint_birthdeath.h
//   Birth-death process for joint driver-response hormone model
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// The joint birth-death process extends the standard birth-death to handle
// coupling between driver and response hormones. Key differences:
//
// 1. Driver pulses: Standard Poisson birth-death process
//    - Constant birth rate across time
//    - Death rate based on likelihood contribution
//
// 2. Response pulses: Inhomogeneous Poisson birth-death process
//    - Birth rate modulated by coupling intensity λ(t) from driver pulses
//    - λ(t) = Σ_j [ρ / √(2πν)] × exp(-(t - τ_j)² / (2ν))
//    - Death rate based on likelihood contribution
//
// After any change to driver pulses (birth, death, or location update),
// all lambda values for response pulses must be recalculated.
//

using namespace Rcpp;


//
// JointBirthDeathProcess class
//   Handles birth-death for coupled driver-response hormones
//
class JointBirthDeathProcess
{

  public:

    // Sample driver hormone pulses (standard birth-death)
    void sample_driver(Patient *driver_patient, int iter);

    // Sample response hormone pulses (birth rate modulated by lambda)
    void sample_response(Patient *driver_patient,
                        Patient *response_patient,
                        const AssociationEstimates &assoc_est,
                        int iter);

    // Exposed for unit testing the response new-pulse draw (log-normal vs
    // natural-scale branch). Public access here does not change normal sampling
    // behavior -- sample_response() drives response births internally.
    void add_new_response_pulse(Patient *driver_patient,
                               Patient *response_patient,
                               const AssociationEstimates &assoc_est,
                               double position);

  private:

    PulseUtils pu;
    BirthDeathProcess standard_bd;  // Reuse standard birth-death for driver

    void remove_pulse(Patient *patient, arma::vec death_rates, int pulse_count);

    // Calculate spatially-varying birth rate for response pulses
    double calculate_lambda_birth_rate(const std::list<PulseEstimates> &driver_pulses,
                                      const AssociationEstimates &assoc_est,
                                      double position);
};


//
// sample_driver()
//   Birth-death for driver hormone (standard process)
//
inline void JointBirthDeathProcess::sample_driver(Patient *driver_patient, int iter) {

  // Driver pulses use standard birth-death process
  standard_bd.sample(driver_patient, false, iter);

}


//
// sample_response()
//   Birth-death for response hormone with coupling to driver
//
inline void JointBirthDeathProcess::sample_response(Patient *driver_patient,
                                                    Patient *response_patient,
                                                    const AssociationEstimates &assoc_est,
                                                    int iter) {

  Rcpp::RNGScope rng_scope;

  int aaa                 = 0;   // Counter for BD iterations
  int max_num_node        = 60;  // Max pulses allowed
  int pulse_count         = 0;
  double S                = 0.0; // Virtual time
  double T                = 1.0; // Stop when S exceeds T
  double total_birth_rate = 0.0; // Total integrated birth rate
  double birth_rate       = 0.0; // Instantaneous birth rate
  double fitstart = response_patient->data.fitstart;
  double fitend   = response_patient->data.fitend;

  int strauss = 1;

  // Step 1. Calculate total birth rate
  // For response hormone, birth rate is modulated by driver pulses
  // We use the average lambda across the window as an approximation
  // More sophisticated: integrate λ(t) over [fitstart, fitend]

  // Simple approximation: sample lambda at midpoint
  double mid_time = (fitstart + fitend) / 2.0;
  double lambda_mid = 0.0;
  for (const auto& driver_pulse : driver_patient->pulses) {
    lambda_mid += assoc_est.kernel(mid_time, driver_pulse.time);
  }

  // Birth rate scaling: base rate * (1 + lambda)
  // This ensures positive birth rate even with no driver pulses
  double base_birth_rate = (double)response_patient->priors.pulse_count;
  total_birth_rate = base_birth_rate * (1.0 + lambda_mid);

  if (strauss == 1) {
    birth_rate = total_birth_rate / (fitend - fitstart);
  }

  // Birth-death loop
  do {

    aaa++;

    // Calculate death rates
    pulse_count = response_patient->get_pulsecount();
    // NOTE: Use false (not true) because response patient stores data in patient.data.concentration
    // and pulses in patient.pulses (not patient.responses). The response_hormone=true mode is for
    // single-patient models where both driver and response are in one Patient object.
    arma::vec partial_likelihood = response_patient->get_partial_likelihood(false);

    arma::vec death_rates(pulse_count);
    if (strauss == 1 && pulse_count > 1) {
      death_rates = partial_likelihood - response_patient->likelihood(false);
    } else if (pulse_count == 1) {
      death_rates(0) = -1e300;
    } else if (pulse_count == 0) {
      death_rates.reset();
    }

    // Calculate total death rate with numerical stability
    double total_death_rate = 0.0;
    double max = 0.0;

    if (death_rates.size() != 0) {

      total_death_rate = death_rates(0);

      for (int i = 1; i < pulse_count; i++) {
        max = (total_death_rate > death_rates(i)) ? total_death_rate : death_rates(i);
        total_death_rate = log(exp(total_death_rate - max) + exp(death_rates(i) - max)) + max;
      }

      for (int i = 0; i < pulse_count; i++) {
        death_rates(i) -= total_death_rate;
        death_rates(i) = exp(death_rates(i));
      }

      if (total_death_rate > 500) {
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

    // Update virtual time
    double lambda = 1 / (total_birth_rate + total_death_rate);
    S += Rf_rexp(lambda);
    if (S > T)      break;
    if (aaa > 5000) break;

    // Calculate probability of birth
    double probability_of_birth = 0.0;
    if (pulse_count <= 1) {
      probability_of_birth = 1.1;
    } else if (pulse_count >= max_num_node) {
      probability_of_birth = -0.1;
    } else {
      probability_of_birth = total_birth_rate / (total_birth_rate + total_death_rate);
    }

    // Birth or death
    if (Rf_runif(0, 1) < probability_of_birth) {

      // Generate new position
      double position = Rf_runif(fitstart, fitend);
      int accept_pos = 1;

      // Strauss accept/reject
      if (strauss == 1) {
        // Calculate lambda at proposed position
        double lambda_at_pos = calculate_lambda_birth_rate(driver_patient->pulses,
                                                          assoc_est, position);

        // Modified birth rate at this position
        double birth_rate_at_pos = birth_rate * (1.0 + lambda_at_pos);

        // Use log-scale arithmetic to avoid numerical overflow/underflow
        // in pow(strauss_repulsion, sum_s_r) for extreme values
        int sum_s_r = response_patient->calc_sr_strauss(position);
        double log_papas_cif = log(birth_rate_at_pos) +
                              sum_s_r * log(response_patient->priors.strauss_repulsion);
        double log_b_ratio = log_papas_cif - log(birth_rate);
        accept_pos = (log(Rf_runif(0, 1)) < log_b_ratio) ? 1 : 0;
      }

      if (accept_pos == 1) {
        add_new_response_pulse(driver_patient, response_patient, assoc_est, position);
      }

    } else {
      remove_pulse(response_patient, death_rates, pulse_count);
    }

  } while (true);

}


//
// add_new_response_pulse()
//   Add a new response pulse with lambda calculated from driver pulses
//
inline void JointBirthDeathProcess::add_new_response_pulse(
    Patient *driver_patient,
    Patient *response_patient,
    const AssociationEstimates &assoc_est,
    double position) {

  Rcpp::RNGScope rng_scope;

  double new_mass  = -1.;
  double new_width = -1.;
  // With Gaussian random effects the per-pulse t-scale (kappa) is fixed at 1;
  // otherwise draw it from the Gamma(2, 0.5) mixing distribution (Student-t).
  double new_tvarscale_mass  = response_patient->gaussian_random_effects ? 1.0 : Rf_rgamma(2, 0.5);
  double new_tvarscale_width = response_patient->gaussian_random_effects ? 1.0 : Rf_rgamma(2, 0.5);
  double new_t_sd_mass  = response_patient->estimates.mass_sd / sqrt(new_tvarscale_mass);
  double new_t_sd_width = response_patient->estimates.width_sd / sqrt(new_tvarscale_width);

  if (response_patient->lognormal_pulses) {
    // Log-normal parameterization (papers): log(theta) ~ N(mu, sigma^2/kappa).
    // mass_mean/width_mean are the means of the LOG values; draw on the log scale
    // and exponentiate. Positivity is automatic, so no reject loop is needed.
    new_mass  = exp(Rf_rnorm(response_patient->estimates.mass_mean, new_t_sd_mass));
    new_width = exp(Rf_rnorm(response_patient->estimates.width_mean, new_t_sd_width));
  } else {
    // Natural-scale truncated-normal parameterization (research option): reject
    // draws until positive.
    while (new_mass < 0) {
      new_mass = Rf_rnorm(response_patient->estimates.mass_mean, new_t_sd_mass);
    }
    while (new_width < 0) {
      new_width = Rf_rnorm(response_patient->estimates.width_mean, new_t_sd_width);
    }
  }

  // Calculate lambda for this pulse from driver pulses
  double new_lambda = calculate_lambda_birth_rate(driver_patient->pulses,
                                                  assoc_est, position);

  // Create new pulse with lambda
  PulseEstimates new_pulse(position, new_mass, new_width,
                          new_tvarscale_mass, new_tvarscale_width,
                          response_patient->estimates.get_decay(),
                          response_patient->data.time,
                          new_lambda);

  response_patient->pulses.push_back(new_pulse);

}


//
// remove_pulse()
//   Remove a pulse from the patient
//
inline void JointBirthDeathProcess::remove_pulse(Patient *patient,
                                                arma::vec death_rates,
                                                int pulse_count) {

  PulseIter pulse = patient->pulses.begin();
  int remove = pu.one_rmultinom(death_rates);
  for (int i = 0; i < remove; i++) pulse++;
  pulse = patient->pulses.erase(pulse);

}


//
// calculate_lambda_birth_rate()
//   Calculate lambda(t) at a specific position from all driver pulses
//
inline double JointBirthDeathProcess::calculate_lambda_birth_rate(
    const std::list<PulseEstimates> &driver_pulses,
    const AssociationEstimates &assoc_est,
    double position) {

  double lambda = 0.0;

  for (const auto& driver_pulse : driver_pulses) {
    lambda += assoc_est.kernel(position, driver_pulse.time);
  }

  return lambda;
}


#endif
