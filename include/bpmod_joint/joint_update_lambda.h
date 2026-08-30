#ifndef GUARD_bpmod_joint_update_lambda_h
#define GUARD_bpmod_joint_update_lambda_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/pulseestimates.h>
#include <bp_datastructures/associationparameters.h>

//
// joint_update_lambda.h
//   Updates lambda (coupling intensity) values for response hormone pulses
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// Lambda (λ) represents the coupling intensity field that determines how
// driver hormone pulses affect the birth rate of response hormone pulses.
// For each response pulse at time t, lambda is calculated as:
//
//   λ(t) = Σ_j [ρ / √(2πν)] × exp(-(t - τ_j)² / (2ν))
//
// where the sum is over all driver pulses at locations τ_j, and ρ and ν
// are the association parameters (cluster size and width).
//
// This function must be called whenever:
//   - Driver pulses are born or die (birth-death process)
//   - Driver pulse locations are updated (MH sampling)
//   - Association parameters (ρ or ν) are updated
//

using namespace Rcpp;


//
// update_lambda_values()
//   Recalculate lambda for all response pulses based on current driver pulses
//
// Parameters:
//   driver_pulses - List of driver hormone pulse estimates
//   response_pulses - List of response hormone pulse estimates (modified in place)
//   assoc_est - Current association parameter estimates (ρ and ν)
//
inline void update_lambda_values(
    const std::list<PulseEstimates> &driver_pulses,
    std::list<PulseEstimates> &response_pulses,
    const AssociationEstimates &assoc_est) {

  // For each response pulse, calculate lambda as sum of driver contributions
  for (auto& response_pulse : response_pulses) {

    double lambda_sum = 0.0;
    double response_time = response_pulse.time;

    // Sum contributions from all driver pulses
    for (const auto& driver_pulse : driver_pulses) {
      double driver_time = driver_pulse.time;
      lambda_sum += assoc_est.kernel(response_time, driver_time);
    }

    // Update lambda field in response pulse
    response_pulse.lambda = lambda_sum;
  }
}


//
// update_single_lambda()
//   Recalculate lambda for a single response pulse
//
// This is more efficient when only one response pulse needs updating
// (e.g., during birth-death or location updates for response pulses)
//
// Parameters:
//   driver_pulses - List of driver hormone pulse estimates
//   response_pulse - Single response pulse to update (modified in place)
//   assoc_est - Current association parameter estimates
//
inline void update_single_lambda(
    const std::list<PulseEstimates> &driver_pulses,
    PulseEstimates &response_pulse,
    const AssociationEstimates &assoc_est) {

  double lambda_sum = 0.0;
  double response_time = response_pulse.time;

  // Sum contributions from all driver pulses
  for (const auto& driver_pulse : driver_pulses) {
    double driver_time = driver_pulse.time;
    lambda_sum += assoc_est.kernel(response_time, driver_time);
  }

  // Update lambda field
  response_pulse.lambda = lambda_sum;
}


#endif
