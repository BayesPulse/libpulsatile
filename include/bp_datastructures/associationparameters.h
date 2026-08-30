#ifndef GUARD_associationparameters_h
#define GUARD_associationparameters_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif

//
// associationparameters.h
//   Data structures for driver-response hormone coupling parameters
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// The association between driver (e.g., LH) and response (e.g., FSH) hormones
// is modeled through a kernel function that determines how driver pulses
// affect the intensity field for response pulses:
//
//   λ(t) = Σ_j [ρ / √(2πν)] × exp(-(t - τ_j)² / (2ν))
//
// where:
//   - λ(t) is the coupling intensity at time t
//   - τ_j are driver pulse locations
//   - ρ is the cluster size parameter (driver effect strength)
//   - ν is the cluster width parameter (temporal spread of effect)
//
// Both ρ and ν are sampled on the log scale for numerical stability.
//

using namespace Rcpp;


//
// AssociationEstimates structure
//
// Holds current estimates of coupling parameters between driver and response hormones.
// These parameters are updated via MCMC and affect the birth rate of response pulses.
//
struct AssociationEstimates {

  double rho;       // Cluster size parameter (on natural scale)
  double nu;        // Cluster width parameter (on natural scale, is variance)
  double log_rho;   // Log of rho (used for MH sampling)
  double log_nu;    // Log of nu (used for MH sampling)

  // Default constructor
  AssociationEstimates()
    : rho(1.0), nu(1.0), log_rho(0.0), log_nu(0.0) { }

  // Constructor with starting values
  AssociationEstimates(double sv_rho, double sv_nu) {
    set_rho(sv_rho);
    set_nu(sv_nu);
  }

  // Setters that maintain consistency between natural and log scales
  void set_rho(double value) {
    rho = value;
    log_rho = log(value);
  }

  void set_nu(double value) {
    nu = value;
    log_nu = log(value);
  }

  void set_log_rho(double value) {
    log_rho = value;
    rho = exp(value);
  }

  void set_log_nu(double value) {
    log_nu = value;
    nu = exp(value);
  }

  // Getters
  double get_rho() const { return rho; }
  double get_nu() const { return nu; }
  double get_log_rho() const { return log_rho; }
  double get_log_nu() const { return log_nu; }

  // Kernel function: Calculate contribution of driver pulse at location tau_j
  // to coupling intensity at time t
  double kernel(double t, double tau_j) const {
    double diff = t - tau_j;
    return (rho / sqrt(2 * M_PI * nu)) * exp(-diff * diff / (2 * nu));
  }
};


//
// AssociationPriors structure
//
// Holds prior distributions for coupling parameters.
// Priors are specified on the log scale (normal priors).
//
struct AssociationPriors {

  // Priors on log(ρ) - normal distribution
  double log_rho_mean;
  double log_rho_var;

  // Priors on log(ν) - normal distribution
  double log_nu_mean;
  double log_nu_var;

  // Default constructor
  AssociationPriors() { }

  // Full constructor
  AssociationPriors(double prior_log_rho_mean,
                   double prior_log_rho_var,
                   double prior_log_nu_mean,
                   double prior_log_nu_var) {
    log_rho_mean = prior_log_rho_mean;
    log_rho_var  = prior_log_rho_var;
    log_nu_mean  = prior_log_nu_mean;
    log_nu_var   = prior_log_nu_var;
  }

  // Constructor from natural scale values (for convenience)
  // Converts to log scale priors
  static AssociationPriors from_natural_scale(
      double rho_mean, double rho_sd,
      double nu_mean, double nu_sd) {

    // Log-normal approximation: if X ~ LogNormal(μ, σ²)
    // then log(X) ~ Normal(μ, σ²)
    // E[X] ≈ exp(μ + σ²/2), Var(X) ≈ exp(2μ + σ²)(exp(σ²) - 1)

    // For rho
    double rho_cv = rho_sd / rho_mean;  // coefficient of variation
    double log_rho_var = log(1 + rho_cv * rho_cv);
    double log_rho_mean = log(rho_mean) - log_rho_var / 2;

    // For nu
    double nu_cv = nu_sd / nu_mean;
    double log_nu_var = log(1 + nu_cv * nu_cv);
    double log_nu_mean = log(nu_mean) - log_nu_var / 2;

    return AssociationPriors(log_rho_mean, log_rho_var,
                            log_nu_mean, log_nu_var);
  }
};


#endif
