#ifndef GUARD_bpmod_joint_draw_association_h
#define GUARD_bpmod_joint_draw_association_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>
#include <bp_datastructures/associationparameters.h>

//
// joint_draw_association.h
//   Metropolis-Hastings samplers for association parameters (ρ and ν)
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// These samplers update the coupling parameters between driver and response hormones:
//   - ρ (rho): Cluster size parameter (driver effect strength)
//   - ν (nu): Cluster width parameter (temporal spread of effect)
//
// Both parameters are sampled on the log scale with normal priors:
//   - log(ρ) ~ N(μ_ρ, σ²_ρ)
//   - log(ν) ~ N(μ_ν, σ²_ν)
//
// The likelihood depends on how well the coupling intensity λ(t) from driver
// pulses explains the observed response pulse pattern.
//

using namespace Rcpp;


//
// JointData structure
//   Container for data needed by association parameter samplers
//
struct JointData {
  Patient *driver_patient;
  Patient *response_patient;
  AssociationPriors *assoc_priors;

  JointData(Patient *driver, Patient *response, AssociationPriors *priors)
    : driver_patient(driver), response_patient(response), assoc_priors(priors) {}
};


//
// Joint_DrawRho
//   MH sampler for cluster size parameter ρ
//
class Joint_DrawRho :
  public ModifiedMetropolisHastings<JointData, AssociationEstimates, double, ProposalVariance>
{

  public:

    Joint_DrawRho(double in_pv,
                  int in_adjust_iter,
                  int in_max_iter,
                  double in_target_ratio,
                  bool verbose,
                  int verbose_iter) :
      ModifiedMetropolisHastings<JointData, AssociationEstimates, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {}

  private:

    std::string get_parameter_name() { return "Association cluster size (ρ)"; }

    bool parameter_support(double val, AssociationEstimates *assoc_est);
    double posterior_function(JointData *data, double proposal, AssociationEstimates *assoc_est);

};


//
// Joint_DrawNu
//   MH sampler for cluster width parameter ν
//
class Joint_DrawNu :
  public ModifiedMetropolisHastings<JointData, AssociationEstimates, double, ProposalVariance>
{

  public:

    Joint_DrawNu(double in_pv,
                 int in_adjust_iter,
                 int in_max_iter,
                 double in_target_ratio,
                 bool verbose,
                 int verbose_iter) :
      ModifiedMetropolisHastings<JointData, AssociationEstimates, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {}

  private:

    std::string get_parameter_name() { return "Association cluster width (ν)"; }

    bool parameter_support(double val, AssociationEstimates *assoc_est);
    double posterior_function(JointData *data, double proposal, AssociationEstimates *assoc_est);

};


//------------------------------------------------------------
// Defined functions for Joint_DrawRho
//------------------------------------------------------------

// parameter_support()
//   ρ must be positive (support on log scale is all reals)
inline bool Joint_DrawRho::parameter_support(double val, AssociationEstimates *assoc_est) {
  return (val > 0.0);
}


// posterior_function()
//   Calculate log posterior ratio for ρ
//
//   Log posterior:
//     log π(ρ | data) ∝ log π(data | ρ) + log π(ρ)
//
//   Prior on log scale:
//     log(ρ) ~ N(μ_ρ, σ²_ρ)
//     => log π(ρ) = log π(log(ρ)) - log(ρ)  [Jacobian adjustment]
//
//   Likelihood contribution:
//     For birth-death process with intensity λ(t), the log-likelihood
//     of observed response pulses at times {t_i} over [a,b] is:
//
//     log L = Σ_i log(λ(t_i)) - ∫_a^b λ(t) dt
//
//     where λ(t) = base_rate × (1 + Σ_j kernel(t, τ_j))
//     and kernel depends on ρ and ν
//
inline double Joint_DrawRho::posterior_function(JointData *data,
                                                double proposal,
                                                AssociationEstimates *assoc_est) {

  double current_rho = assoc_est->get_rho();
  double current_log_rho = assoc_est->get_log_rho();
  double proposal_log_rho = log(proposal);

  double nu = assoc_est->get_nu();

  // Prior ratio on log scale (normal prior on log(ρ))
  double prior_mean = data->assoc_priors->log_rho_mean;
  double prior_var = data->assoc_priors->log_rho_var;

  double log_prior_ratio = -0.5 * ((proposal_log_rho - prior_mean) *
                                   (proposal_log_rho - prior_mean) -
                                   (current_log_rho - prior_mean) *
                                   (current_log_rho - prior_mean)) / prior_var;

  // Jacobian adjustment: sampling ρ with prior on log(ρ)
  // J = |d(log ρ)/dρ| = 1/ρ
  // log(J_proposal / J_current) = log(current_rho / proposal)
  double log_jacobian = log(current_rho) - log(proposal);

  // Likelihood ratio
  // For computational efficiency, we use a simplified approximation:
  // The likelihood depends on λ values at response pulse locations
  // and the integral of λ over the observation window

  double log_likelihood_ratio = 0.0;

  // Sum log(λ) at each response pulse location
  for (const auto& resp_pulse : data->response_patient->pulses) {
    double lambda_current = 0.0;
    double lambda_proposal = 0.0;

    for (const auto& driver_pulse : data->driver_patient->pulses) {
      double diff = resp_pulse.time - driver_pulse.time;
      double exp_term = exp(-diff * diff / (2 * nu));

      // Current: rho / sqrt(2πν) × exp(...)
      lambda_current += (current_rho / sqrt(2 * M_PI * nu)) * exp_term;

      // Proposal: rho* / sqrt(2πν) × exp(...)
      lambda_proposal += (proposal / sqrt(2 * M_PI * nu)) * exp_term;
    }

    // Add base rate (1 + lambda contribution)
    lambda_current = 1.0 + lambda_current;
    lambda_proposal = 1.0 + lambda_proposal;

    log_likelihood_ratio += log(lambda_proposal) - log(lambda_current);
  }

  // Integral term: -∫ λ(t) dt for the response-pulse Poisson intensity
  // μ(t) = c * (1 + λ(t)), where c is the base per-minute response rate. Each
  // coupling kernel (rho / sqrt(2*pi*nu)) * exp(-(t-tau)^2 / (2*nu)) integrates
  // to exactly rho over the window (independent of nu and of the window length),
  // so ∫μ dt = c * (T + n_driver * rho). The c*T part does not depend on rho and
  // cancels in the ratio, leaving a rho-dependent penalty of c * n_driver * rho.
  // c is sourced from the response prior pulse count / window, matching the
  // base birth rate used by the birth-death process (joint_birthdeath.h).
  // NOTE: earlier versions used a coefficient of (time_span * n_driver) and then
  // (n_driver) with an implicit c=1, ~100-1000x too large, which crushed rho
  // toward 0 regardless of the data.
  int n_driver = data->driver_patient->get_pulsecount();
  double T = data->response_patient->data.fitend - data->response_patient->data.fitstart;
  double c = (double)data->response_patient->priors.pulse_count / T;
  double integral_current  = c * n_driver * current_rho;
  double integral_proposal = c * n_driver * proposal;

  log_likelihood_ratio -= (integral_proposal - integral_current);

  // Total log posterior ratio
  return log_prior_ratio + log_jacobian + log_likelihood_ratio;
}


//------------------------------------------------------------
// Defined functions for Joint_DrawNu
//------------------------------------------------------------

// parameter_support()
//   ν must be positive
inline bool Joint_DrawNu::parameter_support(double val, AssociationEstimates *assoc_est) {
  return (val > 0.0);
}


// posterior_function()
//   Calculate log posterior ratio for ν
//
inline double Joint_DrawNu::posterior_function(JointData *data,
                                               double proposal,
                                               AssociationEstimates *assoc_est) {

  double current_nu = assoc_est->get_nu();
  double current_log_nu = assoc_est->get_log_nu();
  double proposal_log_nu = log(proposal);

  double rho = assoc_est->get_rho();

  // Prior ratio
  double prior_mean = data->assoc_priors->log_nu_mean;
  double prior_var = data->assoc_priors->log_nu_var;

  double log_prior_ratio = -0.5 * ((proposal_log_nu - prior_mean) *
                                   (proposal_log_nu - prior_mean) -
                                   (current_log_nu - prior_mean) *
                                   (current_log_nu - prior_mean)) / prior_var;

  // Jacobian adjustment
  double log_jacobian = log(current_nu) - log(proposal);

  // Likelihood ratio
  double log_likelihood_ratio = 0.0;

  for (const auto& resp_pulse : data->response_patient->pulses) {
    double lambda_current = 0.0;
    double lambda_proposal = 0.0;

    for (const auto& driver_pulse : data->driver_patient->pulses) {
      double diff = resp_pulse.time - driver_pulse.time;

      // Current: rho / sqrt(2πν) × exp(-diff²/(2ν))
      lambda_current += (rho / sqrt(2 * M_PI * current_nu)) *
                       exp(-diff * diff / (2 * current_nu));

      // Proposal: rho / sqrt(2πν*) × exp(-diff²/(2ν*))
      lambda_proposal += (rho / sqrt(2 * M_PI * proposal)) *
                        exp(-diff * diff / (2 * proposal));
    }

    lambda_current = 1.0 + lambda_current;
    lambda_proposal = 1.0 + lambda_proposal;

    log_likelihood_ratio += log(lambda_proposal) - log(lambda_current);
  }

  // Integral term: -∫ λ(t) dt. To leading order ∫λ(t)dt = time_span + n_driver*rho,
  // which is INDEPENDENT of nu: each Gaussian kernel integrates to rho regardless
  // of its width because the 1/sqrt(2*pi*nu) normalization exactly cancels the
  // spread. So nu enters the likelihood only through the sum-of-log-lambda term
  // above and contributes no integral term.
  // NOTE: a previous version subtracted time_span * n_driver * rho * sqrt(nu),
  // which carried both the spurious time_span factor and an incorrect sqrt(nu)
  // dependence, collapsing the nu posterior.

  return log_prior_ratio + log_jacobian + log_likelihood_ratio;
}


#endif
