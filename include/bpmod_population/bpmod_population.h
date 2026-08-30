#ifndef GUARD_bpmod_population_h
#define GUARD_bpmod_population_h

//
// bpmod_population.h
//   Main header file for population-level Bayesian pulsatile hormone modeling
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// This file provides the infrastructure for hierarchical Bayesian modeling
// of pulsatile hormone data from multiple subjects. The population model
// estimates:
//
//   Population-level parameters (hyperparameters):
//     - μ_α, μ_ω: population means of subject mean pulse mass/width
//     - θ_b, θ_h: population mean baseline and half-life
//     - υ_α, υ_ω: subject-to-subject SD of mean mass/width
//     - σ_b, σ_h: subject-to-subject SD of baseline/half-life
//     - σ_α, σ_ω: pulse-to-pulse SD of mass/width (shared across subjects)
//     - σ²_e: model error variance
//
//   Subject-level parameters (for each subject s):
//     - μ_α,s, μ_ω,s: subject mean pulse mass/width
//     - θ_b,s, θ_h,s: subject baseline and half-life
//
//   Pulse-level parameters (for each pulse k in subject s):
//     - α_k,s, ω_k,s: individual pulse mass and width
//     - τ_k,s: pulse location (time)
//     - κ_α,k,s, κ_ω,k,s: t-distribution variance scales
//

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif

// Data structures
#include <bp_datastructures/population.h>
#include <bp_datastructures/patient.h>
#include <bp_datastructures/chains.h>

// MCMC infrastructure
#include <bp_mcmc/mh.h>
#include <bp_mcmc/proposalvariance.h>

// Population-level samplers
#include <bpmod_population/pop_draw_means.h>
#include <bpmod_population/pop_draw_baseline_halflife_means.h>
#include <bpmod_population/pop_draw_mean_sds.h>
#include <bpmod_population/pop_draw_baseline_halflife_sds.h>
#include <bpmod_population/pop_draw_pulse_sds.h>
#include <bpmod_population/pop_draw_error.h>

// Single-subject samplers (reused for subject-level and pulse-level updates)
#include <bpmod_singlesubject/birthdeath.h>
#include <bpmod_singlesubject/ss_draw_baselinehalflife.h>
#include <bpmod_singlesubject/ss_draw_fixedeffects.h>
#include <bpmod_singlesubject/ss_draw_randomeffects.h>
#include <bpmod_singlesubject/ss_draw_tvarscale.h>
#include <bpmod_singlesubject/ss_draw_locations.h>

//
// NOTE: The MCMC orchestration loop for the population model will be
// implemented in Phase 2. This header provides the building blocks:
//
// Typical MCMC iteration:
//   1. For each subject:
//      a. Birth-death process for pulse count
//      b. Update subject-level parameters (baseline, halflife, mass_mean, width_mean)
//      c. Update pulse-level parameters (mass, width, location, tvarscales)
//
//   2. Update population-level parameters:
//      a. Population means (μ_α, μ_ω, θ_b, θ_h) - Gibbs
//      b. Subject-to-subject SDs (υ_α, υ_ω, σ_b, σ_h) - MH
//      c. Pulse-to-pulse SDs (σ_α, σ_ω) - MH
//      d. Error variance (σ²_e) - Gibbs
//

#endif
