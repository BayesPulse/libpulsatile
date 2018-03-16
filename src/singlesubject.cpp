#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mh.h"
#include "patient.h"
#include "population.h"
#include "utils.h"
#include "birthdeath.h"
#include "ss_draw_fixedeffects.h"
#include "ss_draw_sdrandomeffects.h"
#include "ss_draw_sdrandomeffects_width.h"
#include "ss_draw_baselinehalflife.h"
#include "ss_draw_locations.h"
#include "ss_draw_randomeffects.h"
#include "ss_draw_randomeffects_width.h"
#include "ss_draw_tvarscale.h"
#include "ss_draw_tvarscale_width.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List singlesubject(Rcpp::NumericVector concentration,
                         Rcpp::NumericVector time,
                         Rcpp::List priors,
                         Rcpp::List proposalvars,
                         Rcpp::List startingvals,
                         int mcmc_iterations,
                         bool verbose,
                         int pv_adjust_iter,
                         int pv_adjust_max_iter,
                         double bivariate_pv_target_ratio,
                         double univariate_pv_target_ratio)
{

  // Check for valid input
  if ( !priors.inherits("bp_priors") ) stop("priors argument must be a bp_priors object");
  if ( !priors.inherits("bp_proposalvariance") ) stop("proposalvars argument must be a bp_proposalvariance object");
  if ( !priors.inherits("bp_startingvals") ) stop("startingvals argument must be a bp_startingvals object");
  if ( concentration.size != time.size ) stop("Time and concentration vectors must be the same size");

  // Create shorter names just for cleaner code appearance
  int adj_iter       = pv_adjust_iter;
  int adj_max        = pv_adjust_max_iter;
  double biv_target  = bivariate_pv_target_ratio;
  double univ_target = univariate_pv_target_ratio;


  // Create patient data object
  PatientData pdone(thistime, conc);

  //Create priors object
  PatientPriors ppsingle(priors["baseline_mean"],
                         priors["baseline_variance"],
                         priors["halflife_mean"],
                         priors["halflife_variance"],
                         priors["mass_mean"],
                         priors["mass_variance"],
                         priors["width_mean"],
                         priors["width_variance"],
                         priors["mass_mean_sd"],
                         priors["width_mean_sd"],
                         priors["error_alpha"],
                         priors["error_beta"],
                         priors["pulse_count"],
                         priors["strauss_repulsion"],
                         priors["strauss_repulsion_range"]);

  // Create estimates object (w/ starting vals)
  PatientEstimates pesingle(priors["baseline"],
                            priors["halflife"],
                            priors["errorsq"],
                            priors["mass_mean"],
                            priors["width_mean"],
                            priors["pulse_count"],
                            priors["mass_sd"],
                            priors["width_sd"]);

  // Create pointers
  PatientData * data = &pdone;
  PatientPriors * priors = &ppsingle;
  PatientEstimates * estimates = &pesingle;

  // Now take all of this and create a Patient object
  Patient pat(data, priors, estimates);
  Patient * patient = &pat;

  //--------------------------------------
  // Create sampler objects
  //--------------------------------------

  // Birth-death process
  BirthDeathProcess birth_death;

  // Modified Metropolis Hastings for fixed effects (mean mass & mean width)
  SS_DrawFixedEffects draw_fixeff_mass(proposalvars["mass_mean"], adj_iter,
                                       adj_max, univ_target, false);
  SS_DrawFixedEffects draw_fixeff_width(proposalvars["width_mean"], adj_iter,
                                        adj_max, univ_target, true);

  // Modified Metropolis Hastings for the standard deviation of the random
  // effects (sd mass & sd width) (patient level estimate)
  SS_DrawSDRandomEffects draw_sd_masses(proposalvars["mass_sd"], adj_iter,
                                        adj_max, univ_target, false);
  SS_DrawSDRandomEffects draw_sd_widths(proposalvars["width_sd"], adj_iter,
                                        adj_max, univ_target, true);

  // Bivariate Modified Metropolis Hastings for the baseline and half-life
  arma::vec bhl_pv = { proposalvars["baseline"], proposalvars["halflife"] };
  SS_DrawBaselineHalflife draw_blhl(bhl_pv, adj_iter, adj_max, biv_target);

  // Modified Metropolis Hastings for pulse locations (pulse level)
  if ( priors["location_prior_type"] == "strauss" ) {
    SS_DrawLocationsStrauss draw_locations(proposalvars["location"], adj_iter,
                                           adj_max, univ_target);
  } else {
    SS_DrawLocationsOS draw_locations(proposalvars["location"], adj_iter,
                                      adj_max, univ_target);
  }
  SS_DrawRandomEffects draw_pmasses(proposalvars["pulse_mass"], adj_iter,
                                    adj_max, univ_target, false);
  SS_DrawRandomEffects draw_pwidths(proposalvars["pulse_width"], adj_iter,
                                    adj_max, univ_target, true);
  SS_DrawTVarScale draw_tvarscale_mass(proposalvars["sdscale_pulse_mass"],
                                       adj_iter, adj_max, univ_target, false);
  SS_DrawTVarScale draw_tvarscale_width(proposalvars["sdscale_pulse_width"],
                                        adj_iter, adj_max, univ_target, true);


  // Sample MMH objects
  for (int i = 0; i < mcmc_iterations; i++) {

    birth_death.sample(patient, false);
    draw_fixeff.sample(patient, &patient->estimates->mass_mean);
    draw_fixeff_widths.sample(patient, &patient->estimates->mass_mean);
    draw_sd_pulse_masses.sample(patient, &patient->estimates->mass_sd, patient);
    draw_sd_pulse_widths.sample(patient, &patient->estimates->mass_sd, patient);
    draw_baselinehalflife.sample(patient, &patient->estimates->baseline_halflife);
    draw_pulse_locations_strauss.sample_pulses(patient);
    draw_pulse_masses.sample_pulses(patient);
    draw_pulse_widths.sample_pulses(patient);
    draw_pulse_tvarscale.sample_pulses(patient);
    draw_pulse_tvarscale_widths.sample_pulses(patient);

    print_diagnostic_output(verbose);
    //std::cout << "Iteration " << i << " Number of pulses = " << patient->pulses.size() << std::endl;
    save_sample(thin);

  }

}
