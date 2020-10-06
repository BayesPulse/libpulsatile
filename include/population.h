#ifndef GUARD_population_h
#define GUARD_population_h

//#include <bp_datastructures/bp_datastructures.h>
//#include <bp_mcmc/bp_mcmc.h>
//#include <bpmod_singlesubject/bpmod_singlesubject.h>

#include <Rcpp.h>
Rcpp::List population_(Rcpp::NumericMatrix concentrations,
                          Rcpp::NumericVector time,
                          Rcpp::CharacterVector location_prior,
                          Rcpp::List inpriors,
                          Rcpp::List proposalvars,
                          Rcpp::List startingvals,
                          int mcmc_iterations,
                          int thin,
                          int burnin,
                          bool verbose,
                          int pv_adjust_iter,
                          int pv_adjust_max_iter,
                          double bivariate_pv_target_ratio,
                          double univariate_pv_target_ratio,
                          bool test_birthdeath = false,
                          bool test_sd_masses = false,
                          bool test_sd_width = false,
                          bool test_s2s_sd_width = false,
                          bool test_s2s_sd_mass = false,
                          bool test_s2s_sd_baseline = false,
                          bool test_s2s_sd_halflife = false,
                          bool test_pop_means_width = false,
                          bool test_pop_means_mass = false,
                          bool test_pop_means_baseline = false,
                          bool test_pop_means_halflife = false,
                          bool test_fixeff_mass = false,
                          bool test_fixeff_width = false,
                          bool test_error = false,
                          bool test_locations = false,
                          bool test_masses = false,
                          bool test_widths = false);




#endif
