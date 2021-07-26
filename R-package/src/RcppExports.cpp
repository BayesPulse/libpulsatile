// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// population_
Rcpp::List population_(Rcpp::NumericMatrix concentrations, Rcpp::NumericVector time, Rcpp::CharacterVector location_prior, Rcpp::List inpriors, Rcpp::List proposalvars, Rcpp::List startingvals, int mcmc_iterations, int thin, int burnin, bool verbose, bool verbose_patient, int verbose_iter, int pv_adjust_iter, int pv_adjust_max_iter, double bivariate_pv_target_ratio, double univariate_pv_target_ratio, Rcpp::List fix_params, Rcpp::NumericVector pat_mass_mean_vec, Rcpp::NumericVector pat_width_mean_vec, Rcpp::NumericVector pat_baseline_vec, Rcpp::NumericVector pat_halflife_vec, Rcpp::NumericVector pulse_count_vec, Rcpp::NumericVector masses_vec, Rcpp::NumericVector width_vec, Rcpp::NumericVector mass_sdscale_vec, Rcpp::NumericVector width_sdscale_vec, Rcpp::NumericVector location_vec);
RcppExport SEXP _bayespulse_population_(SEXP concentrationsSEXP, SEXP timeSEXP, SEXP location_priorSEXP, SEXP inpriorsSEXP, SEXP proposalvarsSEXP, SEXP startingvalsSEXP, SEXP mcmc_iterationsSEXP, SEXP thinSEXP, SEXP burninSEXP, SEXP verboseSEXP, SEXP verbose_patientSEXP, SEXP verbose_iterSEXP, SEXP pv_adjust_iterSEXP, SEXP pv_adjust_max_iterSEXP, SEXP bivariate_pv_target_ratioSEXP, SEXP univariate_pv_target_ratioSEXP, SEXP fix_paramsSEXP, SEXP pat_mass_mean_vecSEXP, SEXP pat_width_mean_vecSEXP, SEXP pat_baseline_vecSEXP, SEXP pat_halflife_vecSEXP, SEXP pulse_count_vecSEXP, SEXP masses_vecSEXP, SEXP width_vecSEXP, SEXP mass_sdscale_vecSEXP, SEXP width_sdscale_vecSEXP, SEXP location_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type concentrations(concentrationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type time(timeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type location_prior(location_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type inpriors(inpriorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type proposalvars(proposalvarsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type startingvals(startingvalsSEXP);
    Rcpp::traits::input_parameter< int >::type mcmc_iterations(mcmc_iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose_patient(verbose_patientSEXP);
    Rcpp::traits::input_parameter< int >::type verbose_iter(verbose_iterSEXP);
    Rcpp::traits::input_parameter< int >::type pv_adjust_iter(pv_adjust_iterSEXP);
    Rcpp::traits::input_parameter< int >::type pv_adjust_max_iter(pv_adjust_max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type bivariate_pv_target_ratio(bivariate_pv_target_ratioSEXP);
    Rcpp::traits::input_parameter< double >::type univariate_pv_target_ratio(univariate_pv_target_ratioSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fix_params(fix_paramsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pat_mass_mean_vec(pat_mass_mean_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pat_width_mean_vec(pat_width_mean_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pat_baseline_vec(pat_baseline_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pat_halflife_vec(pat_halflife_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pulse_count_vec(pulse_count_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type masses_vec(masses_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type width_vec(width_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mass_sdscale_vec(mass_sdscale_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type width_sdscale_vec(width_sdscale_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type location_vec(location_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(population_(concentrations, time, location_prior, inpriors, proposalvars, startingvals, mcmc_iterations, thin, burnin, verbose, verbose_patient, verbose_iter, pv_adjust_iter, pv_adjust_max_iter, bivariate_pv_target_ratio, univariate_pv_target_ratio, fix_params, pat_mass_mean_vec, pat_width_mean_vec, pat_baseline_vec, pat_halflife_vec, pulse_count_vec, masses_vec, width_vec, mass_sdscale_vec, width_sdscale_vec, location_vec));
    return rcpp_result_gen;
END_RCPP
}
// singlesubject_
Rcpp::List singlesubject_(Rcpp::NumericVector concentration, Rcpp::NumericVector time, Rcpp::CharacterVector location_prior, Rcpp::List inpriors, Rcpp::List proposalvars, Rcpp::List startingvals, int mcmc_iterations, int thin, int burnin, bool verbose, int verbose_iter, int pv_adjust_iter, int pv_adjust_max_iter, double bivariate_pv_target_ratio, double univariate_pv_target_ratio, Rcpp::List fix_params, Rcpp::NumericVector masses_vec, Rcpp::NumericVector width_vec, Rcpp::NumericVector mass_tvarscale_vec, Rcpp::NumericVector width_tvarscale_vec, Rcpp::NumericVector location_vec);
RcppExport SEXP _bayespulse_singlesubject_(SEXP concentrationSEXP, SEXP timeSEXP, SEXP location_priorSEXP, SEXP inpriorsSEXP, SEXP proposalvarsSEXP, SEXP startingvalsSEXP, SEXP mcmc_iterationsSEXP, SEXP thinSEXP, SEXP burninSEXP, SEXP verboseSEXP, SEXP verbose_iterSEXP, SEXP pv_adjust_iterSEXP, SEXP pv_adjust_max_iterSEXP, SEXP bivariate_pv_target_ratioSEXP, SEXP univariate_pv_target_ratioSEXP, SEXP fix_paramsSEXP, SEXP masses_vecSEXP, SEXP width_vecSEXP, SEXP mass_tvarscale_vecSEXP, SEXP width_tvarscale_vecSEXP, SEXP location_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type concentration(concentrationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type time(timeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type location_prior(location_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type inpriors(inpriorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type proposalvars(proposalvarsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type startingvals(startingvalsSEXP);
    Rcpp::traits::input_parameter< int >::type mcmc_iterations(mcmc_iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type verbose_iter(verbose_iterSEXP);
    Rcpp::traits::input_parameter< int >::type pv_adjust_iter(pv_adjust_iterSEXP);
    Rcpp::traits::input_parameter< int >::type pv_adjust_max_iter(pv_adjust_max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type bivariate_pv_target_ratio(bivariate_pv_target_ratioSEXP);
    Rcpp::traits::input_parameter< double >::type univariate_pv_target_ratio(univariate_pv_target_ratioSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fix_params(fix_paramsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type masses_vec(masses_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type width_vec(width_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mass_tvarscale_vec(mass_tvarscale_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type width_tvarscale_vec(width_tvarscale_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type location_vec(location_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(singlesubject_(concentration, time, location_prior, inpriors, proposalvars, startingvals, mcmc_iterations, thin, burnin, verbose, verbose_iter, pv_adjust_iter, pv_adjust_max_iter, bivariate_pv_target_ratio, univariate_pv_target_ratio, fix_params, masses_vec, width_vec, mass_tvarscale_vec, width_tvarscale_vec, location_vec));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests();

static const R_CallMethodDef CallEntries[] = {
    {"_bayespulse_population_", (DL_FUNC) &_bayespulse_population_, 27},
    {"_bayespulse_singlesubject_", (DL_FUNC) &_bayespulse_singlesubject_, 21},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayespulse(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
