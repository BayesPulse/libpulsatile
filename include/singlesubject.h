#ifndef GUARD_singlesubject_h
#define GUARD_singlesubject_h

#include <Rcpp.h>
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
                         double univariate_pv_target_ratio);

#endif

