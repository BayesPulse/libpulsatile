# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

singlesubject_ <- function(concentration, time, location_prior, inpriors, proposalvars, startingvals, mcmc_iterations, thin, burnin, verbose, pv_adjust_iter, pv_adjust_max_iter, bivariate_pv_target_ratio, univariate_pv_target_ratio) {
    .Call(`_bayespulse_singlesubject_`, concentration, time, location_prior, inpriors, proposalvars, startingvals, mcmc_iterations, thin, burnin, verbose, pv_adjust_iter, pv_adjust_max_iter, bivariate_pv_target_ratio, univariate_pv_target_ratio)
}

