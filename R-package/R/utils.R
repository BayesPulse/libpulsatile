#
#adjust_proposal_variance <- function(current_psd,
#                                     target_acceptance_ratio,
#                                     current_acceptance_ratio) {
#
#  y <- 1.0 + 100.0 * (current_acceptance_ratio - target_acceptance_ratio) ^ 3
#
#  psd <- current_pv^2
#
#  if (y < 0.9) { 
#    rtn <- psd * 0.9
#  } else if (y > 1.1) {
#    rtn <- psd * 1.1
#  } else {
#    rtn <- psd
#  }
#
#  return(rtn)
#
#}
#
#current_psd <- 1.04881
#target_acceptance_ratio <- 0.3499999999998 
#current_acceptance_ratio <- 0 # not a single acceptance in the "Check tracking iterations and adjusting pv/psd" section
#
#adjust_proposal_variance(current_psd, target_acceptance_ratio, current_acceptance_ratio)
#


#' Warn when burn-in ends before proposal-variance adaptation stops
#'
#' The MH samplers rescale their proposal variances every
#' \code{pv_adjust_iter} iterations until \code{pv_adjust_max_iter}. Draws
#' saved before adaptation stops come from a chain whose transition kernel is
#' still changing, so warn when \code{burnin} does not cover that period. Set
#' \code{options(bayespulse.warn_burnin_adaptation = FALSE)} to silence.
#'
#' @param burnin Number of burn-in iterations.
#' @param pv_adjust_max_iter Last iteration at which proposal variances adapt.
#' @return \code{NULL}, invisibly.
#' @noRd
check_burnin_adaptation <- function(burnin, pv_adjust_max_iter) {
  if (burnin < pv_adjust_max_iter &&
      isTRUE(getOption("bayespulse.warn_burnin_adaptation", TRUE))) {
    warning(
      "burnin (", burnin, ") is less than pv_adjust_max_iter (",
      pv_adjust_max_iter, "): proposal variances are still adapting during ",
      "the saved iterations, so early saved samples are not from a fixed ",
      "MCMC kernel. Consider burnin >= ", pv_adjust_max_iter, ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}
