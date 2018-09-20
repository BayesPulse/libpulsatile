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
