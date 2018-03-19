# 
# #' draw_sd_between_pulse
# #' 
# #' Metropolis Hastings algorithm for drawing the posterior distribution of the
# #' between pulse variability in terms of the standard deviation.
# #' 
# #' @param pulse_list A list of pulses 
# #' @param priors A list of all priors
# #' @param proposal_sd The standard deviation of the proposal distribution.
# #' @return A draw from the posterior distribution.  May end up being a list
# #' or an s3 object also containing an indicator for accepted/rejected.
# #' @import purrr
# #' @import stats
# #' @keywords pulse fit
# #' @examples
# #' this_pulse <- simulate_pulse()
# #' this_spec  <- pulse_spec()
# #' this_fit   <- fit_pulse(.data = this_pulse, iters = 1000, thin = 10,
# #'                         spec = this_spec)
# #' @export
# draw_sd_between_pulses <- function(pulse_list,
#                                    priors,
#                                    proposal_sd, 
#                                    old) {
# 
#   current_value <- pulse_list$parms$sd_mass
#   prior_max <- priors$max_sd$mass
# 
#   # draw new value
#   new_value <- rnorm(1, current_value, proposal_sd)
# 
#   if (new_value < 0 | new_value > prior_max) {
#     return(as_mh_value(current_value, 0))
#   }
# 
#   parms  <- pulse_list$parms
#   pulses <- map(pulse_list$pulses, ~ .x$parameters) %>% purrr::transpose(.) %>%
#     map(~ do.call(c, .x))
# 
#   first_part  <- pulse_list$num_pulses * (log(current_value) - log(new_value))
#   second_part <- 0.5 * ((1 / current_value^2) - (1 / new_value^2))
#   third_part  <- sum(pulses$mass_eta * (pulses$mass - parms$mean_mass)^2)
#   old_int <- sum(pnorm(pulse_list$mean_mass * sqrt(pulses$mass_eta) / current_value, log.p = TRUE))
#   new_int <- sum(pnorm(pulse_list$mean_mass * sqrt(pulses$mass_eta) / new_value, log.p = TRUE))
# 
#   rho <- old_int - new_int + first_part + second_part * third_part
#   rho <- min(0, rho)
# 
#   if (log(runif(1, min = 0, max = 1)) < rho) {
#     return(as_mh_value(new_value, 1))
#   } else {
#     return(as_mh_value(current_value, 0))
#   }
# 
# }
# 
# #' @export
# # Create a struct of this?
# as_mh_value <- function(value, accepted) {
#   structure(list("value" = value, "accepted" = accepted),
#             class = "mh_value")
# }
# 
# 
