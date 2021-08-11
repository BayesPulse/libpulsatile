#' @importFrom dplyr summarise_all
#' 
#' 
#' @export
summary.pulse_fit <- function(object, ...) {
  args <- list(...)
  if(is.null(args$quantiles)) { args$quantiles <- c(.1, .2, .5, .8, .9) }
  pat_params <- object$patient_chain %>%
    select(-c(.data$iteration)) %>%
    summarise_all(quantile, args$quantiles)
  
  cat("Patient parameters\n")
  cat("               ")
  cat(paste0(args$quantiles*100, "%    "))
  cat("\n")
  cat("mass_mean:    ", paste(format(round(pat_params$mass_mean, 3), nsmall = 3), " "), "\n")
  cat("width_mean:   ", paste(format(round(pat_params$width_mean, 2), nsmall = 2), " "), "\n")
  cat("baseline:     ", paste(format(round(pat_params$baseline, 3), nsmall = 3), " "), "\n")
  cat("halflife:     ", paste(format(round(pat_params$halflife, 2), nsmall = 2), " "), "\n")
  cat("model_error:  ", paste(format(round(pat_params$model_error, 3), nsmall = 3), " "), "\n")
  cat("num_pulses:   ", paste(format(round(pat_params$num_pulses, 2), nsmall = 2), " "), "\n")
  
  invisible(pat_params)
}

#' @export
summary.pop_pulse_fit <- function(object, ...) {
  args <- list(...)
  if(is.null(args$quantiles)) { args$quantiles <- c(.1, .2, .5, .8, .9) }
  if(is.null(args$patient)) { args$patient <- 1 }
  # if(!is.null(patient)) {
  pop_params <- object$population_chain %>%
    select(-c(.data$iteration)) %>%
    summarise_all(quantile, args$quantiles)
    
  pat_params <- object$patient_chain[[args$patient]] %>%
    select(-c(.data$iteration, .data$likelihood)) %>%
    summarise_all(quantile, args$quantiles)
  
  cat("Population parameters\n")
  cat("               ")
  cat(paste0(args$quantiles*100, "%    "))
  cat("\n")
  cat("mass_mean:    ", paste(format(round(pop_params$mass_mean, 3), nsmall = 3), " "), "\n")
  cat("mass_s2s_sd:  ", paste(format(round(pop_params$mass_s2s_sd, 3), nsmall = 3), " "), "\n")
  cat("mass_p2p_sd:  ", paste(format(round(pop_params$mass_p2p_sd, 3), nsmall = 3), " "), "\n")
  cat("width_mean:   ", paste(format(round(pop_params$width_mean, 2), nsmall = 2), " "), "\n")
  cat("width_s2s_sd: ", paste(format(round(pop_params$width_s2s_sd, 3), nsmall = 3), " "), "\n")
  cat("width_p2p_sd: ", paste(format(round(pop_params$width_p2p_sd, 3), nsmall = 3), " "), "\n")
  cat("baseline_mean:", paste(format(round(pop_params$baseline_mean, 3), nsmall = 3), " "), "\n")
  cat("baseline_sd:  ", paste(format(round(pop_params$baseline_sd, 3), nsmall = 3), " "), "\n")
  cat("halflife_mean:", paste(format(round(pop_params$halflife_mean, 2), nsmall = 2), " "), "\n")
  cat("halflife_sd:  ", paste(format(round(pop_params$halflife_sd, 3), nsmall = 3), " "), "\n\n")
  
  cat("Patient", names(object$patient_chain)[args$patient], "parameters\n")
  cat("               ")
  cat(paste0(args$quantiles*100, "%    "))
  cat("\n")
  cat("mass_mean:    ", paste(format(round(pat_params$mass_mean, 3), nsmall = 3), " "), "\n")
  cat("width_mean:   ", paste(format(round(pat_params$width_mean, 2), nsmall = 2), " "), "\n")
  cat("baseline:     ", paste(format(round(pat_params$baseline, 3), nsmall = 3), " "), "\n")
  cat("halflife:     ", paste(format(round(pat_params$halflife, 2), nsmall = 2), " "), "\n")
  cat("model_error:  ", paste(format(round(pat_params$model_error, 3), nsmall = 3), " "), "\n")
  cat("num_pulses:   ", paste(format(round(pat_params$num_pulses, 2), nsmall = 2), " "), "\n")
  
  invisible(list("pop_quantiles" = pop_params,
              "pat_quantiles" = pat_params))
}

