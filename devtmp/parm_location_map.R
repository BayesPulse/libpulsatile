


patient_mean_mass <-
  c("patient->priors->mass_mean"     = "priors->fe_mean[0]",
    "patient->priors->mass_variance" = "priors->fe_variance[0]",
    "patient->estimates->mass"       = "parms->theta[0]",
    "proposal"                       = "theta[0]")

pulse_mass <-
  c("patient->pulses->mass"           = "node->theta[0]",
    "patient->pulses->tvarscale_mass" = "node->eta[0]",
    "patient->estimates->mass_sd"     = "parms->re_sd[0]")






