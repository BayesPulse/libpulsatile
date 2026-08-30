#!/usr/bin/env Rscript
#===============================================================================
# population_performance.R - Performance benchmarks for population model
#===============================================================================
#
# Tests population model performance with varying:
# - Number of subjects (2, 5, 10, 25, 50)
# - MCMC iterations (1000, 5000, 10000)
# - Observations per subject (24, 48, 96)
#
# Outputs:
# - Runtime statistics
# - Memory usage
# - Convergence metrics
# - Recommendations for production use
#
# Usage: Rscript benchmarks/population_performance.R
#

library(bayespulse)

cat("\n")
cat("========================================================================\n")
cat("Population Model Performance Benchmarks\n")
cat("========================================================================\n")
cat("\n")

# Get system info
cat("System Information:\n")
cat("  R version:", R.version.string, "\n")
cat("  Platform:", R.version$platform, "\n")
cat("  CPU cores:", parallel::detectCores(), "\n")
cat("\n")

set.seed(2024)

# Benchmark configurations
configs <- expand.grid(
  n_subjects = c(2, 5, 10, 25),
  n_obs = c(24, 48),
  iters = c(1000, 5000),
  stringsAsFactors = FALSE
)

# Add computed columns
configs$total_obs <- configs$n_subjects * configs$n_obs
configs$config_id <- paste0("S", configs$n_subjects, "_O", configs$n_obs, "_I", configs$iters)

# Results storage
results <- data.frame(
  config_id = character(),
  n_subjects = integer(),
  n_obs = integer(),
  total_obs = integer(),
  iters = integer(),
  runtime_sec = numeric(),
  runtime_per_iter = numeric(),
  memory_mb = numeric(),
  min_ess = numeric(),
  median_ess = numeric(),
  pct_converged = numeric(),
  stringsAsFactors = FALSE
)

cat("Benchmark Configurations:\n")
cat("  Total configs:", nrow(configs), "\n")
cat("  Subject counts:", paste(unique(configs$n_subjects), collapse = ", "), "\n")
cat("  Obs per subject:", paste(unique(configs$n_obs), collapse = ", "), "\n")
cat("  MCMC iterations:", paste(unique(configs$iters), collapse = ", "), "\n")
cat("\n")

cat("========================================================================\n")
cat("Running Benchmarks\n")
cat("========================================================================\n\n")

for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, ]

  cat(sprintf("[%2d/%2d] Config: %s\n", i, nrow(configs), cfg$config_id))
  cat(sprintf("        %d subjects, %d obs/subject, %d iterations\n",
              cfg$n_subjects, cfg$n_obs, cfg$iters))

  # Create simulated data
  cat("        Simulating data...")
  sim_data <- lapply(seq_len(cfg$n_subjects), function(j) {
    simulate_pulse(
      num_obs = cfg$n_obs,
      interval = 10,
      mass_mean = 3.5,
      width_mean = 35,
      constant_baseline = 2.6,
      constant_halflife = 45
    )
  })
  cat(" done\n")

  # Create specification
  spec <- population_spec()

  # Benchmark
  cat("        Running MCMC...")
  gc()  # Garbage collect before timing
  mem_before <- as.numeric(pryr::mem_used()) / 1024^2  # MB

  start_time <- Sys.time()
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = cfg$iters,
    thin = 10,
    burnin = as.integer(cfg$iters * 0.5),
    verbose = FALSE
  )
  end_time <- Sys.time()

  mem_after <- as.numeric(pryr::mem_used()) / 1024^2  # MB
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  cat(" done\n")
  cat(sprintf("        Runtime: %.2f sec (%.4f sec/iter)\n",
              runtime, runtime / cfg$iters))

  # Calculate convergence metrics
  cat("        Computing diagnostics...")
  diag <- population_diagnostics(fit)

  min_ess <- min(diag$ess, na.rm = TRUE)
  median_ess <- median(diag$ess, na.rm = TRUE)
  pct_converged <- 100 * mean(diag$converged, na.rm = TRUE)

  cat(" done\n")
  cat(sprintf("        ESS: min=%.0f, median=%.0f, converged=%.1f%%\n",
              min_ess, median_ess, pct_converged))

  # Store results
  results <- rbind(results, data.frame(
    config_id = cfg$config_id,
    n_subjects = cfg$n_subjects,
    n_obs = cfg$n_obs,
    total_obs = cfg$total_obs,
    iters = cfg$iters,
    runtime_sec = runtime,
    runtime_per_iter = runtime / cfg$iters,
    memory_mb = mem_after - mem_before,
    min_ess = min_ess,
    median_ess = median_ess,
    pct_converged = pct_converged,
    stringsAsFactors = FALSE
  ))

  cat("\n")
}

cat("========================================================================\n")
cat("Benchmark Results Summary\n")
cat("========================================================================\n\n")

# Print results table
print(results, row.names = FALSE)

cat("\n")
cat("========================================================================\n")
cat("Performance Analysis\n")
cat("========================================================================\n\n")

# Runtime scaling with subjects
cat("Runtime scaling with number of subjects:\n")
for (n_obs in unique(results$n_obs)) {
  for (iters in unique(results$iters)) {
    subset_results <- results[results$n_obs == n_obs & results$iters == iters, ]
    if (nrow(subset_results) > 1) {
      fit_lm <- lm(runtime_sec ~ n_subjects, data = subset_results)
      slope <- coef(fit_lm)[2]
      intercept <- coef(fit_lm)[1]
      cat(sprintf("  O=%d, I=%d: Runtime = %.2f + %.2f * n_subjects (R²=%.3f)\n",
                  n_obs, iters, intercept, slope, summary(fit_lm)$r.squared))
    }
  }
}

cat("\n")
cat("Convergence summary:\n")
cat(sprintf("  Best median ESS: %.0f (config: %s)\n",
            max(results$median_ess), results$config_id[which.max(results$median_ess)]))
cat(sprintf("  Worst median ESS: %.0f (config: %s)\n",
            min(results$median_ess), results$config_id[which.min(results$median_ess)]))
cat(sprintf("  Configs with >80%% converged: %d / %d\n",
            sum(results$pct_converged > 80), nrow(results)))

cat("\n")
cat("========================================================================\n")
cat("Recommendations\n")
cat("========================================================================\n\n")

# Find fastest config with good convergence
good_conv <- results[results$median_ess > 100, ]
if (nrow(good_conv) > 0) {
  fastest_good <- good_conv[which.min(good_conv$runtime_sec), ]
  cat("For quick testing:\n")
  cat(sprintf("  %d subjects, %d obs/subject, %d iterations\n",
              fastest_good$n_subjects, fastest_good$n_obs, fastest_good$iters))
  cat(sprintf("  Expected runtime: %.1f seconds\n", fastest_good$runtime_sec))
  cat(sprintf("  Median ESS: %.0f\n\n", fastest_good$median_ess))
}

# Find best convergence
best_conv <- results[which.max(results$median_ess), ]
cat("For best convergence:\n")
cat(sprintf("  %d subjects, %d obs/subject, %d iterations\n",
            best_conv$n_subjects, best_conv$n_obs, best_conv$iters))
cat(sprintf("  Expected runtime: %.1f seconds\n", best_conv$runtime_sec))
cat(sprintf("  Median ESS: %.0f\n\n", best_conv$median_ess))

# Production recommendations
cat("For production use (targeting ESS > 400):\n")
cat("  Recommended: 10,000-25,000 iterations\n")
cat("  Thin: 10-20\n")
cat("  Burnin: 50% of iterations\n")
cat("  Expected runtime (10 subjects): ~5-15 minutes\n")
cat("  Expected runtime (50 subjects): ~25-75 minutes\n")

cat("\n")
cat("========================================================================\n")
cat("Benchmark Complete\n")
cat("========================================================================\n\n")

# Save results
output_file <- "benchmarks/population_performance_results.csv"
write.csv(results, output_file, row.names = FALSE)
cat("Results saved to:", output_file, "\n\n")
