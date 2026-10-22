#-------------------------------------------------------------------------------
# sampler_equivalence.R
#
# Compares two builds of bayespulse -- typically a reference (e.g. master) and a
# candidate -- on identical seeded fits of the single-subject, population and
# joint models. For each scenario it reports the best-of-`reps` wall time of
# each build and the largest relative difference between their chains.
#
# A performance change that leaves the sampler's arithmetic unchanged should
# give a max relative difference at floating-point rounding level (~1e-16):
# every accept/reject and birth/death decision is then identical, so the
# chains agree draw for draw.
#
# Install each build into its own library first, e.g.
#   R CMD INSTALL -l /tmp/lib_ref  <reference checkout>/R-package
#   R CMD INSTALL -l /tmp/lib_new  <candidate checkout>/R-package
# then run
#   Rscript sampler_equivalence.R /tmp/lib_ref /tmp/lib_new [iters] [reps]
#
# Each fit runs in a fresh R process (the two builds share a package name).
#-------------------------------------------------------------------------------

args  <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript sampler_equivalence.R <lib_ref> <lib_new> ",
       "[iters] [reps]")
}
libs  <- c(reference = args[1], candidate = args[2])
iters <- if (length(args) >= 3) as.integer(args[3]) else 2000L
reps  <- if (length(args) >= 4) as.integer(args[4]) else 3L

scenarios <- c("single_n144", "single_n288", "single_n576",
               "population_3x144", "joint_n144")

# Code run in the child process: fit one scenario, save its numeric chains
child <- function(lib, scenario, iters, out) {
  suppressMessages(library(bayespulse, lib.loc = lib))
  options(bayespulse.warn_burnin_adaptation = FALSE)
  numeric_chains <- function(fit) {
    res <- list()
    walk <- function(x, path) {
      if (is.data.frame(x) || is.matrix(x)) {
        x <- as.data.frame(x)
        res[[path]] <<- as.matrix(x[vapply(x, is.numeric, logical(1))])
      } else if (is.list(x)) {
        nms <- names(x)
        if (is.null(nms)) nms <- rep("", length(x))
        for (k in seq_along(x)) walk(x[[k]], paste0(path, "/", k, nms[k]))
      }
    }
    walk(unclass(fit), "")
    res
  }
  quiet <- function(expr) invisible(utils::capture.output(expr))
  if (startsWith(scenario, "single")) {
    n <- as.integer(sub("single_n", "", scenario))
    set.seed(1)
    sim  <- simulate_pulse(num_obs = n, interval = 10, ipi_mean = 12)
    spec <- pulse_spec(prior_mean_pulse_count = round(n / 12))
    set.seed(2024)
    time <- system.time(quiet(
      fit <- fit_pulse(sim$data, spec = spec, iters = iters, thin = 1,
                       burnin = 0)))[["elapsed"]]
  } else if (scenario == "population_3x144") {
    set.seed(123)
    sims <- lapply(1:3, function(i) {
      simulate_pulse(num_obs = 144, interval = 10, ipi_mean = 12,
                     mass_mean = 3.5, width_mean = 35)
    })
    spec <- population_spec(pulse_distribution = "truncnorm",
                            prior_mass_mean_mean = 3.5,
                            prior_width_mean_mean = 35,
                            prior_baseline_mean_mean = 2.6,
                            prior_halflife_mean_mean = 45,
                            prior_mean_pulse_count = 12)
    set.seed(7)
    time <- system.time(quiet(
      fit <- fit_pulse_population(data = sims, spec = spec, iters = iters,
                                  thin = 1, burnin = 0)))[["elapsed"]]
  } else {
    sim <- simulate_pulse_joint(rho = 1.5, nu = 100, response_pulse_count = 10,
                                num_obs = 144, interval = 10,
                                pulse_distribution = "lognormal", seed = 21)
    set.seed(11)
    time <- system.time(quiet(
      fit <- fit_pulse_joint(sim$driver_data, sim$response_data,
                             spec = joint_spec(prior_response_pulse_count = 10),
                             iters = iters, thin = 1, burnin = 0)))[["elapsed"]]
  }
  saveRDS(list(time = time, chains = numeric_chains(fit)), out)
}

run_child <- function(lib, scenario, iters) {
  out    <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  writeLines(c(paste0("child <- ", paste(deparse(child), collapse = "\n")),
               sprintf("child(%s, %s, %dL, %s)", deparse(lib),
                       deparse(scenario), iters, deparse(out))), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), script,
                    stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("fit failed for ", scenario, " with ", lib)
  readRDS(out)
}

max_rel_diff <- function(a, b) {
  stopifnot(identical(names(a), names(b)))
  worst <- 0
  for (nm in names(a)) {
    if (!identical(dim(a[[nm]]), dim(b[[nm]]))) return(Inf)
    if (length(a[[nm]])) {
      d <- abs(a[[nm]] - b[[nm]]) / (1 + abs(a[[nm]]))
      worst <- max(worst, d, na.rm = TRUE)
    }
  }
  worst
}

results <- data.frame()
for (scenario in scenarios) {
  # Long series are slow under the reference; fewer iterations keep it bounded
  n_iter <- if (scenario == "single_n576") max(100L, iters %/% 4L) else iters
  best <- c(reference = Inf, candidate = Inf)
  fits <- list()
  for (r in seq_len(reps)) {
    for (b in names(libs)) {
      fits[[b]] <- run_child(libs[[b]], scenario, n_iter)
      best[[b]] <- min(best[[b]], fits[[b]]$time)
    }
  }
  results <- rbind(results, data.frame(
    scenario       = scenario,
    iters          = n_iter,
    reference_sec  = best[["reference"]],
    candidate_sec  = best[["candidate"]],
    speedup        = best[["reference"]] / best[["candidate"]],
    max_rel_diff   = max_rel_diff(fits$reference$chains, fits$candidate$chains)
  ))
  print(utils::tail(results, 1), row.names = FALSE)
}

cat("\n")
print(transform(results,
                reference_sec = round(reference_sec, 2),
                candidate_sec = round(candidate_sec, 2),
                speedup       = round(speedup, 1),
                max_rel_diff  = signif(max_rel_diff, 3)),
      row.names = FALSE)
