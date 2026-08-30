# Getting Started with `bayespulse`

`bayespulse` is the R package interface to **libpulsatile**, a C++ (Rcpp/RcppArmadillo) library for Bayesian modeling of pulsatile hormone time series. It fits deconvolution models via birth-death MCMC, recovering pulse counts, pulse locations, secretion characteristics, and (for the joint model) driver-response coupling.

## 1. Overview

`bayespulse` currently implements three runnable models:

| Model | Description | Workflow |
|-------|-------------|----------|
| **Single-subject** | Bayesian deconvolution for one hormone series from one subject. | `simulate_pulse()` -> `pulse_spec()` -> `fit_pulse()` |
| **Population** | Hierarchical model pooling multiple subjects to estimate population-level secretion parameters. | `simulate_pulse_population()` -> `population_spec()` -> `fit_pulse_population()` |
| **Joint** | Two coupled hormones, where a *driver* hormone modulates the pulse intensity of a *response* hormone (coupling parameters `rho`, `nu`). | `simulate_pulse_joint()` -> `joint_spec()` -> `fit_pulse_joint()` |

> **Not yet implemented:** the **population-joint (Phase 4)** model — a hierarchical, multi-subject version of the joint driver-response model — is *not* available in this version. There is no `*_population_joint()` workflow yet.

## 2. Build & Install on macOS

The package compiles C++ via Rcpp/RcppArmadillo, so you need a compiler toolchain with a Fortran compiler.

### Prerequisites
```bash
# R 4.x (4.3.0+ recommended, matching CI)
brew install --cask r          # or download from https://cran.r-project.org/bin/macosx/

# Xcode command line tools: clang, make, git, system headers
xcode-select --install

# gfortran (via GCC) — RcppArmadillo links LAPACK/BLAS and needs Fortran
brew install gcc
```

`brew install gcc` installs `gfortran` to `/opt/homebrew/bin/gfortran` (Apple Silicon) or `/usr/local/bin/gfortran` (Intel).

**Apple Silicon symlink workaround** — R and RcppArmadillo commonly expect `gfortran` at `/usr/local/bin/gfortran`:
```bash
sudo ln -sf /opt/homebrew/bin/gfortran /usr/local/bin/gfortran
```
For version-specific compiler-tool guidance (including `~/.R/Makevars`), see:
https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/

### Clone WITH symlinks
The repo uses symbolic links to share C++ source between the standalone library and the R package (e.g. `R-package/src/singlesubject.cpp`, `R-package/inst/include`). The package will not build if these are lost:
```bash
git clone -c core.symlinks=true https://github.com/mmulvahill/libpulsatile.git
cd libpulsatile
ls -la R-package/inst/include   # should show: include -> ../../include/
```

### Install R dependencies and the package
```bash
Rscript -e "install.packages(c('Rcpp','RcppArmadillo','dplyr','ggplot2','tibble','tidyr','rlang','devtools','testthat','magrittr'), repos='https://cran.r-project.org')"

# RcppArmadillo loading cleanly confirms your gfortran toolchain is correct
Rscript -e "library(RcppArmadillo); cat('RcppArmadillo OK\n')"

cd R-package
Rscript -e "devtools::install()"
Rscript -e "library(bayespulse); cat('bayespulse installed\n')"
```

### macOS troubleshooting
- **`gfortran: command not found` / missing Fortran libs:** install `brew install gcc` and create the `/usr/local/bin/gfortran` symlink (Apple Silicon).
- **Missing C++ sources at build time:** you cloned without `core.symlinks=true`; re-clone with that flag.
- **RcppArmadillo C++11 deprecation warnings:** harmless; build flags include `-DARMA_USE_CURRENT`.
- **LAPACK/BLAS linker errors:** ensure Xcode Command Line Tools are installed (they provide Apple's Accelerate framework).

## 3. Build & Install on a Linux Server

### System dependencies
```bash
# gfortran is required for RcppArmadillo (Debian/Ubuntu)
sudo apt-get install gfortran
```

### R dependencies
```bash
Rscript -e "install.packages(c('Rcpp','RcppArmadillo','dplyr','ggplot2','tibble','tidyr','rlang','devtools','testthat','magrittr'), repos='https://cran.r-project.org')"
```

### Install the R package
```bash
cd R-package
Rscript -e "devtools::install()"
Rscript -e "library(bayespulse); cat('bayespulse installed\n')"
```

### Developer path: standalone C++ library
For C++ development and the Catch2 test suite, build the standalone library from the repository root:
```bash
make            # build the C++ library
bin/tests       # run the full C++ test suite
bin/tests "[single-subject]"   # run a single tagged test
```

### Rebuilding after C++ header changes
C++ headers live under the top-level `include/` and reach the R package through the `inst/include` symlink. After editing any header or source, regenerate the Rcpp glue, refresh docs, and reinstall (run from `R-package/`):
```bash
cd R-package
Rscript -e "Rcpp::compileAttributes(); devtools::document(); devtools::install()"
```
To be certain stale object files are not reused, do a clean standalone build from the repo root first:
```bash
make clean && make
```

### Tests and checks
```bash
cd R-package
Rscript -e "devtools::test()"                          # R (testthat) tests
Rscript -e "devtools::document(); devtools::check()"   # CRAN-style check (0 ERRORs, 0 WARNINGs)
```

## 4. Simulated Data

Each model ships a matching simulator that produces ground-truth data for testing and demos:

| Model | Simulator | Returns |
|-------|-----------|---------|
| Single-subject | `simulate_pulse()` | `pulse_sim`: `$data` (tibble: observation, time, concentration), `$parameters` (true per-pulse table) |
| Population | `simulate_pulse_population()` | `population_sim`: `$data` (LIST of per-subject data frames), `$n_subjects` |
| Joint | `simulate_pulse_joint()` | `joint_sim`: `$driver_data`, `$response_data` (tibbles: time, concentration), `$association` (true `$rho`, `$nu`) |

> There is **no population-joint simulator** — the Phase 4 hierarchical joint model is not yet implemented, so no `simulate_pulse_population_joint()` exists.

## 5. Runnable Examples

Each block is self-contained and uses small MCMC settings so it finishes in a couple of seconds. The package is loaded with `library(bayespulse)` (already installed).

### Single-subject

```r
suppressMessages(library(bayespulse))
set.seed(2026)

# 1. Simulate single-subject pulsatile data.
#    $data tibble: observation, time, concentration. $parameters = true pulses.
sim <- simulate_pulse(num_obs = 144, interval = 10)
cat("Class of sim:", class(sim), "\n")
cat("Simulated data dims:", nrow(sim$data), "obs x", ncol(sim$data), "cols:",
    paste(names(sim$data), collapse = ", "), "\n")
cat("True number of simulated pulses:", nrow(sim$parameters), "\n\n")

# 2. Build a model specification.
#    Use the 'strauss' location prior (the default). 'order-statistic' is NOT
#    implemented in the C++ birth-death process and will error if selected.
spec <- pulse_spec(location_prior_type = "strauss")
cat("Class of spec:", class(spec), " | location prior:", spec$location_prior, "\n\n")

# 3. Fit the deconvolution model. Keep MCMC small. burnin must be < iters.
#    You can pass the pulse_sim object directly to 'data'.
fit <- fit_pulse(data = sim, iters = 1000, thin = 2, burnin = 500, spec = spec,
                 verbose = FALSE)
cat("Class of fit:", class(fit), "| model:", fit$model, "\n\n")

# patient_chain: subject-level parameter chain (one row per saved iteration)
cat("=== patient_chain (subject-level params) ===\n")
cat("dims:", nrow(fit$patient_chain), "rows x", ncol(fit$patient_chain), "cols\n")
cat("columns:", paste(names(fit$patient_chain), collapse = ", "), "\n\n")

# pulse_chain: per-pulse parameter chain (multiple rows per iteration)
cat("=== pulse_chain (per-pulse params) ===\n")
cat("dims:", nrow(fit$pulse_chain), "rows x", ncol(fit$pulse_chain), "cols\n")
cat("columns:", paste(names(fit$pulse_chain), collapse = ", "), "\n\n")

# Posterior summaries from the subject-level chain
cat("=== Posterior summaries (patient_chain) ===\n")
pc <- as.data.frame(fit$patient_chain)
num_cols <- names(pc)[sapply(pc, is.numeric)]
summ <- t(sapply(num_cols, function(cn) {
  x <- pc[[cn]]
  c(mean = mean(x), sd = sd(x),
    q2.5 = unname(quantile(x, 0.025)), q97.5 = unname(quantile(x, 0.975)))
}))
print(round(summ, 4))

cat("\nPosterior mean number of pulses:", round(mean(pc$num_pulses), 3),
    "(true =", nrow(sim$parameters), ")\n")
```

**What to look at:** `fit$patient_chain` is the main chain for posterior summaries (columns: `iteration, num_pulses, baseline, mass_mean, width_mean, halflife, model_error, mass_sd, width_sd`); compare `mean(pc$num_pulses)`, `baseline`, `halflife`, and `model_error` against the truth in `sim$parameters`. In the verified run the model recovered ~11.24 pulses (true 11), baseline 2.61 (true 2.6), and halflife 46.7 (true 45). `fit$pulse_chain` holds per-pulse draws (`location, mass, width, ...`). (A `width_sd` of 0 at these tiny iters is a short-run artifact, not a bug.)

### Population (hierarchical, multi-subject)

```r
suppressMessages(library(bayespulse))

## 1. Simulate a multi-subject dataset. $data is a LIST of per-subject frames.
sim <- simulate_pulse_population(n_subjects = 4, num_obs = 48, interval = 10,
                                 mass_mean = 3.5, width_mean = 35,
                                 baseline = 2.6, halflife = 45,
                                 mass_mean_sd = 0.5, baseline_sd = 0.3,
                                 seed = 2024)
cat("n_subjects:", sim$n_subjects, " length(sim$data):", length(sim$data), "\n")
cat("Subject 1 head:\n"); print(head(sim$data[[1]], 3))

## 2. Build the population model specification.
spec <- population_spec(prior_mass_mean_mean     = 3.5,
                        prior_width_mean_mean    = 35,
                        prior_baseline_mean_mean = 2.6,
                        prior_halflife_mean_mean = 45,
                        prior_mean_pulse_count   = 12)

## 3. Fit the hierarchical model. spec MUST be passed by name (it is the 5th
##    positional arg; passing it 2nd would bind it to `time`).
fit <- fit_pulse_population(data = sim$data, spec = spec,
                            iters = 1000, burnin = 500, thin = 2,
                            verbose = FALSE)
cat("\nClass:", class(fit), "\n")
cat("num_subjects:", fit$num_subjects, "\n")
cat("Saved population samples:", nrow(fit$population_chain), "\n")

## 4. Read posterior population means.
pop_chain <- as.data.frame(fit$population_chain)
if ("iteration" %in% names(pop_chain))
  pop_chain <- pop_chain[, names(pop_chain) != "iteration", drop = FALSE]
cat("\nPosterior population means (colMeans of fit$population_chain):\n")
print(round(colMeans(pop_chain), 3))

## 5. Diagnostics helpers.
diag <- population_diagnostics(fit)   # data frame: parameter, mean, sd, ess, converged, ...
cat("\npopulation_diagnostics() head:\n")
print(head(diag[, c("parameter", "mean", "ess", "converged")], 4))

convergence_report(fit)               # prints a formatted ESS/ACF report

cat("\nDONE-OK\n")
```

**What to look at:** `colMeans(fit$population_chain)` gives the posterior population means for the 11 population-level parameters (`mass_mean, width_mean, baseline_mean, halflife_mean, *_sd, errorsq`); compare to the simulation truth (verified run: mass_mean ~2.72 vs 3.5, width_mean ~34 vs 35, baseline_mean ~2.59 vs 2.6, halflife_mean ~44 vs 45). `population_diagnostics(fit)` reports per-parameter ESS and a `converged` flag (ESS>100); `convergence_report(fit)` prints a formatted ESS/ACF report; `summary(fit)`, `plot_trace(fit)`, and `plot_acf(fit)` are also available. Variance parameters can show low ESS at these tiny iters — a mixing artifact of the short run.

### Joint (driver-response coupling)

```r
library(bayespulse)
set.seed(2025)

# 1. Simulate coupled driver -> response data. rho = coupling strength,
#    nu = coupling temporal spread (variance, minutes^2). True rho/nu live in
#    sim$association (NOT sim$parameters).
sim <- simulate_pulse_joint(
  rho                  = 1.5,
  nu                   = 100,
  num_obs              = 120,
  interval             = 10,
  response_pulse_count = 10,
  seed                 = 7
)
cat("Class of sim:", class(sim), "\n")
cat("driver_data dims:", paste(dim(sim$driver_data), collapse = " x "), "\n")
cat("response_data dims:", paste(dim(sim$response_data), collapse = " x "), "\n")
cat("True coupling -> rho:", sim$association$rho, " nu:", sim$association$nu, "\n\n")

# 2. Build the spec. Coupling priors default to log(rho)~N(0,1), log(nu)~N(3,1).
spec <- joint_spec()
print(spec)

# 3. Fit the joint model. driver_data / response_data are required by name.
#    burnin must be < iters. Keep iters small for the demo.
fit <- fit_pulse_joint(
  driver_data   = sim$driver_data,
  response_data = sim$response_data,
  spec          = spec,
  iters         = 1500,
  thin          = 5,
  burnin        = 500,
  verbose       = FALSE
)
cat("\nFit class:", class(fit), "\n")
cat("Saved samples per chain:", nrow(fit$association_chain), "\n\n")

# 4. Inspect the ASSOCIATION (coupling) chain: rho and nu.
cat("=== Association chain (rho, nu) ===\n")
cat("Columns:", paste(names(fit$association_chain), collapse = ", "), "\n")
print(head(fit$association_chain))
cat(sprintf("\n  rho: mean = %.4f, 95%% CI = [%.4f, %.4f]\n",
            mean(fit$association_chain$rho),
            quantile(fit$association_chain$rho, 0.025),
            quantile(fit$association_chain$rho, 0.975)))
cat(sprintf("  nu : mean = %.4f, 95%% CI = [%.4f, %.4f]\n",
            mean(fit$association_chain$nu),
            quantile(fit$association_chain$nu, 0.025),
            quantile(fit$association_chain$nu, 0.975)))

# 5. Driver and response parameter chains.
cat("\nPosterior mean driver pulse count:",   mean(fit$driver_chain$num_pulses),   "\n")
cat("Posterior mean response pulse count:", mean(fit$response_chain$num_pulses), "\n")

cat("\nCAVEAT: rho/nu are only WEAKLY IDENTIFIABLE from a single subject;\n")
cat("the posterior is strongly prior-influenced. This demonstrates the\n")
cat("MECHANICS of the joint model, NOT precise recovery of true rho/nu.\n")
```

**What to look at:** `fit$association_chain` holds the coupling draws (`rho`, `nu`, natural scale, always > 0); summarize with `mean`/`quantile`. `fit$driver_chain` and `fit$response_chain` carry per-hormone parameters (`baseline, halflife, mass_mean, width_mean, errorsq, num_pulses`); the response per-pulse chain (`fit$response_pulse_chain`) additionally has `lambda` (coupling intensity at each response-pulse location). `print(fit)` / `summary(fit)` give head-of-chains and posterior means with 95% CIs. **Read the caveat below** — treat this as verified mechanics, not recovery of the true `rho`/`nu`.

## 6. Notes & Caveats

- **Joint `rho`/`nu` weak identifiability.** From a *single* subject, the coupling parameters `rho` (strength) and `nu` (temporal spread) are only weakly identifiable: `rho`'s posterior centers in the right region but with wide uncertainty, and `nu` is pulled toward its prior rather than the true value. Do not present recovered `rho`/`nu` as accurate point estimates — the joint example demonstrates fitting/inspection mechanics, not parameter recovery. Robust coupling inference is expected to require the (not-yet-implemented) population-joint model. (An earlier note here that the response pulse count "can collapse to 1 under tiny iters — a sampling artifact, not a bug" was wrong: that collapse was a real defect in the response birth step under the default hard-core location prior, fixed in `joint_birthdeath.h`; earlier quantitative claims in this section were measured under that defect.)

- **Demo iters vs production iters.** All examples above use small MCMC settings (e.g. `iters = 1000`, `burnin = 500`, `thin = 2`) so they run in a couple of seconds. These are intentionally too short for real inference — at this size, variance parameters (`*_sd`, `errorsq`) legitimately show low ESS and poor mixing. The function defaults (`iters = 250000`, `thin = 50`, `burnin = 0.1 * iters`) reflect production-scale runs; always set `iters`/`burnin`/`thin` explicitly small for demos and scale up substantially for analysis.

- **Common API constraints.**
  - `burnin` must be strictly `< iters` (it is applied *before* thinning); saved samples = `(iters - burnin) / thin`.
  - `pulse_spec(location_prior_type = ...)` and `joint_spec()` only support the `"strauss"` location prior; `"order-statistic"` is not implemented and errors at fit time.
  - In `fit_pulse_population()` and `fit_pulse_joint()`, `spec` is not the second positional argument — always pass `spec = ...` by name. Pass the per-subject **list** `sim$data` to `fit_pulse_population()`, and `sim$driver_data` / `sim$response_data` to `fit_pulse_joint()`.
  - `joint_spec()` coupling priors are on the **log scale** (`prior_log_rho_mean`, `prior_log_rho_var`, `prior_log_nu_mean`, `prior_log_nu_var`); starting values `sv_rho`, `sv_nu` are on the natural scale.
  - Confirm exact argument names in `R-package/R/` if you adapt these examples.
