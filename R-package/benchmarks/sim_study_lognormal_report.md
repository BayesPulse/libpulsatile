# Simulation study: log-scale pulse-width SD recovery (papers-default fit)

Reproducible script: `benchmarks/sim_study_lognormal_width_sd.R`.

## Goal

Show that the papers-default single-subject fit -- `{lognormal random effects,
Uniform(0, .) SD prior, Gaussian pulses (kappa fixed at 1)}`, i.e. the
out-of-the-box `pulse_spec()` -- recovers a **known log-scale pulse-width SD**
(`sigma_width`) and pulse-mass SD (`sigma_mass`) from log-normal-consistent
simulated data, with a finite, non-divergent trace (contrast the old
natural-scale parameterization, whose pulse widths space-filled / ran away).

## Data-generating process

`simulate.R` now exposes a first-class `pulse_distribution` argument. Under
`pulse_distribution = "lognormal"` the per-pulse random effects are drawn

```
mass  = exp(rnorm(1, mass_mean,  mass_sd))   # log-scale N, no kappa, no truncation
width = exp(rnorm(1, width_mean, width_sd))
```

so `mass_mean / mass_sd / width_mean / width_sd` are LOG-scale parameters --
exactly the fit parameterization (`log theta ~ N(mu, sigma^2)` with kappa = 1).
The legacy natural-scale truncated-Student-t mixture remains the default
(`pulse_distribution = "truncnorm"`), so existing simulate/fit tests are
unaffected.

**Verification of consistency.** With truths `mu_width = 3.0`, `sigma_width = 0.5`,
`mu_mass = 1.2`, `sigma_mass = 0.5` (seed 20260706, `num_obs = 240`, ~20 pulses),
the realized pulse table has empirical log-width SD 0.435 and log-mass SD 0.452
-- i.e. `log(width)` and `log(mass)` from the simulator are Gaussian with the
requested mean/SD. The simulator's log-normal path is consistent with the fit.

## Study design

- Single long single-subject series: `num_obs = 240`, `interval = 10` (40 h),
  `ipi_mean = 12` -> ~20 pulses (enough pulse-to-pulse width observations to
  identify `sigma_width`), constant half-life 45 and baseline 2.6,
  `error_var = 0.005`.
- Truths (log scale): `mu_mass = 1.2`, `sigma_mass = 0.5`, `mu_width = 3.0`,
  `sigma_width = 0.5`.
- Fit: `pulse_spec(location_prior_type = "strauss", prior_mean_pulse_count = np)`
  (all other settings the papers defaults), 60000 iters, thin 20, burnin 20000
  (2000 saved draws). `set.seed(20260706)` -- reproducible.

(The order-statistic location prior is not implemented in the birth-death step,
so the Strauss prior is used, matching the default.)

## Results (with the retuned default proposal variance)

Posterior mean and 95% credible interval vs. the log-scale truth
(ESS = effective sample size from the initial-positive-sequence estimator).
These are the numbers the current `pulse_spec()` default produces (after the
proposal-variance retuning described below):

| parameter   | truth | post. mean | 95% CI            | covers | ESS  |
|-------------|-------|------------|-------------------|--------|------|
| width_sd    | 0.50  | 0.443      | [0.055, 1.214]    | yes    | 129  |
| mass_sd     | 0.50  | 0.423      | [0.302, 0.603]    | yes    | 1840 |
| width_mean  | 3.00  | 2.997      | [2.388, 3.494]    | yes    | 279  |
| mass_mean   | 1.20  | 1.232      | [1.022, 1.430]    | yes    | 627  |
| num_pulses  | 20    | 20.09      | (range 20-21)     | exact  | --   |

**All four log-scale parameters' 95% CIs cover the truth**, the pulse count is
recovered exactly, and `width_sd` is finite throughout (range 0.018-2.768, no
runaway) -- the fit recovers the generating parameters and the width trace is
bounded and non-divergent (the qualitative win over the old space-filling
natural-scale widths). `width_sd`'s posterior mean (0.443) sits close to both
the log-scale truth (0.50) and this dataset's realized empirical log-width SD
(0.435).

## Retuning the individual-pulse proposal variance

An initial run with the pre-retune default exposed a mixing problem, which was
diagnosed and fixed.

### Symptom

With `pv_indiv_pulse_width = 0.25` (the value carried over from the log-scale
mean/SD proposals), width mixing was poor while mass mixing was excellent:

| parameter  | ESS @ pv=0.25 | ESS @ pv=9 | post.mean @0.25 -> @9 |
|------------|---------------|------------|-----------------------|
| width_sd   | 36            | **129**    | 0.341 -> 0.443        |
| width_mean | 49            | **279**    | 3.061 -> 2.997        |
| mass_sd    | 2000          | 1840       | 0.422 -> 0.423        |
| mass_mean  | 856           | 627        | 1.237 -> 1.232        |

Coverage held in both cases (every CI covered truth), but at `pv=0.25` the
`width_sd` posterior was diffuse with heavy mass near 0 ([0.007, 1.005]) and its
autocorrelation time was ~55 thinned draws.

### Diagnosis

The individual-pulse random-effect update is a **symmetric random walk on the
NATURAL-scale pulse value** (`ss_draw_randomeffects.h`), even under the
log-normal parameterization. The natural pulse width is `exp(width_mean) ~ 20`,
while the natural pulse mass is `exp(mass_mean) ~ 3.3`. A proposal variance of
0.25 (step SD 0.5) is well matched to mass (~15% of its magnitude) but far too
small for width (~2.5% of its magnitude), so individual widths barely moved.
Because `width_sd`'s conditional draw depends on the scatter of the (sticky)
latent log-widths, it inherited the slow mixing. Proposal-variance adaptation
only rescales multiplicatively (+/-10% per 500 iters, frozen at iter 25000), so
from a start of 0.25 it needs ~37 adjustment steps (~18500 iters) just to reach
a width-appropriate scale -- nearly the whole adaptation budget.

### Fix

The log-scale `ln_defaults` in `pulse_spec()` (and, by the same argument,
`population_spec()` and `joint_spec()`) now scale the individual-pulse proposal
variances to the natural pulse magnitude `exp(mean)` instead of to log-scale
sizes:

- `pv_indiv_pulse_width`: `0.25 -> 9`  (width ~ exp(3) ~ 20)
- `pv_indiv_pulse_mass`:  `0.25 -> 0.5` (mass ~ exp(1.2) ~ 3.3)

The fixed-effect mean/SD proposals (`pv_mean_*`, `pv_sd_*`) are unchanged -- they
act on log-scale parameters and were already correctly scaled. The `truncnorm`
defaults are untouched, so the legacy path is unaffected. The retune lifts
`width_sd` ESS 3.6x and `width_mean` ESS 5.7x with no cost to mass mixing.

## Caveat: pulse width is the weakest-identified quantity

Even after retuning, pulse **width** mixes more slowly than pulse mass
(`width_sd` ESS 129 vs `mass_sd` ESS 1840 at equal cost), and its posterior is
wide and right-skewed. This is a genuine identifiability property of
deconvolution -- a pulse's mass (area) is far better constrained by the
concentration curve than its width (shape), so the pulse-to-pulse width SD is
estimated with real, irreducible uncertainty. The study therefore demonstrates
*unbiased recovery with a finite, mixing trace*, not a tight width_sd posterior;
tightening it further would require more pulses (a longer series) or more
iterations, not a different proposal tuning.

## Reproducing

```
cd R-package
Rscript benchmarks/sim_study_lognormal_width_sd.R   # ~15-20 min single fit
```

The script simulates the data (`pulse_distribution = "lognormal"`, seed 20260706),
fits with the default `pulse_spec()`, and prints the recovery table above.
