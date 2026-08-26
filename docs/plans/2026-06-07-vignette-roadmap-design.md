# bayespulse Vignette Roadmap — Design

**Date:** 2026-06-07
**Status:** Approved design; implementation plan to follow (writing-plans).
**Decisions locked in:** full portfolio roadmap · fix-and-fold the existing WIP vignette · simulation-only data for now · progressive-ladder structure · document future ideas (I–O) but do **not** plan or implement them yet.

---

## 1. Context & goals

Build a vignette portfolio that fully demonstrates `bayespulse` to R users, serving **both** biostatisticians (methods, MCMC, diagnostics, study design) and endocrinologists (hormone biology, interpretation, coupling).

### Verified capability map (exercised end-to-end this session)

All three model families run today:

- **Single-subject deconvolution** — `simulate_pulse` → `pulse_spec` → `fit_pulse` → `bp_trace` / `bp_posteriors` / `bp_location_posterior` / `bp_predicted` / `predict.pulse_fit`.
- **Population / hierarchical** — `simulate_pulse_population` → `population_spec` → `fit_pulse_population` → `population_diagnostics` / `summary.population_fit` / `convergence_report`.
- **Joint driver→response** — `simulate_pulse_joint` → `joint_spec` → `fit_pulse_joint` → `summary.joint_fit`. Couples a driver hormone's pulses to the response hormone's pulse *occurrence* via coupling strength `rho` and timescale `nu`.
- **MCMC diagnostics** — `effective_sample_size`, `gelman_rubin`, `plot_trace`, `plot_acf`, `convergence_report`.
- **Chain extraction** — `patient_chain`, `pulse_chain`, `chains`.

### Audiences and decisions

- Serve biostatisticians and endocrinologists in one coherent ladder; flagship joint-model vignette deliberately framed for endocrinologists.
- **Simulation-only** data for now (`simulate_pulse*`); no real dataset bundled (deferred).
- The orphaned `vignettes/my-vignette.Rmd` is **fixed and folded** into vignette A (salvage its intro/math, drop WIP notes, fix runtime).

---

## 2. Constraints

- **BD-MCMC is slow.** The existing vignette's default `fit_pulse()` (250,000 iters) ran in **~27 minutes** for a single subject. Live MCMC during `R CMD build` is infeasible.
- **Mixing is imperfect** (half-life especially). Low-iteration runs are visibly under-converged (median ESS ≈ 6–64 at 1k–5k iters), so we must *not* showcase toy under-mixed runs.
- **Vignettes aren't wired into the package today.** No `R-package/vignettes/` directory, empty `Suggests:`, no `VignetteBuilder` field. `my-vignette.Rmd` is never built or checked.
- **CRAN build-time limits** cap total vignette build time; each must build fast.
- **Weakly-identified pulse parameterization** (per `docs/lit_review.md`) means reasonably informative priors matter — an honest topic for the priors/diagnostics vignettes, not something to hide.

---

## 3. Infrastructure & runtime strategy

### Wiring (step zero of implementation)

- Create `R-package/vignettes/` (follow the repo's existing symlink convention where applicable).
- Add `VignetteBuilder: knitr` to `DESCRIPTION`; add `knitr`, `rmarkdown`, `ggplot2` (and plot deps) to `Suggests:`.
- Retire/replace root `vignettes/my-vignette.Rmd` after salvaging its content into vignette A.

### Precompute (`.Rmd.orig`) pattern — applies to every vignette with a non-trivial fit

1. Keep a `NN-name.Rmd.orig` that runs the **real, well-converged** MCMC (production `iters`) offline.
2. Knit it once locally to `NN-name.Rmd` with results/figures embedded; save fitted objects to `inst/extdata/*.rds`.
3. Ship the rendered `.Rmd` plus a `precompute.R` driver script and an in-file note explaining the pattern.
4. `R CMD build` then runs **no live MCMC** (builds in seconds) while readers still see genuine converged posteriors.

This is the established CRAN-friendly approach (used by sf, ggplot2, etc.) and resolves both the 27-minute problem and the "low-iter looks bad" problem. Each vignette spec below separates **precompute cost** (offline only) from **build cost** (≈ 0).

---

## 4. The eight current vignettes (A–H) — progressive ladder

Ordering is encoded via `\VignetteIndexEntry` numbering (`01-…` through `08-…`).

| # | Slug | Vignette | Key functions | Primary hook |
|---|------|----------|---------------|--------------|
| A | `01-getting-started` | Single-subject deconvolution *(folds in `my-vignette.Rmd`)* | `simulate_pulse`, `pulse_spec`, `fit_pulse`, `bp_*` | The on-ramp; cached fit builds in seconds. |
| B | `02-anatomy-of-the-model` | The generative process | `simulate_pulse`, `plot.pulse_sim` | See baseline + secretion + elimination + error; parameter sweeps. |
| C | `03-priors-and-specification` | Priors & starting values | `pulse_spec`, `population_spec` | Physiological meaning of priors; informative vs. vague; Strauss prior. |
| D | `04-convergence-diagnostics` | Assessing convergence | `gelman_rubin`, `effective_sample_size`, `plot_trace`, `plot_acf`, `convergence_report` | Honest mixing treatment; multi-chain R-hat; remedies. |
| E | `05-population-model` | Hierarchical / population | `simulate_pulse_population`, `population_spec`, `fit_pulse_population`, `population_diagnostics` | Partial pooling; subject vs. population estimates; shrinkage. |
| F | `06-joint-driver-response` | ★ Flagship: coupled hormones | `simulate_pulse_joint`, `joint_spec`, `fit_pulse_joint`, `summary.joint_fit` | GnRH→LH / ACTH→cortisol; "is the response driven by the driver?" |
| G | `07-simulation-study` | Parameter recovery & study design | `simulate_pulse` + repeated `fit_pulse`, `effective_sample_size` | Recover truth; CI coverage; sampling-design / power angle. |
| H | `08-working-with-chains` | From fit to custom inference | `patient_chain`, `pulse_chain`, `chains`, `predict.pulse_fit`, `bp_predicted` | Tidy per-pulse posteriors; posterior-predictive checks. |

### Per-vignette detail

**A · `01-getting-started` — Single-subject deconvolution** *(salvages `my-vignette.Rmd`)*
Audience: both. Beats: what pulsatile hormones are → simulate one series → minimal spec → fit → read the four diagnostic plots → "you just deconvolved a hormone series." Keep salvaged Section-2 math; drop Rachel's TODO notes. Precompute: one converged single-subject fit (~25–50k iters, a few min offline).

**B · `02-anatomy-of-the-model` — The generative process**
Audience: both (intuition). Beats: decompose an observed series into baseline + Gaussian-pulse secretion + exponential elimination + error; small-multiples sweeping mass, width, half-life, inter-pulse interval. No MCMC. Precompute: none (instant; can run live).

**C · `03-priors-and-specification` — Specifying priors & starting values**
Audience: both, leans biostat. Beats: physiological meaning of each prior; informative vs. vague (prior-predictive plots); Strauss location prior (`gamma`/`range`); starting values & proposal variances; why reasonably informative priors matter (identifiability). Answers the WIP vignette's open "how much prior detail?" question. Precompute: light — prior-predictive draws; optionally one short prior-sensitivity fit.

**D · `04-convergence-diagnostics` — Assessing convergence**
Audience: biostatisticians. Beats: good vs. bad mixing (use the real half-life mixing issue as the honest worked example); multi-chain R-hat; ESS targets; remedies (longer runs, thinning, tighter priors). Reuses A's fit. Precompute: 2–4 chains at production length (minutes offline).

**E · `05-population-model` — Hierarchical / population modeling**
Audience: both. Beats: why pool subjects (partial pooling / borrowing strength); simulate a small cohort; fit; compare subject- vs. population-level posteriors; shrinkage visualization; cohort framing. Precompute: one population fit, ~5–10 subjects, production length (minutes offline).

**F · `06-joint-driver-response` — ★ Flagship: coupled hormones**
Audience: endocrinologists especially. Beats: the biology (GnRH→LH, ACTH→cortisol, LH→testosterone); simulate a coupled pair with known `rho`/`nu`; fit; interpret coupling strength (`rho`) and timescale (`nu`); headline inference — "is the response genuinely driven by the driver?" — with a null/uncoupled contrast (`rho≈0`) so readers see what *no* coupling looks like. Precompute: two joint fits (coupled + null), production length (minutes offline).

**G · `07-simulation-study` — Parameter recovery & study design**
Audience: biostatisticians + study designers. Beats: simulate-from-known-truth across replicates; bias and credible-interval coverage; then "how dense/long must sampling be to recover pulses reliably?" (vary `num_obs`/`interval`). Trust-building + power-analysis angle. Precompute: a modest replicate grid (heaviest offline job; cache results table as `.rds`).

**H · `08-working-with-chains` — From fit to custom inference**
Audience: power users / biostatisticians. Beats: extract per-pulse and per-parameter posterior draws as tidy data; roll your own summaries/plots; posterior-predictive checks. Reuses A's cached fit (no new MCMC). Precompute: none beyond reusing A.

### Ladder dependency graph

- **A** is the trunk.
- **D** and **H** reuse A's cached fit.
- **B** and **C** are conceptual prerequisites but computationally independent.
- **E** and **F** are the two "advanced model" branches.
- **G** stands alone.

---

## 5. Build sequence

1. **Infrastructure** (Section 3 wiring + precompute scaffolding).
2. **Core ladder:** A → B → C → D.
3. **Advanced + tooling:** E, F, G, H (any order; F prioritized as the flagship showcase).
4. **Roadmap stub:** `99-inference-roadmap.Rmd` (catalogue only).

---

## 6. Future appendix (I–O) — catalogued, NOT planned

Documented here and in a `99-inference-roadmap.Rmd` stub as short entries (one-line concept · gating backend dependency · `docs/lit_review.md` reference). No outlines, code, or estimates. Explicitly **"not scheduled; revisit after the current approach is solid."**

**Inference engines (speed/mixing) — gated on new C++/sampler work, cited to `docs/lit_review.md`:**
- **I** — SMC / generator-process model (HormoneBayes-style ON/OFF switching; eliminates transdimensional sampling).
- **J** — Compressed-sensing warm-start (L1/FOCUSS init to cut burn-in); HMC/NUTS; variational inference.

**Model extensions — gated on model + backend changes:**
- **K** — Multi-hormone networks (>2 nodes; e.g., GnRH→LH→testosterone cascade).
- **L** — Time-varying / circadian baseline `B_i(t)` (e.g., cortisol diurnal pattern).
- **M** — Covariate regression on pulse parameters (age / BMI / disease status).

**Inference tooling — gated on likelihood/IC support:**
- **N** — Model comparison / selection (WAIC / DIC across pulse-count or coupling models).

**Data — gated on the deferred "bundle real data" decision:**
- **O** — Real-data clinical case studies (PCOS, hypogonadism, Cushing's).

---

## 7. Out of scope / explicitly deferred

- Implementation or planning of I–O (revisit after the current approach is solid).
- Bundling real hormone datasets.
- Alternative samplers (SMC/HMC/VI) and the `docs/lit_review.md` research directions.

---

## 8. Update log — 2026-08-25 (post-merge with master)

The branch was cut before PRs #18 and #19 landed on master. Merging them in
invalidated several assumptions recorded above. Amendments, in order of impact:

**The parameterization changed, and with it every prior number in every
vignette.** `pulse_spec()` now defaults to `pulse_distribution = "lognormal"`,
`sd_prior = "uniform"`, `student_t_pulses = FALSE`. Pulse mass and width priors
are on the **log** scale (width additionally on the log of a *variance*). All of
A, D, and the real-data vignette were passing natural-scale numbers and were
silently requesting priors an order of magnitude off. They now set only
data-specific arguments (baseline, half-life, pulse count, error, location prior)
and leave mass/width to the tuned defaults.

**`simulate_pulse()` and `pulse_spec()` disagree by default.** The simulator
still defaults to `"truncnorm"` while the spec defaults to `"lognormal"`, so any
simulate-then-fit vignette must pass `pulse_distribution = "lognormal"` to the
simulator explicitly. A and D now do. B deliberately stays on the natural-scale
simulator, since it never fits a model and natural scale reads more clearly.

**The 27-minute runtime constraint in §2 is obsolete.** A 250,000-iteration
single-subject fit now takes about 50 seconds. The precompute pattern is still
worth keeping — CRAN build limits are unchanged and per-vignette build time
should stay near zero — but "BD-MCMC is slow" is no longer the binding
constraint it was when this was written, and vignettes can afford production
chain lengths freely.

**The mixing story moved.** §2 and vignette D were written around half-life as
the problem parameter, with the pulse-width SD frozen. #18 un-froze the SDs, so
`width_sd` is now estimated and reported. At production length the weak parameters
are **half-life and baseline** (they trade off against each other in the troughs);
mass and width mix well. D was rewritten accordingly and now also demonstrates
`identifiability_check()`.

**Real data needs informative priors; simulated data does not.** The real-data
vignette will not fit sensibly on package defaults — width collapses below the
sampling resolution and half-life ESS falls into the low teens. It now sets a
tight width prior and a tighter half-life prior, and says plainly why. This is
the §2 "weakly-identified parameterization" bullet showing up in practice, and it
is arguably the most useful thing that vignette teaches.

**Decision reversed: a real-data vignette exists.** §1 recorded a
simulation-only decision with real data deferred. `v99-real-data-lh` fits
`datasets::lh` (Diggle's 48-sample series). No dataset is bundled — it comes from
base R's **datasets** package — so the "don't ship a dataset" part of the
decision stands; the "simulation only" part does not. It is numbered `v99` to sit
outside the A–H ladder. Item **O** in §6 (real-data clinical case studies) is
still deferred; this does not satisfy it.

**Vignette numbering.** Files are `v01`–`v04` and `v99`, not the `01`–`08` in the
table above. E–H remain unwritten.

### Resolved: observations logged against the pre-merge branch

Two problems recorded while reviewing the original rendered vignettes:

1. *"`width_sd` is fixed at 10.000, no variation at all."* — Real bug. Fixed by
   #18, which un-freezes and estimates the pulse-to-pulse SDs. Verified after the
   merge: `width_sd` now varies across draws with a healthy ESS. The claim in A's
   prose that the width SD "is held fixed in this version of the package" has been
   removed.
2. *"Pulse width appears to not do anything"* (vignette B). — **Not a bug.**
   `width_mean` is on the variance scale, so B's sweep of 15/35/60 was really a
   sweep of pulse SD 3.9/5.9/7.7 minutes, all below the 10-minute sampling
   interval and therefore genuinely indistinguishable in the output. B now sweeps
   4/100/625 (SD 2/10/25 minutes) and explains the variance scale. The underlying
   point — that pulses narrower than the sampling grid cannot be resolved — was
   worth surfacing and is now made explicitly in B, C, and the real-data vignette.

### Deferred package fix (logged 2026-06-07)

Extend the convergence diagnostics — `plot_trace()`, `plot_acf()`, `convergence_report()`, and `population_diagnostics()` — to also accept **single-subject `pulse_fit`** objects. They are currently population-only (they `stop()` on non-`population_fit` and read `fit$population_chain` / `fit$num_subjects`). This was discovered while planning vignette D, which works around the gap by using only the model-agnostic tools (`effective_sample_size`, `gelman_rubin`, `bp_trace`). When fixed, vignette D can adopt the turnkey suite directly. Likely approach: a `pulse_fit` branch that pulls the chain from `patient_chain(fit)` with `num_subjects = 1`, plus tests and Roxygen updates. Not scheduled.
