# Annotated bibliography: alternatives and extensions to BD-MCMC deconvolution

**Date:** 2026-06-07
**Status:** Background reading behind the future-appendix entries (I-O) in
`2026-06-07-vignette-roadmap-design.md`. Method summaries only; nothing here
is scheduled or endorsed work.
**Citations:** All seven DOIs verified against Crossref (2026-08-25);
first author, year, and title match in each case.

---

## 1. HormoneBayes: SMC in place of BD-MCMC (2024)

Voliotis, Abbara, Prague, Veldhuis, et al. avoid transdimensional sampling by
reformulating the generative model [1]:

- Instead of estimating the number/location of pulses (requiring BD-MCMC), the
  underlying hypothalamic ON/OFF switching process is modeled as a two-state
  Markov chain with transition rates tau_ON and tau_OFF; the pulse count
  emerges from the ON/OFF process rather than being a sampled discrete
  parameter.
- The observed LH is driven by both a pulsatile signal P_t and a basal signal
  B_t.
- Inference uses a Gibbs sampler where:
  - Sequential Monte Carlo (SMC) with ancestral sampling handles the latent
    state trajectories (H_t, B_t)
  - Simplified manifold MALA (sMMALA) samples the continuous parameters
    (k, d, f) using gradient information
  - Adaptive Metropolis-Hastings handles (tau_ON, tau_OFF)

The model uses 5 parameters plus latent states (vs. a per-pulse parameter
structure). The authors report convergence across healthy men,
pre/post-menopausal women, PCOS, and HA cohorts; an open-source C++
implementation is available. They cite Johnson et al. — the BD-MCMC approach
libpulsatile implements — as prior art.

## 2. Compressed sensing / sparse optimization (2014, 2022)

Faghih, Dahleh, Brown et al. reformulate pulse deconvolution as a sparse
recovery problem without MCMC [2]:

- Cortisol secretion is modeled as a 2nd-order linear ODE with sparse pulsatile
  inputs.
- Solved via coordinate descent: the FOCUSS algorithm (iteratively reweighted
  L1 minimization) plus generalized cross-validation for regularization.
- Assumes 15-22 secretory events over 24 hours; GCV balances sparsity against
  residual error.
- Reports R^2 > 0.92 on real 24-hour cortisol data sampled every 10 minutes.
- Produces point estimates only; no posterior uncertainty.

Amin, Faghih et al. (2022) extended this to multi-hormone sparse system
identification (leptin-cortisol dynamics) using state-space models with sparse
recovery [3].

## 3. Variational inference for point-process deconvolution (2020)

Shibue & Komaki use marked point processes as latent variables for calcium
imaging deconvolution (structurally analogous to hormone-pulse deconvolution)
and solve with variational inference rather than MCMC [4], simultaneously
estimating event times, amplitudes, and shapes. VI yields approximate
posteriors at reported MCMC-competitive quality with lower compute cost.

## 4. HMC for birth-death-sampling models (2024)

Shao, Magee, Suchard developed scalable gradients enabling Hamiltonian Monte
Carlo for episodic birth-death-sampling models, reporting 10-200x efficiency
gains over standard MCMC [5], based on a linear-time gradient computation.
The application is phylodynamics, but the mathematical structure (birth-death
process with observations) is closely related.

## 5. Normalizing flows to accelerate MCMC (2022)

Gabrie, Rotskoff, Vanden-Eijnden formalize Monte Carlo augmented with
normalizing flows [6]: a normalizing flow is learned as an adaptive proposal
distribution for MCMC, targeting the posterior's modes, with limited prior
data needed for training.

## 6. Posterior-based proposals (2019)

Pooley et al. introduce posterior-based proposals (PBPs) for accelerating MCMC,
reporting them "significantly faster than or competitive with existing
methods" across various model types [7]. Information from the posterior
informs the proposal distribution.

---

## References

[1] Voliotis M, Abbara A, Prague JK, Veldhuis JD, Dhillo WS,
Tsaneva-Atanasova K. "HormoneBayes: A novel Bayesian framework for the
analysis of pulsatile hormone dynamics." PLOS Computational Biology 20,
e1011928 (2024). doi:10.1371/journal.pcbi.1011928

[2] Faghih RT, Dahleh MA, Adler GK, Klerman EB, Brown EN. "Deconvolution of
Serum Cortisol Levels by Using Compressed Sensing." PLoS ONE 9, e85204 (2014).
doi:10.1371/journal.pone.0085204

[3] Amin MR, Pednekar DD, Azgomi HF, van Wietmarschen H, Aschbacher K,
Faghih RT. "Sparse System Identification of Leptin Dynamics in Women With
Obesity." Frontiers in Endocrinology 13, 769951 (2022).
doi:10.3389/fendo.2022.769951

[4] Shibue R, Komaki F. "Deconvolution of calcium imaging data using marked
point processes." PLoS Computational Biology 16, e1007650 (2020).
doi:10.1371/journal.pcbi.1007650

[5] Shao Y, Magee AF, Vasylyeva TI, Suchard MA. "Scalable gradients enable
Hamiltonian Monte Carlo sampling for phylodynamic inference under episodic
birth-death-sampling models." PLOS Computational Biology 20, e1011640 (2024).
doi:10.1371/journal.pcbi.1011640

[6] Gabrie M, Rotskoff GM, Vanden-Eijnden E. "Adaptive Monte Carlo augmented
with normalizing flows." PNAS 119, e2109420119 (2022).
doi:10.1073/pnas.2109420119

[7] Pooley CM, Bishop SC, Doeschl-Wilson A, Marion G. "Posterior-based
proposals for speeding up Markov chain Monte Carlo." Royal Society Open
Science 6, 190619 (2019). doi:10.1098/rsos.190619
