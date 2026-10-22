# Sampler performance: running-sum likelihood

Reproducible script: `benchmarks/sampler_equivalence.R`.

## Summary

Each likelihood evaluation used to re-sum the contributions of all K pulses
over all n observations: O(K·n) per Metropolis-Hastings proposal and
O(K²·n) per birth-death pass (K partial likelihoods). Birth-death accounted
for ~75% of run time. The sampler now keeps a running sum of pulse
contributions and updates it incrementally (see `Patient::sync_pulse_sum()`),
which gives **3.3x-13x faster fits with unchanged chains**. The speedup grows
with series length because the removed term scaled as K²·n.

## Correctness: the chains are unchanged

The reference build is the original algorithm with only the two bug fixes on
this branch (initial-pulse contribution, RNG replay -- see below), so it draws
the same random stream. On identical seeded fits the largest relative
difference anywhere in any chain is at floating-point rounding level: every
accept/reject and birth/death decision is identical and the chains agree draw
for draw.

Additional checks:

- A debug build (`-DBP_VERIFY_LIKELIHOOD`) compares every fast-path likelihood
  and partial likelihood with the original full recomputation during real
  fits of all three models and stops on any mismatch beyond 1e-9 relative;
  none occurred. Its chains are bit-identical to the normal build's.
- Catch tests (`tests/likelihood_cache_tests.cpp`) drive a 400-step random
  sequence of proposals, restores, decay changes, births, deaths and
  error-variance changes, checking every quantity against the reference after
  each step, plus exact restores, Patient copies, and the slot cache.

## Results

Linux x86-64 container, 4 cores, R 4.3.3, `-O2`. Best of 3 runs;
`thin = 1`, `burnin = 0`. `max_rel_diff` is the largest
`|a - b| / (1 + |a|)` over every numeric chain.

Against the reference (original algorithm + bug fixes):

| scenario                     | iters | reference (s) | new (s) | speedup | max_rel_diff |
|------------------------------|------:|--------------:|--------:|--------:|-------------:|
| single subject, 144 obs      |  2000 |          3.30 |    0.69 |    4.8x |      2.6e-17 |
| single subject, 288 obs      |  2000 |         20.95 |    2.72 |    7.7x |      1.3e-17 |
| single subject, 576 obs      |   500 |         44.53 |    3.43 |   13.0x |      2.0e-17 |
| population, 3 x 144 obs      |  2000 |          9.97 |    3.00 |    3.3x |      6.7e-18 |
| joint driver-response, 144   |  2000 |          3.57 |    1.03 |    3.5x |      2.0e-17 |

Against `master` before this branch (speed only -- the bug fixes change the
random stream, so the chains legitimately differ):

| scenario                     | iters | master (s) | new (s) | speedup |
|------------------------------|------:|-----------:|--------:|--------:|
| single subject, 144 obs      |  2000 |       3.34 |    0.62 |    5.4x |
| single subject, 288 obs      |  2000 |      21.62 |    2.63 |    8.2x |
| single subject, 576 obs      |   500 |      49.92 |    3.27 |   15.3x |
| population, 3 x 144 obs      |  2000 |      12.68 |    3.36 |    3.8x |
| joint driver-response, 144   |  2000 |       4.11 |    1.12 |    3.7x |

## Where the time went

Instruction counts (callgrind, single subject, 288 obs, 150 iterations,
sampler only):

| step                                                        | instructions | vs. original |
|-------------------------------------------------------------|-------------:|-------------:|
| original                                                    |       16.2 B |         1.0x |
| running sum, two-slot pulse cache                           |        4.9 B |         3.3x |
| reuse birth-death rates after a rejected birth              |        3.4 B |         4.7x |
| single-pass contribution kernel                             |        3.1 B |         5.2x |
| sparse partials, `pnorm()` short-circuit, undo of rejected  |        2.1 B |         7.8x |
| proposals with memoized likelihood                          |              |              |

(The RNG fix below also removes ~180 `GetRNGstate()`/`PutRNGstate()` pairs
per iteration, about 0.4 ms.)

What remains is mostly per-evaluation arithmetic: `log()` in the likelihood
(~40%) and `pnorm()`/`exp()` for newly proposed pulse parameters (~20%).

## Bug fixes found along the way

Both change the random stream, so they are separate commits and are included
in the reference build above.

1. **Initial pulse contribution** -- the `Patient` constructor evaluated the
   initial pulse's contribution at the log-concentrations instead of the
   observation times. Because the contribution cache is keyed only on the
   pulse parameters, the wrong vector persisted until the first baseline/
   half-life or pulse proposal, so the first birth-death step of every
   single-subject fit used an incorrect likelihood.
2. **Replayed random numbers** -- each sampler opened an `Rcpp::RNGScope`,
   and nested scopes re-sync the R RNG state (`GetRNGstate()` reloads the
   state saved by the last `PutRNGstate()`). Draws made outside a scope were
   therefore rewound and drawn again by the next sampler: the single-subject
   error-variance Gibbs draw shared its uniforms with the next birth-death
   step, and the population mean/error Gibbs draws with the following MH
   proposals. `BPRNGScope` syncs only at the outermost level.
