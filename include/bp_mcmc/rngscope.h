#ifndef GUARD_bp_mcmc_rngscope_h
#define GUARD_bp_mcmc_rngscope_h

#include <Rcpp.h>

//
// rngscope.h
//   BPRNGScope: a depth-counted replacement for Rcpp::RNGScope.
//
//   Rcpp::RNGScope calls GetRNGstate()/PutRNGstate() (copying .Random.seed in
//   and out of R) every time one is constructed and destroyed, even when
//   nested. That had two consequences for the samplers, which open a scope on
//   every draw:
//
//   1. Cost: roughly 2us per scope, a large share of an MCMC iteration.
//   2. Reused random numbers: GetRNGstate() reloads the state saved by the
//      last PutRNGstate(), so any draw made *between* nested scopes (e.g. the
//      error-variance Gibbs draw, which opens no scope) was rewound and its
//      uniforms were drawn again by the next sampler that opened one.
//
//   BPRNGScope only syncs the R RNG state at the outermost level, so nested
//   scopes are free and every draw advances one continuous stream. The model
//   entry points (singlesubject_(), population_(), jointsinglesubject_())
//   open one for the whole run; samplers called directly (e.g. from the C++
//   tests) still sync on every call.
//

inline int & bp_rng_scope_depth() {
  static int depth = 0;
  return depth;
}

class BPRNGScope {
  public:
    BPRNGScope()  { if (bp_rng_scope_depth()++ == 0) GetRNGstate(); }
    ~BPRNGScope() { if (--bp_rng_scope_depth() == 0) PutRNGstate(); }
    BPRNGScope(const BPRNGScope &) = delete;
    BPRNGScope & operator=(const BPRNGScope &) = delete;
};

#endif
