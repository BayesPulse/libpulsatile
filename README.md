# libpulsatile

[![Build Status](https://travis-ci.org/BayesPulse/libpulsatile.svg?branch=master)](https://travis-ci.org/BayesPulse/libpulsatile)
[![codecov](https://codecov.io/gh/BayesPulse/libpulsatile/branch/master/graph/badge.svg)](https://codecov.io/gh/BayesPulse/libpulsatile)


The development repository for the C++ backend to the unified Bayesian pulsatile
hormone modeling algorithm used in the poppulsatile R package.

This library extends (and refactors) the single-subject model to multiple
subjects and up to two-hormones per patient (driver and response hormones).

The package is currently in development and not yet functional, but feel free to
take a look around.

---

**Installation and dependencies**

This library is intended to be used from within an R package, and therefore
relies on R's RNGs and datastructures.  The list of dependencies for the
development phase is rather specific and will be reduced later.

- C++11
- clang4.0 (on Mac for OpenMP support)
- R (>= 3.4.3)
- R packages
  - RInside (development version via GitHub)
  - Rcpp (CRAN)
  - RcppArmadillo (CRAN)

Once these dependencies are satisfied, build the executable and the test
executable with `make`.  To run the tests, run `./bin/tests` or `./bin/tests -s`
to view the detailed test results.

