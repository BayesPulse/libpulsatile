# libpulsatile

[![Build Status](https://travis-ci.org/BayesPulse/libpulsatile.svg?branch=master)](https://travis-ci.org/BayesPulse/libpulsatile)
[![codecov](https://codecov.io/gh/BayesPulse/libpulsatile/branch/master/graph/badge.svg)](https://codecov.io/gh/BayesPulse/libpulsatile)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version-ago/BayesPulse)](https://cran.r-project.org/package=BayesPulse) 
[![CRAN\_Download\_Badge](http://cranlogs.r-pkg.org/badges/BayesPulse)](https://cran.r-project.org/package=BayesPulse)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)<Paste>
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

The development repository for the C++ backend to the unified Bayesian pulsatile
hormone modeling algorithm used in the poppulsatile R package.

This library extends (and refactors) the single-subject model to multiple
subjects and up to two-hormones per patient (driver and response hormones).

The package is currently in development and not yet functional, but feel free to
take a look around.

---

**Dependencies**

This library is intended to be used from within an R package, and therefore
relies on R's RNGs and datastructures.  The list of dependencies for the
development phase is rather specific and will be reduced later.

- C++11
- clang (>= 4.0) (on Mac for OpenMP support)
  - libclang-x.x-dev (ubuntu)
  - libomp-dev (Ubuntu)
- R (>= 3.4.3)
- R packages
  - [Rcpp](https://cran.r-project.org/package=Rcpp)
  - [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo)
  - [RInside](https://cran.r-project.org/package=RInside)

**Installing R package**

The process for compiling and installing the R package is a bit more complicated
than normal due to the repository focusin on C++ code, with the R package
structure in a subdirectory.

First clone the repository and enter the R-package directory

```{sh}
git clone git@github.com:BayesPulse/libpulsatile.git
cd libpulsatile/R-package
```

Then open R, run the build_package.R script, and install with devtools

```{R}
source("./build_package.R", echo=TRUE)
Rcpp::compileAttributes()
devtools::document()
devtools::install()
```

**Building the C++ binaries**

*This is only necessary if you are running the binary directly.* To build the C++
binary, first ensure the dependencies are satisfied, then build the executable
and the test executable with `make`.  To run the tests, run `./bin/tests` or
`./bin/tests -s` to view the detailed test results.

```{sh}
make clean
make
./bin/tests -s
```


# Notes

When compiled outside of the R package, the library is generally dependent on RInside for
random number generators. The chains class also uses R object types and is
geared towards structuring the R return object. Gist is, use the R package and C++
only is for dev work.


