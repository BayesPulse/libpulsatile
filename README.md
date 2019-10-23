# libpulsatile

[![Build Status](https://travis-ci.org/BayesPulse/libpulsatile.svg?branch=master)](https://travis-ci.org/BayesPulse/libpulsatile)
[![codecov](https://codecov.io/gh/BayesPulse/libpulsatile/branch/master/graph/badge.svg)](https://codecov.io/gh/BayesPulse/libpulsatile)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version-ago/BayesPulse)](https://cran.r-project.org/package=BayesPulse) 
[![CRAN\_Download\_Badge](http://cranlogs.r-pkg.org/badges/BayesPulse)](https://cran.r-project.org/package=BayesPulse)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)<Paste>
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

The development repository for the C++ backend to the unified Bayesian pulsatile
hormone modeling algorithm used in the bayespulse R package.

This library extends (and refactors) the single-subject model to multiple
subjects and up to two-hormones per patient (driver and response hormones).

The package is currently in development, with only the single subject model currently functioning. Feel free to
take a look around and test it out. 

---

## Dependencies

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


## Important Note for Cloning this repo

This project structure relies on symbolic links.  If you are working on Linux or
Mac then you should have not problem cloning and working with this repo.  If,
however, you are working on a Windows machine there are a few additional steps
you will need to take.

1. Windows Vista or newer with NTFS file system, not FAT.
2. You need administrator rights or at least `SeCreateSymbolicLinkPrivilege`
   privilege
3. git bash version 2.10.2 or later.  It will be helpful to install with
   `core.symlinks` option turned on.

In the git bash shell clone the repo.  (The example below uses SSH, change the
URL as needed for https.)

    git clone -c core.symlinks=true git@github.com:BayesPulse/libpulsatile.git

If you are using Windows and you are not sure that you cloned the repo as noted
above then please reclone the repo!


## Installing R package

**From the command line**
After cloning the repository navigate to the `R-package` directory and use the
provided `Makefile`. Note that `gawk` is required on MacOS.


```{sh}
cd <path-to-libpulsatile>/R-package
make install
``` 

**Within an Interactive R Session**

```{r}
setwd("<path-to-libpulsatile>/R-package")
devtools::load_all(recompile = TRUE)
devtools::document()
devtools::install()
```

**Using RStudio**
Either use the instructions for installing the package in an interactive R
session or open the project via the file `bayespulse.Rproj` and use the build
and install menu options.  The project file should tell RStudio to use the
`Makefile`.

## Building the C++ binaries

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

When compiled outside of the R package, the library is dependent on RInside for
random number generators. The chains class also uses R object types and is
geared towards structuring the R return object. Gist is, use the R package and C++
only is for dev work.

