# bayespulse

[![CI (hosted)](https://github.com/mmulvahill/libpulsatile/actions/workflows/ci-hosted.yml/badge.svg?branch=master)](https://github.com/mmulvahill/libpulsatile/actions/workflows/ci-hosted.yml)
[![CI (self-hosted)](https://github.com/mmulvahill/libpulsatile/actions/workflows/ci-selfhosted.yml/badge.svg?branch=master)](https://github.com/mmulvahill/libpulsatile/actions/workflows/ci-selfhosted.yml)
[![Coverage](https://codecov.io/gh/mmulvahill/libpulsatile/branch/master/graph/badge.svg)](https://codecov.io/gh/mmulvahill/libpulsatile)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle: Maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)

**Bayesian deconvolution of pulsatile hormone concentration data**

A unified Bayesian framework for analyzing pulsatile hormone secretion patterns in both single-subject and population studies. The package implements birth-death MCMC for trans-dimensional pulse detection and hierarchical models for multi-subject analysis.

## Features

- **Single Subject Model**: Deconvolution analysis for individual hormone time series
- **Population Model**: Hierarchical Bayesian model for analyzing multiple subjects simultaneously
  - Population-level parameters (means and variances)
  - Subject-specific parameters with shrinkage
  - Individual pulse characteristics
- **Joint Hormone Model**: Coupled analysis of driver-response hormone pairs
  - Driver pulses influence response pulse occurrence
  - Coupling strength (ρ) and decay rate (ν) estimation
  - Simultaneous deconvolution of both hormones
- **MCMC Diagnostics**: Comprehensive convergence assessment tools
  - Effective Sample Size (ESS) calculation
  - Gelman-Rubin diagnostics
  - Trace plots and autocorrelation functions
- **Performance Benchmarking**: Linear scaling with subjects and iterations
- **Flexible Priors**: Customizable priors for all model parameters

## Installation

### Prerequisites

**R packages:**
```r
install.packages(c("Rcpp", "RcppArmadillo", "dplyr", "ggplot2",
                   "tibble", "tidyr", "rlang", "devtools"))
```

**System libraries:**
- **gfortran** (required for RcppArmadillo)
  - Ubuntu/Debian: `sudo apt-get install gfortran`
  - macOS: `brew install gfortran`
  - For macOS troubleshooting: https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/

### Installing the Package

**From source (recommended for development):**
```r
# Clone the repository (see Platform Support below for Windows)
git clone git@github.com:mmulvahill/libpulsatile.git
cd libpulsatile/R-package

# Install using devtools
R -e "devtools::install()"
```

**Using the Makefile:**
```bash
cd libpulsatile/R-package
make install
```

## Quick Start

### Single Subject Analysis

```r
library(bayespulse)

# Simulate pulsatile hormone data
sim_data <- simulate_pulse(
  num_obs = 24,
  interval = 10,
  mass_mean = 3.5,
  width_mean = 35
)

# Fit the model
fit <- fit_pulse(
  data = sim_data$data,
  iters = 10000,
  thin = 10,
  burnin = 5000
)

# View results
summary(fit)
plot(fit)
```

### Population Analysis

```r
# Simulate data for multiple subjects
sim_data <- lapply(1:5, function(i) {
  simulate_pulse(num_obs = 24, interval = 10)
})

# Create model specification
spec <- population_spec(
  prior_mass_mean_mean = 3.5,
  prior_width_mean_mean = 42,
  prior_mean_pulse_count = 12
)

# Fit population model
fit <- fit_pulse_population(
  data = sim_data,
  spec = spec,
  iters = 10000,
  thin = 10,
  burnin = 5000
)

# Examine convergence
summary(fit)
convergence_report(fit)
plot_trace(fit)
plot_acf(fit)
```

## Performance

Runtime scales roughly linearly in iterations and observations. Measured
on Apple Silicon (2026-08), at the production defaults of 250,000
iterations:

- **Single subject** (60 observations): under a minute
- **Population, 5 subjects** (720 total observations): ~11 minutes
- **Joint driver-response** (2 x 144 observations): a few minutes

Benchmark scripts live in [R-package/benchmarks/](R-package/benchmarks/);
the simulation-study vignette (`v07`) shows how to time and scale fits for
your own design.

## Platform Support

### Supported Platforms
- **Linux**: Fully supported (Ubuntu 20.04+, tested in CI/CD)
- **macOS**: Fully supported (Intel and Apple Silicon, tested in CI/CD)

### Windows Support
**Status**: Not currently supported

Windows support will be required for CRAN submission and is planned for future development. The primary challenge is the symbolic link structure used in the dual-build system (C++ library + R package).

**For Windows users attempting to build from source:**
The project uses symbolic links between the C++ library and R package. You'll need:
1. Windows Vista+ with NTFS file system
2. Administrator rights or `SeCreateSymbolicLinkPrivilege`
3. Git Bash 2.10.2+
4. Clone with: `git clone -c core.symlinks=true git@github.com:mmulvahill/libpulsatile.git`

Note: Even with proper symbolic link setup, Windows builds are untested and may encounter other issues.

## Development

This package has a dual-build system:
- **R package** (`R-package/`): Primary user interface (production)
- **C++ library** (root): Standalone build for C++ development and testing

### For Developers

See [CLAUDE.md](CLAUDE.md) for comprehensive development guidelines including:
- Code style and documentation standards
- Pre-commit quality assurance checklist
- Common issues and troubleshooting
- CI/CD pipeline details

**Quick pre-commit check:**
```bash
cd R-package
Rscript -e "Rcpp::compileAttributes(); devtools::document()"
Rscript -e "devtools::test()"
Rscript -e "devtools::check()"  # Must pass with 0 errors, 0 warnings
```

### Running Tests

**R tests:**
```bash
cd R-package
Rscript -e "devtools::test()"
```

**C++ tests:**
```bash
make clean && make
./bin/tests
```

## Citation

If you use this software in your research, please cite:

```bibtex
@software{bayespulse,
  title = {bayespulse: Bayesian Analysis of Pulsatile Hormone Data},
  author = {Mulvahill, Matthew and Carlson, Nichole},
  year = {2025},
  url = {https://github.com/mmulvahill/libpulsatile}
}
```

## License

GPL (>= 2)

## Authors

- **Matthew Mulvahill** - *Creator and Maintainer* - matt@mulvahill.com
- **Nichole Carlson** - *Author and Thesis Supervisor*

## Project Status

Active development with both single-subject and population models fully implemented. The package is undergoing final testing and documentation refinement before CRAN submission.

### Recent Updates
- ✅ Population model implementation with hierarchical structure
- ✅ Comprehensive MCMC convergence diagnostics
- ✅ Performance benchmarking and optimization
- ✅ Extensive test coverage
- 🔄 CRAN submission preparation (in progress)
- 📋 Windows platform support (planned)

## References

For more information on the statistical methodology, see:
- Carlson, N.E. (2018). "Bayesian Deconvolution of Pulsatile Hormone Concentration Data"
- Johnson, M.L. et al. (2008). "Deconvolution analysis of hormone data"

## Support

- **Issues**: https://github.com/mmulvahill/libpulsatile/issues
- **Documentation**: Run `?fit_pulse` or `?fit_pulse_population` in R
- **Developer Guide**: See [CLAUDE.md](CLAUDE.md)
