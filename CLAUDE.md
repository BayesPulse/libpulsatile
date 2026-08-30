# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Architecture

**libpulsatile** is a C++ library with R package bindings for Bayesian pulsatile hormone modeling. The project has a dual-build system:

1. **Standalone C++ library** (root directory) - For development and testing C++ components
2. **R package** (`R-package/` directory) - Production interface using Rcpp/RcppArmadillo

### Core Components
- **Data Structures** (`include/bp_datastructures/`): Patient, Population, Chains, PulseEstimates
- **MCMC Algorithms** (`include/bp_mcmc/`): Metropolis-Hastings, Birth-Death process, Proposal variance adaptation
- **Single Subject Model** (`include/bpmod_singlesubject/`): Deconvolution for individual subjects
- **Population Model** (`include/bpmod_population/`): Hierarchical model for multiple subjects
- **R Interface** (`R-package/src/`): Rcpp exports including `singlesubject_()` and `population_()`
  - `fit_pulse()`: Single subject model fitting
  - `fit_pulse_population()`: Population model fitting
  - `population_spec()`: Population model specification
  - `population_diagnostics()`: MCMC convergence diagnostics
  - `convergence_report()`, `plot_trace()`, `plot_acf()`: Diagnostic visualization

### Symbolic Links
The project uses symbolic links between the main library and R package. Windows users must clone with `git clone -c core.symlinks=true` and have appropriate permissions.

## Dependencies
- **R packages**:
  - Required: Rcpp, RcppArmadillo, dplyr, ggplot2, tibble, tidyr, rlang
  - Development: testthat, remotes, magrittr, devtools
  - Benchmarking: pryr (for memory profiling)
  - Optional: RInside (for standalone C++ builds)
- **System libraries**: gfortran (required for RcppArmadillo)
- Install R dependencies: `Rscript -e "install.packages(c('Rcpp', 'RcppArmadillo', 'dplyr', 'ggplot2', 'tibble', 'tidyr', 'rlang', 'testthat', 'remotes', 'magrittr', 'devtools', 'pryr'))"`

## Build and Test Commands
- **Build library**: `make`
- **Build R package**: `cd R-package && make`
- **Run C++ tests**: `bin/tests` or single test: `bin/tests [test-name]`
- **Run R tests**: `cd R-package && Rscript -e 'devtools::test()'`
- **Check R package**: `cd R-package && Rscript -e 'devtools::check()'`
- **Install R package**: `cd R-package && Rscript -e 'devtools::install()'`
- **Generate Rcpp exports**: `cd R-package && Rscript -e "Rcpp::compileAttributes()"`
- **Rebuild documentation**: `cd R-package && Rscript -e "devtools::document()"`

## Quality Assurance / Pre-commit Checklist

**CRITICAL: Always run these checks locally BEFORE pushing to GitHub to avoid CI/CD failures**

### Before Every Commit to R Package Code:

1. **Rebuild documentation** (especially after modifying Roxygen comments):
   ```bash
   cd R-package
   Rscript -e "Rcpp::compileAttributes(); devtools::document()"
   ```

2. **Run R CMD check** (this is what CI/CD runs - catch issues early):
   ```bash
   cd R-package
   Rscript -e "devtools::check()"
   ```
   - **Must complete with 0 ERRORs, 0 WARNINGs**
   - NOTEs are acceptable but should be reviewed
   - Common issues to watch for:
     - Undocumented function parameters
     - Missing `@importFrom` declarations
     - Non-ASCII characters (use `\uXXXX` Unicode escapes for Greek letters)
     - Empty `@section` tags in Roxygen documentation
     - NSE variable bindings (use `utils::globalVariables()`)

3. **Run all tests**:
   ```bash
   cd R-package
   Rscript -e "devtools::test()"
   ```
   - All tests must pass
   - Add new tests for new functionality

4. **Build and install** (verify package builds cleanly):
   ```bash
   cd R-package
   Rscript -e "devtools::install()"
   ```

### Before Creating a Pull Request:

1. Complete all steps above
2. Review all modified files for:
   - Code quality and adherence to style guide
   - Appropriate documentation (Roxygen2 for R, comments for C++)
   - Test coverage for new functionality
3. Ensure commit messages are descriptive
4. Verify CI/CD will pass by running local equivalents

### Quick Pre-commit Command:
```bash
cd R-package && \
  Rscript -e "Rcpp::compileAttributes(); devtools::document()" && \
  Rscript -e "devtools::test()" && \
  Rscript -e "devtools::check()"
```

This comprehensive check typically takes 2-5 minutes but saves 30+ minutes of CI/CD debugging.

## Code Style Guidelines

### General Style
- **Headers**: Use include guards: `GUARD_filename_h`
- **Naming**: ClassNames, function_names, local_variables, CONSTANTS
- **Formatting**: 2-space indentation, 80-character line length
- **C++ Standards**: C++14 (required by current RcppArmadillo; set via `-std=c++14` in the Makefile and `CXX_STD = CXX14` in `R-package/src/Makevars`), use Rcpp/RcppArmadillo for R integration
- **R Style**: Follow tidyverse style guide, use `<-` for assignment
- **Imports**: Group imports by library (Rcpp, std, project headers)

### Documentation
- **R functions**: Use Roxygen2 with `#'` comments
  - Document ALL parameters with `@param`
  - Include `@return` for return values
  - Use `@examples` for usage demonstrations
  - Add `@importFrom` for all external functions used
  - **Never** leave empty `@section` tags - either fill them or remove them
- **C++ code**: Header comments explaining purpose and algorithm
- **Inline comments**: Explain complex logic, not obvious operations

### CRAN Compliance (Critical for R Package)
- **Character encoding**: Use ONLY ASCII in R code and NAMESPACE
  - Greek letters: Use Unicode escapes (μ → `\u03BC`, α → `\u03B1`, ω → `\u03C9`, etc.)
  - Find non-ASCII: `tools::showNonASCIIfile("R/your_file.R")`
- **NSE variables**: Declare global variables for ggplot2/dplyr NSE
  ```r
  utils::globalVariables(c("iteration", "value", "lag"))
  ```
- **Imports**: Explicitly declare ALL imported functions
  - Add `@importFrom package function` to Roxygen documentation
  - Never rely on implicit imports from dependencies
- **No warnings**: R CMD check must complete with 0 WARNINGs, 0 ERRORs

### Error Handling
- **C++**: `Rcpp::stop()` for errors with informative messages
- **R**: `stopifnot()` for input validation, `stop()` for runtime errors
- **Tests**: Catch2 v2.13.10 framework for C++, testthat for R

## Development Workflow

### General Workflow
1. Make code changes following style guidelines
2. Add/update tests for new functionality
3. Run pre-commit checks (see Quality Assurance section above)
4. Commit with descriptive messages
5. Push only after local checks pass

### Model Architecture
- **Single subject model**: `singlesubject_()` Rcpp export in `R-package/src/singlesubject.cpp`
  - C++ implementation in `include/bpmod_singlesubject/`
  - R interface: `fit_pulse()`
- **Population model**: `population_()` Rcpp export in `R-package/src/population.cpp`
  - C++ implementation in `include/bpmod_population/`
  - R interface: `fit_pulse_population()`, `population_spec()`
  - Diagnostics: `population_diagnostics()`, `convergence_report()`

### Adding New Features
- C++ development should follow the existing modular structure in `include/` directories
- When adding new functionality:
  1. Add corresponding tests in both `tests/` (C++) and `R-package/tests/testthat/` (R)
  2. Update Roxygen2 documentation for any R-visible functions
  3. Regenerate documentation: `Rcpp::compileAttributes(); devtools::document()`
  4. Run R CMD check to verify no warnings introduced
- Use `make clean && make` to ensure clean builds after header changes

### Benchmarking New Features
- Performance benchmarks go in `R-package/benchmarks/` (excluded from package builds via `.Rbuildignore`)
- Use `pryr::mem_used()` for memory profiling
- Document benchmark results in markdown reports
- Track performance metrics over time for regression detection

## CI/CD Pipeline

**GitHub Actions** workflows in `.github/workflows/` run on every push and PR:

### Workflows
- **`ci-hosted.yml`**: Build and test on GitHub-hosted runners with the
  current R release (`r-version: 'release'`; the old 4.3.0 pin had aged out
  of CRAN's binary windows, forcing slow source installs). R dependencies are
  cached (`actions/cache` on the user library, keyed on resolved R version +
  DESCRIPTION hash) and Ubuntu installs come as Posit P3M binaries — do not
  hardcode `repos=` in install steps or the binary repo is bypassed.
  This is where pull requests (including fork PRs) run — free on a public repo
  and safely isolated from the self-hosted runner.
  - `pull_request` (master/develop): full Ubuntu + macOS matrix
  - `push` (master/develop): **macOS only** — Linux pushes build on the
    self-hosted runner instead (see below); macOS can't move there
  - Runs `R CMD check` + testthat; must pass with 0 errors, 0 warnings
- **`ci-selfhosted.yml`**: Linux build/test on the `maroonbells` self-hosted
  runner. Triggers: `push` (master/develop), `workflow_dispatch`, `release` —
  **never `pull_request`** (hard security rule: fork PR code must not execute
  on the self-hosted runner). No hosted Linux fallback on push, so an offline
  runner leaves the Linux check pending rather than red.
- **`coverage.yml`**: Generates code coverage reports using `covr` package
- **`claude-review.yml`**: Automated code review assistant

### CI/CD Best Practices
- **Run checks locally first** - The pre-commit checklist above mirrors what CI/CD runs
- **Fix warnings immediately** - Don't push code that fails R CMD check
- **CI/CD is verification, not validation** - Validate locally before pushing
- **Typical CI/CD run time**: 5-15 minutes (much slower than local checks)

### Technical Details
- **RInside handling**: CI attempts to install RInside from GitHub, falls back to CRAN, or skips C++ standalone build if unavailable
- **R package focus**: CI prioritizes R package build/test which uses `NORINSIDE` flag and is more robust
- Replaces legacy Travis CI configuration

## Common Issues

### Build and Compilation
- **gfortran not found**: Install gfortran for your platform
  - macOS with Homebrew: `brew install gfortran` then `sudo ln -sf /opt/homebrew/bin/gfortran /usr/local/bin/gfortran`
  - For macOS general guidance: https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/
- **RInside dependency**: Use `-DNORINSIDE` flag to skip RInside dependency
- **Windows symbolic links**: Ensure links are properly created during clone with `git clone -c core.symlinks=true`
- **ARM64 Mac compatibility**: Using Catch2 v2.13.10 with native Apple Silicon support
- **Ubuntu/glibc compatibility**: Modern Catch2 resolves SIGSTKSZ compilation errors
- **RcppArmadillo warnings**: `-DARMA_USE_CURRENT` flag suppresses C++11 deprecation warnings

### R CMD check Warnings/Errors

**Undocumented parameters**:
```
checking Rd files ... WARNING
Undocumented arguments in documentation object 'my_function'
  'param1' 'param2'
```
Fix: Add `@param param1 Description` to Roxygen documentation

**Non-ASCII characters**:
```
checking R files for non-ASCII characters ... WARNING
Found the following files with non-ASCII characters:
  R/my_file.R
```
Fix: Replace Greek letters with Unicode escapes (μ → `\u03BC`, α → `\u03B1`)
Find them: `tools::showNonASCIIfile("R/my_file.R")`

**Empty Roxygen sections**:
```
Warning: Empty @section tag
```
Fix: Remove empty `@section` tags or fill them with content

**Missing imports**:
```
checking R code for possible problems ... NOTE
my_function: no visible global function definition for 'acf'
```
Fix: Add `@importFrom stats acf` to function documentation

**NSE variable bindings** (ggplot2/dplyr):
```
checking R code for possible problems ... NOTE
Undefined global functions or variables:
  iteration value
```
Fix: Add `utils::globalVariables(c("iteration", "value"))` at top of file

**Non-standard directories**:
```
Non-standard file/directory found at top level:
  'benchmarks'
```
Fix: Add `^benchmarks$` to `.Rbuildignore`

### Testing Framework
- Upgraded from Catch v2.0.1 (2017) to Catch2 v2.13.10 (2022) for improved reliability and platform support
