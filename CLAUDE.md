# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package does

`nlmixr2est` provides the core nonlinear mixed-effects (NLME) estimation routines for the `nlmixr2` ecosystem. It implements multiple estimation methods for PK/PD modeling: SAEM, FOCEI, FOCE, FO/FOI, Laplace, AGQ, and NLM-family methods (nlminb, nlm, bobyqa, newuoa, n1qn1, lbfgsb3c, optim, uobyqa, nls). It depends heavily on `rxode2` for ODE solving — models are compiled to rxode2 C code at runtime.

## Common commands

```r
# Load package for interactive development
devtools::load_all()

# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-focei-1.R")

# Build documentation (regenerates NAMESPACE via roxygen2)
devtools::document()

# Full R CMD check
devtools::check()
```

The package requires compilation (`NeedsCompilation: yes`). C++17 is required (set in `src/Makevars`). After changing C++ files, `devtools::load_all()` recompiles.

## Architecture

### Estimation dispatch

Entry point is `nlmixr2()` in `R/nlmixr2.R`, which calls `nlmixr2Est()` — an S3 generic dispatching on the `est` argument string. Each estimation method has its own R file (`R/saem.R`, `R/focei.R`, `R/fo.R`, `R/foce.R`, `R/laplace.R`, `R/agq.R`, `R/nlm.R`, `R/nlminb.R`, `R/bobyqa.R`, etc.) plus a corresponding control function (e.g., `saemControl()`, `foceiControl()`).

### UI → rxode2 pipeline

Before estimation, model UI objects are translated to rxode2-compatible representations via `rxPipeline()` (`R/rxPipeline.R`) and method-specific `*RxUiGet` handlers (`R/saemRxUiGet.R`, `R/saemRxUiGetModel.R`, `R/splitModelRxUiGet.R`, `R/nlmeRxUiGet.R`). The `ini` block sets initial estimates; mu-referencing (where population parameters appear linearly in ETAs) is detected and handled specially for SAEM efficiency.

### Pre-processing hooks

`R/hook.R` provides a hook system for pre-processing steps that run before estimation: bounded parameter transforms (`R/preProcessBoundedTransform.R`), covariate handling (`R/preProcessCovariatesPresent.R`), data/UI validation (`R/preProcessDataUi.R`), zero-omega fixing (`R/preProcessZeroOmega.R`), and literal fix processing (`R/preProcessLiteralFix.R`).

### C++ layer (`src/`)

Core numerical work is in C++17 with RcppArmadillo/RcppEigen:
- `saem.cpp` — SAEM algorithm kernel
- `nlm.cpp` / `innerPoint.cpp` — NLM-family optimization
- `cwres.cpp`, `npde.cpp`, `ires.cpp` — residual calculations
- `censResid.cpp`, `censEst.cpp` — M2/M3/M4 censoring methods
- `nearPD.cpp` — nearest positive-definite matrix
- `uninformativeEtas.cpp`, `shrink.cpp` — ETA diagnostics

### Post-fit objects

Fit results are `nlmixr2FitData` objects (a data frame subclass). Post-fit accessors use `nmObjGet` S3 dispatch (`R/nmObjGet.R`, `R/nmObjHandle.R`). Residuals are added lazily via `addCwres()` and `addNpde()`.

### Test fixtures

`tests/testthat/helper-zzz-fits.R` pre-fits models and caches them so individual test files can reference fitted objects without re-running estimation. Fixture `.rds` files live in `tests/testthat/fixtures/`.

## Key conventions

- `mu2` referencing (`R/mu2.R`) detects and validates mu-referenced parameters; violations produce warnings, not errors, for SAEM
- IOV (inter-occasion variability) support is in `R/iov.R` with special handling in the mu-referencing hooks
- The `sharedControl()` function (`R/sharedControl.R`) defines options common across all estimation methods
- `R/nlmixr2Est.R` contains the central S3 dispatch table — look here to understand what methods exist
