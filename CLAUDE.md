# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this package does

`nlmixr2est` provides the core nonlinear mixed-effects (NLME) estimation
routines for the `nlmixr2` ecosystem. It implements multiple estimation
methods for PK/PD modeling: SAEM, FOCEI, FOCE, FO/FOI, Laplace, AGQ, and
NLM-family methods (nlminb, nlm, bobyqa, newuoa, n1qn1, lbfgsb3c, optim,
uobyqa, nls). It depends heavily on `rxode2` for ODE solving – models
are compiled to rxode2 C code at runtime.

## Common commands

``` r

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

# Set the number of threads for parallel processing (e.g., SAEM)
rxode2::setRxThreads(threads = 4)  # Use 4 threads
```

The package requires compilation (`NeedsCompilation: yes`). C++17 is
required (set in `src/Makevars`). After changing C++ files,
[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
recompiles.

## Architecture

### Estimation dispatch

Entry point is
[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md)
in `R/nlmixr2.R`, which calls
[`nlmixr2Est()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Est.md)
– an S3 generic dispatching on the `est` argument string. Each
estimation method has its own R file (`R/saem.R`, `R/focei.R`, `R/fo.R`,
`R/foce.R`, `R/laplace.R`, `R/agq.R`, `R/nlm.R`, `R/nlminb.R`,
`R/bobyqa.R`, etc.) plus a corresponding control function (e.g.,
[`saemControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md),
[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md)).

### UI -\> rxode2 pipeline

Before estimation, model UI objects are translated to rxode2-compatible
representations via `rxPipeline()` (`R/rxPipeline.R`) and
method-specific `*RxUiGet` handlers (`R/saemRxUiGet.R`,
`R/saemRxUiGetModel.R`, `R/splitModelRxUiGet.R`, `R/nlmeRxUiGet.R`). The
`ini` block sets initial estimates; mu-referencing (where population
parameters appear linearly in ETAs) is detected and handled specially
for SAEM efficiency.

### Pre-processing hooks

`R/hook.R` provides a hook system for pre-processing steps that run
before estimation: bounded parameter transforms
(`R/preProcessBoundedTransform.R`), covariate handling
(`R/preProcessCovariatesPresent.R`), data/UI validation
(`R/preProcessDataUi.R`), zero-omega fixing (`R/preProcessZeroOmega.R`),
and literal fix processing (`R/preProcessLiteralFix.R`).

### C++ layer (`src/`)

Core numerical work is in C++17 with RcppArmadillo/RcppEigen: -
`saem.cpp` – SAEM algorithm kernel - `nlm.cpp` / `innerPoint.cpp` –
NLM-family optimization - `cwres.cpp`, `npde.cpp`, `ires.cpp` – residual
calculations - `censResid.cpp`, `censEst.cpp` – M2/M3/M4 censoring
methods - `nearPD.cpp` – nearest positive-definite matrix -
`uninformativeEtas.cpp`, `shrink.cpp` – ETA diagnostics

### Registering C/C++ (`.Call`) routines – `src/init.c` is MANUAL

`src/init.c` hand-maintains the `R_registerRoutines` `callMethods[]`
table (Rcpp’s generated `R_init_nlmixr2est` in `RcppExports.cpp` is NOT
used – init.c owns `R_init_nlmixr2est`). Each `.Call` entry lists an
explicit prototype and a HARDCODED argument count,
e.g. `{"_nlmixr2est_foo", (DL_FUNC) &_nlmixr2est_foo, 5}`. So when you
add a `[[Rcpp::export]]` function OR change an existing one’s argument
count, you MUST do BOTH: 1. `Rcpp::compileAttributes(".")` – regenerates
`src/RcppExports.cpp` + `R/RcppExports.R`; and 2. hand-edit `src/init.c`
– add/update the forward `SEXP ...(SEXP, ...)` prototype AND the
`callMethods[]` entry’s arg count.

Skipping step 2 gives a runtime `.Call` error like
`Incorrect number of arguments (N), expecting M for '_nlmixr2est_foo'`
(the loaded registration table still has the old arity), NOT a compile
error – so it survives a clean rebuild. A header change still needs
`rm -f src/*.o src/*.so` before rebuilding.

### Post-fit objects

Fit results are `nlmixr2FitData` objects (a data frame subclass).
Post-fit accessors use `nmObjGet` S3 dispatch (`R/nmObjGet.R`,
`R/nmObjHandle.R`). Residuals are added lazily via
[`addCwres()`](https://nlmixr2.github.io/nlmixr2est/reference/addCwres.md)
and
[`addNpde()`](https://nlmixr2.github.io/nlmixr2est/reference/addNpde.md).

### Test fixtures

`tests/testthat/helper-zzz-fits.R` pre-fits models and caches them so
individual test files can reference fitted objects without re-running
estimation. Fixture `.rds` files live in `tests/testthat/fixtures/`.

### Test thread policy

`tests/testthat.R` keeps CI and CRAN from oversubscribing core-limited
runners (which historically caused 6h timeouts / exit-143 “the runner
has received a shutdown signal”):

- **testthat workers**: a single worker on CI (`CI=true`) or CRAN;
  everywhere else testthat manages `Config/testthat/parallel`
  (`parallel: true`, `edition: 3` in `DESCRIPTION`) normally.
- **rxode2 / data.table within-solve threads**: capped to 2 on CRAN
  only; on CI and locally rxode2 manages its own threads.
- **Do not pass a non-parallel reporter** to `test_check()` – testthat
  silently falls back to serial when
  `reporter$capabilities$parallel_support` is FALSE
  (e.g. `LocationReporter`), which defeats the parallel config.

## Key conventions

- `mu2` referencing (`R/mu2.R`) detects and validates mu-referenced
  parameters; violations produce warnings, not errors, for SAEM
- IOV (inter-occasion variability) support is in `R/iov.R` with special
  handling in the mu-referencing hooks
- The `sharedControl()` function (`R/sharedControl.R`) defines options
  common across all estimation methods
- `R/nlmixr2Est.R` contains the central S3 dispatch table – look here to
  understand what methods exist

## Comment and documentation style

- Avoid AI-typical verbosity: no multi-paragraph explanations in code
  comments, no worked examples or numeric anecdotes embedded in
  comments, no restating what the code obviously does. One line is
  usually enough; only explain non-obvious *why*.
- Roxygen descriptions should be short (a sentence or two) – keep every
  `@param`/`@return` tag, but don’t pad `@details` with tutorial-style
  prose.
- `NEWS.md` entries should be terse, single-bullet, past-tense, similar
  to the existing entries under older version headers
  (e.g. `# nlmixr2est 6.0.1` / `6.0.0`) – not multi-line narrations of
  root cause and mechanism.
- No Unicode in source, comments, or docs (em dashes, smart quotes,
  arrows, etc.) – CRAN policy expects plain ASCII. Use `--`, `'`/`"`, or
  `->` instead.

## Coding Conventions

- **R Style**: All variable and function names in R must be in
  `camelCase` (e.g., `simData`, `fitSaem`, `oneCompartmentMix`). Do not
  use `snake_case` or `dot.case` for R objects.
- **C++ Style**: Follow Rcpp/Armadillo conventions. Use `snake_case` or
  `camelCase` as appropriate.

## Generating rxode2 model text (augmented/inner/outer models)

When emitting rxode2 model text programmatically (e.g. the augmented
sensitivity models in `R/foceiCovAnalytic.R` / `R/foceiGradAnalytic.R`,
or the inner/pred models in `R/focei.R`):

- **Name every generated output variable `rx_<name>_`** (leading `rx_`,
  trailing `_`), matching the inner model’s `rx_pred_` / `rx__sens_..._`
  convention. Two reasons: (1) the user can define a model variable with
  any plain name (e.g. `rvar`, `predf`, `f1`), so an un-prefixed
  generated name can collide and silently corrupt the model; (2)
  **rxode2 parses a top-level `x = <constant>` as an initial value** (x
  becomes a parameter, dropping out of the solve columns) UNLESS the
  name is of the `rx_<name>_` form – those keep a constant assignment as
  a real lhs output column. So a constant sensitivity
  (e.g. `d(R)/d(eta)=0` for additive error) comes back as a column of
  that constant, no special-casing needed. Use `rx_predf_`,
  `rx_f1_<dir>`, `rx_rvarf_`, etc. – never `predf`, `f1_...`, `rvar`,
  `rxPredF` (a plain `rxFoo` without the trailing `_` is still parsed as
  an init).
- **Use `=` only for a variable you will read back** (an output column
  of the solve). Use `~` (suppressed assignment) for any intermediate
  you compute but do NOT output – extra `=` output columns shift the
  positional column layout and add solve overhead.
- **Declare the theta/eta inputs AND the model covariates (`ui$allCovs`)
  with `param(...)`** so the solve parameter order is fixed and
  positional (values are keyed by name, but pinning order keeps the
  input layout stable across builds). Covariates are supplied from the
  event data at solve time but must still be declared in `param(...)`.
