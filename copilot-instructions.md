# Copilot instructions for nlmixr2est

Purpose - Provide concise guidance for future Copilot sessions in this
repository (R package with C++ sources).

Build, test, and lint commands - Load package (compiles C/C++):
devtools::load_all() - Run full test suite: devtools::test() - Run a
single test file (example):
testthat::test_file(“tests/testthat/test-focei-1.R”) - Build
documentation (regenerates NAMESPACE): devtools::document() - Full R CMD
check: devtools::check() or from shell: R CMD check . - Note: package
has NeedsCompilation: yes and requires C++17 (see src/Makevars).
devtools::load_all() is the usual dev workflow to recompile native code.

High-level architecture (big picture) - Entry points: R/nlmixr2.R -\>
nlmixr2Est() (S3 dispatch). See R/nlmixr2Est.R for the dispatch table. -
Estimation methods: each method has an R file and control function
(e.g., R/saem.R + R/saemControl.R, R/focei.R + R/foceiControl.R). - UI →
rxode2 pipeline: model UI objects are translated to rxode2 via
R/rxPipeline.R and method-specific *RxUiGet handlers (R/saemRxUiGet.R,
R/splitModelRxUiGet.R, R/nlmeRxUiGet.R). - Pre-processing hooks:
R/hook.R registers pre-processing steps. Implementations live in
R/preProcess*.R (bounded transforms, covariates, data/ui validation,
zero-omega handling, literal fix processing). - Core numerical layer:
C++ in src/ (saem.cpp, nlm.cpp, inner.cpp, cwres.cpp, nearPD.cpp, etc.).
These expose R interfaces via Rcpp exports (RcppExports.cpp /
R/RcppExports.R). - Post-fit and accessors: Fits are nlmixr2FitData
objects; accessors and lazy computations are implemented in
R/nmObjGet\*.R and addCwres()/addNpde(). - Tests and fixtures:
helper-zzz-fits.R pre-fits models and stores .rds fixtures under
tests/testthat/fixtures to avoid re-running long estimations.

Key conventions and repository-specific patterns - “mu2”
(mu-referencing): Special handling and detection in R/mu2.R. SAEM has
optimizations that depend on mu-referencing; validate via mu2 code and
tests. - Preprocess hooks: Use R/hook.R to add/remove pre-processing
logic. Search for preProcess\* files to see available hooks. - Control
objects: sharedControl() centralizes options common to estimators;
method-specific control functions delegate to it. - Native code
workflow: C++ changes require recompilation via devtools::load_all()
(C++17, see src/Makevars). Pay attention to Rcpp and Armadillo/Eigen
usage. - Tests: Many tests rely on cached fixtures
(tests/testthat/fixtures). When adding tests that require long fits, add
fixtures to avoid CI slowness. - Naming: Estimation methods map 1:1 to
R/\*.R files (saem.R, focei.R, fo.R, nlm.R, etc.). Look there to
understand behavior and control options.

Files to read first (high signal) - CLAUDE.md (this repo) — high-level
notes about architecture and commands - README.md / README.Rmd — package
overview - R/nlmixr2Est.R — S3 dispatch and method registry -
R/rxPipeline.R and R/saemRxUiGet.R — UI → rxode2 translation - R/hook.R
and R/preProcess*.R — pre-processing hook system - src/saem.cpp,
src/nlm.cpp, src/inner.cpp — core numeric kernels -
tests/testthat/helper-zzz-fits.R and tests/testthat/fixtures/* — test
fixtures and expectations

Other AI assistant configs - CLAUDE.md present — include its guidance
when reasoning about architecture and commands. -
.claude/settings.local.json exists (local CLAUDE settings); no other
assistant rule files detected.

Quick tips for Copilot sessions - Prefer reading: R/nlmixr2Est.R,
R/rxPipeline.R, R/hook.R before proposing estimator changes. - For any
change touching src/, ensure compilation via devtools::load_all() and
run the smallest test that covers the change. - Avoid re-running long
SAEM fits in CI; use tests/fixtures and helper-zzz-fits.R when
introducing heavy tests.

Contact - If uncertain about a design decision, point to
R/sharedControl.R and the method-specific control files for precedent.
