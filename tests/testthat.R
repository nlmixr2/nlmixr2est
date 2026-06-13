library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)

# Parallel-test resource policy.
#  * On CRAN: deterministic and small -- at most 2 cores (CRAN's limit and the
#    <10 min budget).  Most heavy fits live in nmTest({}) blocks that skip on
#    CRAN anyway.
#  * Otherwise (CI, local dev): use the whole machine.  detectCores() auto-
#    adapts -- 4 on the GitHub windows/ubuntu runners today, 32 on a dev box --
#    so we never have to hard-code the runner size.
.onCran <- !identical(Sys.getenv("NOT_CRAN"), "true")
.cores <- if (.onCran) {
  2L
} else {
  .n <- suppressWarnings(parallel::detectCores())
  if (is.na(.n) || .n < 1L) 2L else as.integer(.n)
}
# testthat's worker count = getOption("Ncpus") %||% Sys.getenv("TESTTHAT_CPUS") %||% 2.
# Set both so the worker count is what we intend regardless of the CI image.
options(Ncpus = .cores)
Sys.setenv(TESTTHAT_CPUS = .cores)

# Keep within-solve threads small in the main process (each parallel worker does
# the same via helper-aaa-threads.R, which also isolates the per-worker rxode2
# model-compile directory so concurrent compiles don't race).
setRxThreads(2L)
setDTthreads(2L)
Sys.setenv(OMP_NUM_THREADS = "2")
Sys.setenv(MKL_NUM_THREADS = "2")
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  rxode2::rxUnloadAll(set = FALSE)
}

test_check("nlmixr2est", stop_on_failure = FALSE,
           reporter = testthat::LocationReporter)
