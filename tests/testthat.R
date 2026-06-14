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

# Hard-cap every thread pool BEFORE test_check() spawns the parallel workers, so
# the workers inherit these in their environment.  Without this each worker
# multiplies the host's cores: the CI runners link the threaded BLAS
# (libopenblas-pthread), which by default busy-waits on one thread PER CORE for
# every matrix op, so N workers x (cores) BLAS threads oversubscribe and hang the
# runner -- it goes unresponsive and the job is killed (exit 143).  This is the
# real cause of the ubuntu-devel/oldrel failures; the faster ubuntu-release and
# reference-BLAS dev machines squeak through.  Parallelism comes from the workers,
# so BLAS is single-threaded per worker.
setRxThreads(2L)
setDTthreads(2L)
Sys.setenv(OMP_NUM_THREADS = "2")          # OpenMP (rxode2 inner loop)
Sys.setenv(MKL_NUM_THREADS = "2")          # Intel MKL, if linked
Sys.setenv(OPENBLAS_NUM_THREADS = "1")     # threaded OpenBLAS (Linux CI runners)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")   # Accelerate/vecLib (macOS)
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  rxode2::rxUnloadAll(set = FALSE)
}

test_check("nlmixr2est", stop_on_failure = FALSE,
           reporter = testthat::LocationReporter)
