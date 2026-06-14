library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)

# Parallel-test resource policy.  testthat runs N worker processes, each doing
# within-solve threading; if  N_workers x threads_per_worker  exceeds the CPUs,
# the runner agent is starved of CPU and the job is killed ("the hosted runner
# lost communication ... starves it for CPU/Memory" -> exit 143 / 6h-cancel).
# Keep the product ~= nproc:
#   * within-solve threads     = 2 (rxode2 OpenMP inner loop)
#   * testthat workers (Ncpus) = 1 on CRAN, else floor(nproc / 2)
.onCran  <- !identical(Sys.getenv("NOT_CRAN"), "true")
.threads <- 2L
.cores   <- if (.onCran) {
  1L
} else {
  # An explicit override wins over detection.  parallel::detectCores() reports
  # the *host* CPU count -- correct on GitHub's full-VM runners, but it
  # over-reports inside CPU-limited containers/cgroups that cannot see the
  # cpuset, so containerized reproductions of CI (and constrained runners) can
  # pin the real core allotment instead of the host total.
  .ov <- suppressWarnings(as.integer(Sys.getenv("NLMIXR2_TESTTHAT_CPUS", "")))
  if (!is.na(.ov) && .ov >= 1L) {
    .ov
  } else {
    .n <- suppressWarnings(parallel::detectCores())
    if (is.na(.n) || .n < 1L) 1L else max(1L, as.integer(.n) %/% 2L)
  }
}
# testthat's worker count = getOption("Ncpus") %||% Sys.getenv("TESTTHAT_CPUS") %||% 2.
options(Ncpus = .cores)
Sys.setenv(TESTTHAT_CPUS = .cores)

# Thread caps.  setRxThreads()/setDTthreads() are runtime APIs and take effect
# immediately.  The BLAS/OpenMP *environment* variables below are read by their
# libraries only at process startup (OpenBLAS is loaded via libRblas.so when R
# itself starts), so Sys.setenv() here is TOO LATE for this (main) process -- it
# only reaches the parallel workers, which are spawned after this point and
# inherit the environment.  The main process's BLAS must therefore be capped in
# the environment BEFORE R launches; CI does this in env: (see R-CMD-check.yaml).
# Leaving OpenBLAS uncapped makes it spin one thread per core, pinning the runner
# and triggering exit-143 ("the runner has received a shutdown signal").
setRxThreads(.threads)
setDTthreads(.threads)
Sys.setenv(OMP_NUM_THREADS        = as.character(.threads))  # rxode2 OpenMP inner loop
Sys.setenv(OMP_WAIT_POLICY        = "passive")               # idle OpenMP threads sleep, not spin
Sys.setenv(GOMP_SPINCOUNT         = "0")                     # libgomp: no busy-wait before sleeping
Sys.setenv(MKL_NUM_THREADS        = "1")                     # BLAS pools (inherited by workers)
Sys.setenv(OPENBLAS_NUM_THREADS   = "1")                     # threaded OpenBLAS (inherited by workers)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")                     # Accelerate/vecLib (inherited by workers)
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  rxode2::rxUnloadAll(set = FALSE)
}

# Do NOT pass a non-parallel reporter (e.g. LocationReporter): testthat silently
# falls back to SERIAL when reporter$capabilities$parallel_support is FALSE, which
# defeats Config/testthat/parallel and moves all work into the BLAS-uncapped main
# process.  The default check reporter supports parallel.
test_check("nlmixr2est", stop_on_failure = FALSE)
