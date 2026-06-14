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

# Within-solve thread caps.  setRxThreads()/setDTthreads() are runtime APIs, so
# they take effect immediately.  BLAS/OpenMP *environment* caps are deliberately
# NOT set here: those libraries read their env vars only at process startup
# (OpenBLAS is loaded via libRblas.so when R itself starts), so Sys.setenv()
# would be too late to have any effect.  They are set in the environment before R
# launches instead -- see env: in .github/workflows/R-CMD-check.yaml.  Without
# that cap OpenBLAS spins one thread per core, which together with the solve
# threads saturates the CI runner and the job is killed with exit 143.
setRxThreads(.threads)
setDTthreads(.threads)
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  rxode2::rxUnloadAll(set = FALSE)
}

# Do NOT pass a non-parallel reporter (e.g. LocationReporter): testthat silently
# falls back to SERIAL when reporter$capabilities$parallel_support is FALSE, which
# defeats Config/testthat/parallel and moves all work into the BLAS-uncapped main
# process.  The default check reporter supports parallel.
test_check("nlmixr2est", stop_on_failure = FALSE)
