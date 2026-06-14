library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)

# Test-suite resource policy.  Each test process does within-solve threading
# (rxode2 OpenMP) plus BLAS; if (workers x within-solve threads) saturate every
# CPU, the GitHub runner's heartbeat agent is starved and the job is killed ("the
# runner has received a shutdown signal" -> exit 143).  Keep within-solve threads
# at 2 but run a SINGLE testthat worker, so the suite uses ~2 cores and leaves the
# rest for the runner agent.  (BLAS is capped to 1 thread/process in the
# environment -- see env: in .github/workflows/R-CMD-check.yaml; OpenBLAS reads
# that only at process startup, so it must be set there, not via Sys.setenv.)
# Bump the worker count with NLMIXR2_TESTTHAT_CPUS on a machine with spare cores.
.threads <- 2L                                            # within-solve threads
.ov      <- suppressWarnings(as.integer(Sys.getenv("NLMIXR2_TESTTHAT_CPUS", "")))
.cores   <- if (!is.na(.ov) && .ov >= 1L) .ov else 1L    # testthat workers (default 1)
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
