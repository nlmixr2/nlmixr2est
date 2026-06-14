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
  .n <- suppressWarnings(parallel::detectCores())
  if (is.na(.n) || .n < 1L) 1L else max(1L, as.integer(.n) %/% 2L)
}
# testthat's worker count = getOption("Ncpus") %||% Sys.getenv("TESTTHAT_CPUS") %||% 2.
options(Ncpus = .cores)
Sys.setenv(TESTTHAT_CPUS = .cores)

# Apply the caps BEFORE test_check() spawns the workers, so they inherit them.
# Within-solve threads = 2; BLAS pools single-threaded per worker (parallelism is
# workers x the 2 solve threads); idle OpenMP threads sleep instead of spinning
# (spinning pins every core and starves the runner agent).
setRxThreads(.threads)
setDTthreads(.threads)
Sys.setenv(OMP_NUM_THREADS        = as.character(.threads))  # rxode2 OpenMP inner loop
Sys.setenv(OMP_WAIT_POLICY        = "passive")               # idle OpenMP threads sleep, not spin
Sys.setenv(GOMP_SPINCOUNT         = "0")                     # libgomp: no busy-wait before sleeping
Sys.setenv(MKL_NUM_THREADS        = "1")                     # BLAS pools: 1 thread per worker
Sys.setenv(OPENBLAS_NUM_THREADS   = "1")                     # threaded OpenBLAS (Linux CI)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")                     # Accelerate/vecLib (macOS)
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  rxode2::rxUnloadAll(set = FALSE)
}

test_check("nlmixr2est", stop_on_failure = FALSE,
           reporter = testthat::LocationReporter)
