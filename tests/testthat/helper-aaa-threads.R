# Sourced first (sorts before the other helper-*.R) in the main test process
# AND in every parallel testthat worker, before any fixture is built or test
# runs.  Two jobs:
#
# 1. Keep within-solve threading small so (worker count) x (threads-per-solve)
#    does not oversubscribe the CPUs.  With Config/testthat/parallel the
#    parallelism comes from the workers, so each worker solves with few threads.
#
# 2. Give each process its OWN rxode2 model-compile directory.  rxode2 resolves
#    its compile/cache dir once and caches it in the internal `.rxTempDir0`;
#    by default that is a single shared path (R_user_dir("rxode2", "cache")).
#    Parallel workers all compiling models into that one directory race and
#    corrupt each other's intermediate files, which surfaces as
#    "error building model" / "something went wrong in compilation".  Point
#    rxTempDir at a per-process directory and clear the cached value so the new
#    location takes effect on the next (i.e. first real) model compile.
local({
  try(rxode2::setRxThreads(2L), silent = TRUE)
  try(data.table::setDTthreads(2L), silent = TRUE)
  .dir <- file.path(tempdir(), paste0("rxcompile-", Sys.getpid()))
  dir.create(.dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(rxTempDir = .dir)
  # rxode2 caches the resolved temp dir during .onLoad; reset it so the
  # per-process directory above is picked up.  Guarded: if a future rxode2
  # changes this internal, we simply fall back to the default behaviour.
  try(utils::assignInNamespace(".rxTempDir0", NULL, ns = "rxode2"), silent = TRUE)
  invisible(try(rxode2::rxTempDir(), silent = TRUE))
})
