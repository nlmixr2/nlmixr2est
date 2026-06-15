library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)

# Test thread policy (kept deliberately simple):
#   * testthat workers: a SINGLE worker on CI or CRAN, so the suite does not
#     oversubscribe a shared / core-limited runner; everywhere else testthat
#     manages Config/testthat/parallel normally.
#   * rxode2 (and data.table) within-solve threads: capped to 2 on CRAN, per
#     CRAN's two-core policy; on CI and locally rxode2 manages its own threads.
.onCran <- !identical(Sys.getenv("NOT_CRAN"), "true")
.onCI   <- isTRUE(as.logical(Sys.getenv("CI", "false")))

if (.onCI || .onCran) {
  options(Ncpus = 1L)
  Sys.setenv(TESTTHAT_CPUS = "1")
}
if (.onCran) {
  setRxThreads(2L)
  setDTthreads(2L)
}
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  rxode2::rxUnloadAll(set = FALSE)
}

test_check("nlmixr2est", stop_on_failure = FALSE)
