library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)

# Test thread policy (kept deliberately simple):
#   * testthat workers: SERIAL (in-process) on CI or CRAN; everywhere else
#     testthat manages Config/testthat/parallel normally.  Serial rather than
#     a single parallel worker: with 1 worker parallel mode adds nothing but
#     the callr message pipe, and that pipe base64-serializes every non-success
#     test event -- an *erroring* test whose backtrace inlines a large object
#     (any do.call(f, list(<fit/data>)) frame) ships the whole object through
#     it at ~18x its size in the orchestrator (measured), and a single >1.5GB
#     payload kills the whole run with "Error in gsub: result string is too
#     long" (R's 2GB string cap).  This is what "the runner has received a
#     shutdown signal" / runner-OOM windows+devel failures were.  Serial mode
#     has no pipe, and streams output so a hung/failing file is identifiable
#     from testthat.Rout.
#   * rxode2 (and data.table) within-solve threads: capped to 2 on CRAN, per
#     CRAN's two-core policy; on CI and locally rxode2 manages its own threads.
.onCran <- !identical(Sys.getenv("NOT_CRAN"), "true")
.onCI <- isTRUE(as.logical(Sys.getenv("CI", "false")))

if (.onCI || .onCran) {
  options(Ncpus = 1L)
  Sys.setenv(TESTTHAT_CPUS = "1")
  Sys.setenv(TESTTHAT_PARALLEL = "FALSE")
}
if (.onCran) {
  setRxThreads(2L)
  setDTthreads(2L)
}
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  rxode2::rxUnloadAll(set = FALSE)
}

# -------------------------------------------------------------------------
# CI test selection lives in .github/workflows/ci.yaml, not here: the fast tier
# runs tests/ci-core.txt and the nightly full tier runs every file, both
# sharded by .github/scripts/ci-shards.R and run through testthat::test_dir().
# R CMD check (NOT_CRAN=false) runs the CRAN-visible subset below.  Coverage
# sets NLMIXR2EST_TEST_FILTER to the files it can afford.
# -------------------------------------------------------------------------
.filter <- NULL
if (nzchar(Sys.getenv("NLMIXR2EST_TEST_FILTER"))) {
  .filter <- Sys.getenv("NLMIXR2EST_TEST_FILTER")
} else if (.onCran) {
  # CRAN-visible files that compile models; they pushed CRAN's Windows check
  # past its 10 min limit.  CI's NOT_CRAN test shards still run them.
  .cranSkip <- c(
    "focei-lincmt-carry-trans",
    "focei-lincmt-carry-ll",
    "nlmixr2output",
    "focei-lincmt-carry-jump",
    "focei-lincmt-carry-emit-fd",
    "imp-cov",
    "focei-lincmt-carry-eligible-gates",
    "iov-method-gate",
    "focei-lincmt-carry-emit",
    "rxuideparse",
    "cov-analytic",
    "lint-superassign",
    "saem-iov-twolevel",
    "focei-lincmt-carry-eligible",
    "focei-family-control",
    "cov",
    "focei-outer-trust"
  )
  .filter <- paste0("^(?!(", paste(.cranSkip, collapse = "|"), ")$)")
}
# Locally .filter stays NULL -> run everything.

# stop_on_failure MUST stay TRUE (the default).  It was set to FALSE in
# 4af6a7eba as a debugging aid ("Be more verbose in testthat to track-down
# future seg faults"), and with it R CMD check never errors on a test failure
# -- so `checking tests` passed no matter what.  The run on main at 722e5775
# reported SUCCESS while its three shards printed FAIL 8, FAIL 2 and FAIL 3;
# thirteen real failures had reached main with a green badge the whole way.
# If a segfault ever needs that verbosity again, set it for the one run
# (R_TESTTHAT_... / an env guard), never in the committed harness.
test_check("nlmixr2est", filter = .filter, perl = TRUE)
