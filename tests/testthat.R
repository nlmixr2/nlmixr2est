library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)

# Test thread policy (kept deliberately simple):
#   * testthat workers: a SINGLE worker on CI or CRAN, so the suite does not
#     oversubscribe a shared / core-limited runner; everywhere else testthat
#     manages Config/testthat/parallel normally.
#   * rxode2 (and data.table) within-solve threads: capped to 2 on CI and CRAN.
#     An uncapped thread pool oversubscribes the core-limited hosted runners --
#     it multiplies both CPU contention and the per-thread solve-buffer memory,
#     which produced exit-143 "runner has received a shutdown signal" kills
#     mid-suite (seen consistently on the oldrel-1 leg).  2 matches CRAN's
#     two-core policy.  Locally rxode2 manages its own threads.
.onCran <- !identical(Sys.getenv("NOT_CRAN"), "true")
.onCI   <- isTRUE(as.logical(Sys.getenv("CI", "false")))

if (.onCI || .onCran) {
  options(Ncpus = 1L)
  Sys.setenv(TESTTHAT_CPUS = "1")
  setRxThreads(2L)
  setDTthreads(2L)
}
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  rxode2::rxUnloadAll(set = FALSE)
}

# -------------------------------------------------------------------------
# CI test partitioning
#
# The full suite takes many hours on a single-worker CI runner, so push/PR
# R-CMD-check runs only the "essential" subset -- every test file EXCEPT the
# slow ones listed in .slowBatches below (still much more than CRAN).
#
# The slow files run separately in the slow-tests workflow, split into
# batches that each stay under ~1h and run one-at-a-time (non-overlapping).
# That workflow sets NLMIXR2EST_TEST_BATCH=<n> to run only batch n's files.
#
# Names are the test file basename with the leading "test-" and trailing
# ".R" removed (what testthat's `filter` matches).  When a test file grows
# past a few minutes, move it into one of the batches here.
# -------------------------------------------------------------------------
# Batches are sized from measured single-worker run times so each stays well
# under an hour; the slow-tests workflow runs them one at a time.
.slowBatches <- list(
  # batch 1
  c("focei-wang2007-boxcox", "focei-wang2007-combined", "vpcSim"),
  # batch 2
  c("focei-wang2007-lognormal", "cov-analytic", "focei-wang2007-power"),
  # batch 3
  c("focei-wang2007-boxcox-half", "nlm-cens", "issue-429",
    "focei-wang2007-bounded"),
  # batch 4
  c("impmap", "matexp", "mufocei", "focei-wang2007-yeojohnson",
    "focei-wang2007-boxcox-lnorm", "nlme", "focei-fast-grad"),
  # batch 5
  c("focei-llik", "iov", "nlm-adjoint", "saem-mix", "posthoc", "ar-est",
    "mu-family", "vae-fit", "focei-wang2007-basic", "vae-neonatal",
    "vae-errmodel", "table-cmt", "vae-covariate")
)
.slowAll <- unlist(.slowBatches)

.batch <- Sys.getenv("NLMIXR2EST_TEST_BATCH")

.filter <- NULL
if (nzchar(.batch)) {
  # Slow-batch mode: run ONLY this batch's slow files.
  .b <- suppressWarnings(as.integer(.batch))
  if (is.na(.b) || .b < 1L || .b > length(.slowBatches)) {
    stop(sprintf("NLMIXR2EST_TEST_BATCH=%s out of range (1..%d)",
                 .batch, length(.slowBatches)))
  }
  .files <- .slowBatches[[.b]]
  if (length(.files) == 0L) {
    .filter <- "^$"  # empty batch: match nothing
  } else {
    .filter <- paste0("^(", paste(.files, collapse = "|"), ")$")
  }
} else if (.onCI && !.onCran && length(.slowAll) > 0L) {
  # Essential subset on push/PR CI: everything EXCEPT the slow files.
  .filter <- paste0("^(?!(", paste(.slowAll, collapse = "|"), ")$)")
}
# Locally (and on CRAN) .filter stays NULL -> run everything.

# Progress side-channel (CI only): flush the name of each test file/test to a
# file at a stable path ($RUNNER_TEMP) as it starts.  R CMD check buffers the
# test .Rout and discards it when the runner is killed mid-suite (exit 143
# "the runner has received a shutdown signal"), so that buffer never shows the
# culprit.  An `always()` workflow step prints this file, whose last line is the
# exact test that was running when the runner died.
# Default reporter is the check reporter ("check"); on CI wrap it in a
# MultiReporter alongside a tiny progress reporter that flushes to $RUNNER_TEMP.
# If testthat's (unexported) CheckReporter is unavailable, fall back to the
# plain default so the suite still runs.
.reporter <- testthat::check_reporter()
if (.onCI) {
  .progLog <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()),
                        "testthat-progress.log")
  try(cat(sprintf("[%s] progress log start\n",
                  format(Sys.time(), "%H:%M:%OS2")),
          file = .progLog), silent = TRUE)
  .ProgReporter <- R6::R6Class(
    "ProgReporter", inherit = testthat::Reporter,
    public = list(
      logmsg = function(msg) {
        try(cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%OS2"), msg),
                file = .progLog, append = TRUE), silent = TRUE)
      },
      start_file = function(filename) self$logmsg(paste("FILE", filename)),
      start_test = function(context, test) self$logmsg(paste("   test:", test))
    )
  )
  .checkRep <- tryCatch(
    getFromNamespace("CheckReporter", "testthat")$new(),
    error = function(e) NULL)
  if (!is.null(.checkRep)) {
    .reporter <- testthat::MultiReporter$new(
      reporters = list(.checkRep, .ProgReporter$new()))
  }
}

test_check("nlmixr2est", stop_on_failure = FALSE, filter = .filter,
           perl = TRUE, reporter = .reporter)
