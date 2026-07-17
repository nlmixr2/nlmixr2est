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
.onCI   <- isTRUE(as.logical(Sys.getenv("CI", "false")))

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
  c("focei-wang2007-boxcox", "focei-wang2007-combined", "vpcSim",
    "qrpem-slow", "focei-foce-plus"),
  # batch 2
  c("focei-wang2007-lognormal", "cov-analytic", "focei-wang2007-power",
    "fsaem", "cov-condition"),
  # batch 3
  c("focei-wang2007-boxcox-half", "nlm-cens", "issue-429", "issue-470",
    "focei-wang2007-bounded", "saem-loglik", "mu-timevarying", "saem-nearpd",
    "saem-nonmutheta", "saem-sharedinner", "focei-theta-reset-bounds"),

  # batch 4
  c("impmap", "matexp", "mfocei", "focei-wang2007-yeojohnson",
    "focei-wang2007-boxcox-lnorm", "nlme", "focei-fast-grad"),
  # batch 5
  c("focei-llik", "iov", "nlm-adjoint", "saem-mix", "saem-mix-regress", "posthoc", "ar-est",
    "mu-family", "mu-plain-fit", "vae-fit", "focei-wang2007-basic",
    "vae-neonatal", "vae-errmodel", "table-cmt", "vae-covariate"),
  # batch 6 -- heaviest remaining files on the single-worker CI runner
  # (VAE internals + a few slow structural tests), moved out of the essential
  # push/PR subset to trim its wall time / reclamation exposure.
  c("vae-encoder", "vae-train", "vae-decoder", "vae-elbo", "vae-inner",
    "vae-fixbounds", "vae-parhist", "vae-iov", "split", "unary-mu", "timing"),
  # batch 7 -- advi (variational inference) multi-iteration fits
  c("advi-repro", "advi-focei-agreement", "advi-neonatal", "advi-fullrank",
    "advi-fullbayes"),
  # batch 8 -- nonparametric (npag/npb) fit-based validation.  These set up the
  # FOCEi inner problem and run full NPAG cycles / independent solves, so they are
  # much slower than the essential npag unit tests (dispatch/ipm/grid, which stay
  # in the push/PR subset) and run weekly only.
  c("npag-psi", "npag-cycle", "npag-fit", "npb-fit", "npag-bimodal", "npag-fixed",
    "npag-error-models", "npag-mixture", "npag-general-lik", "npag-muexpand", "npag-golden")
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

test_check("nlmixr2est", stop_on_failure = FALSE, filter = .filter, perl = TRUE)
