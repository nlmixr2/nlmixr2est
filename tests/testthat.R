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
    "qrpem-slow"),
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
    "vae-errmodel", "table-cmt", "vae-covariate"),
  # batch 6 -- heaviest remaining files on the single-worker CI runner
  # (VAE internals + a few slow structural tests), moved out of the essential
  # push/PR subset to trim its wall time / reclamation exposure.
  c("vae-encoder", "vae-train", "vae-decoder", "vae-elbo", "vae-inner",
    "vae-fixbounds", "vae-parhist", "vae-iov", "split", "unary-mu", "timing"),
  # batches 7-11 -- est="rpem".  The quick core-functionality files (the C++ E-step /
  # M-step units, the LL-model build, the dispatch/control smoke, and a basic end-to-end
  # fit: rpem-cpp-estep, rpem-cpp-mstep, rpem-llik-model, rpem-est, rpem-fit) stay in the
  # essential push/PR subset; every other (multi-iteration, fit-based) rpem file runs
  # weekly here.  Sized from measured single-worker times to balance the batches.
  # batch 7
  c("rpem-mix-pow", "rpem-cens", "rpem-pow", "rpem-comb"),
  # batch 8
  c("rpem-multi-pow", "rpem-covariate", "rpem-mix-multiparam", "rpem-impinflate",
    "rpem-tbs", "rpem-parallel"),
  # batch 9
  c("rpem-fisher", "rpem-multi-comb", "rpem-mix-tbs", "rpem-llik", "rpem-boundfix",
    "rpem-lnorm"),
  # batch 10
  c("rpem-multi", "rpem-mix", "rpem-mix-percomp", "rpem-iov", "rpem-multi-lnorm",
    "rpem-parhist", "rpem-focei-agreement"),
  # batch 11
  c("rpem-mix-guards", "rpem-multi-tbs", "rpem-mix-comb", "rpem-jump", "rpem-struct",
    "rpem-prop"),
  # batch 12 -- advi (variational inference) multi-iteration fits
  c("advi-repro", "advi-focei-agreement", "advi-neonatal", "advi-fullrank",
    "advi-fullbayes")
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
