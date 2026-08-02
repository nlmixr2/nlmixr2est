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
# Platforms that only need to prove the package BUILDS and checks (macOS and
# Windows in R-CMD-check.yaml) run the CRAN-visible subset instead of the full
# NOT_CRAN suite.  The same `checking tests` step costs ~19 min on an Ubuntu
# runner and ~5h45m on macOS/Windows -- model compilation dominates there -- which
# blew GitHub's hard 6-hour job limit and left those jobs cancelled, i.e. giving
# no signal at all.  Ubuntu still runs the full essential subset.
if (isTRUE(as.logical(Sys.getenv("NLMIXR2EST_TEST_CRAN_ONLY", "false")))) {
  Sys.setenv("NOT_CRAN" = "false") # nolint
}

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
    "cov-condition", "agq-cov", "cov-decouple-saimp"),
  # batch 3
  c("focei-wang2007-boxcox-half", "nlm-cens", "issue-429", "issue-470",
    "focei-wang2007-bounded", "saem-loglik", "mu-timevarying", "saem-nearpd",
    "saem-nonmutheta", "focei-theta-reset-bounds",
    "saem-cov-analytic", "focei-shi21-bounds", "splitbolus-interp"),

  # batches 4-7 -- the former batches 4 (13 files) and 5 (15 files), each split
  # in two.  Both TIMED OUT on the weekly runner (run 30688874991, 2026-08-01):
  # "the runner has received a shutdown signal", exit 143, after ~3h06m of test
  # time, so they reported NO assertion results at all.  A meaningful slice of
  # the suite therefore had no CI signal, which is how a wrong nlm adjoint
  # gradient and several stale assertions survived here unnoticed.  Each
  # focei-wang2007-* file is kept in a separate batch -- one is ~10 error models
  # x 6 fits (ODE and solved-form) and dominates whatever batch it lands in.
  # batch 4
  c("impmap", "matexp", "mfocei", "focei-wang2007-yeojohnson",
    "nlme", "focei-fast-grad", "lincmt-ode-fit"),
  # batch 5
  c("focei-wang2007-boxcox-lnorm", "nlme-cov", "agq-fast-grad",
    "focei-ll-fast-grad-fit", "focei-fast-methods-fit", "odeswap-fit"),
  # batch 6
  c("focei-llik", "iov", "iov-zero-eta", "saem-mix", "saem-mix-regress",
    "posthoc", "ar-est", "mu-family"),
  # batch 7
  c("mu-plain-fit", "vae-fit", "focei-wang2007-basic", "vae-neonatal",
    "vae-errmodel", "table-cmt", "vae-covariate-selection"),
  # batch 8 -- heaviest remaining files on the single-worker CI runner
  # (VAE internals + a few slow structural tests), moved out of the essential
  # push/PR subset to trim its wall time / reclamation exposure.
  c("vae-encoder", "vae-train", "vae-decoder", "vae-elbo", "vae-inner",
    "vae-fixbounds", "vae-parhist", "vae-iov", "vae-grad-fit", "vae-ll-grad-fit",
    "vae-l0learn-fit", "vae-hockey-fit", "split", "unary-mu", "timing", "bounded-transform"),
  # batch 9 -- emvi/fbvi (variational inference) multi-iteration fits, plus the
  # cross-method omega off-diagonal fit checks (vae/emvi/npag/npb)
  c("vi-repro", "vi-focei-agreement", "vi-neonatal", "vi-fullrank",
    "vi-fullbayes", "vi-stan", "omega-offdiag", "augpred", "vae-residopt"),
  # batch 10 -- nonparametric (npag/npb) fit-based validation.  These set up the
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
  #
  # NLMIXR2EST_TEST_SHARD="i/n" splits that subset across n parallel jobs.  A
  # hosted runner reclaims a long job mid-run ("exit code 143 / The runner has
  # received a shutdown signal"), and runner speed varies enormously -- the same
  # `checking tests` step measured 27m on one ubuntu-release runner and was still
  # going at 4h15m on another before being reclaimed.  Sharding keeps every job
  # short enough to finish without dropping a single test file.
  .shard <- Sys.getenv("NLMIXR2EST_TEST_SHARD")
  if (nzchar(.shard)) {
    .p <- strsplit(.shard, "/", fixed = TRUE)[[1]]
    .i <- suppressWarnings(as.integer(.p[1])); .n <- suppressWarnings(as.integer(.p[2]))
    if (length(.p) != 2L || is.na(.i) || is.na(.n) || .n < 1L || .i < 1L || .i > .n) {
      stop(sprintf("NLMIXR2EST_TEST_SHARD=%s must be \"i/n\" with 1 <= i <= n", .shard))
    }
    # Deal the files out round-robin over a SORTED list so the split is stable
    # across jobs and platforms, and every file lands in exactly one shard.
    .all <- sort(sub("^test-", "", sub("\\.R$", "",
                                       basename(Sys.glob(file.path("testthat", "test-*.R"))))))
    .ess <- setdiff(.all, .slowAll)
    # An empty list would make every shard match nothing and report a green run
    # having tested nothing at all -- fail loudly instead.
    if (length(.ess) == 0L) {
      stop("NLMIXR2EST_TEST_SHARD set but no test files found in ./testthat")
    }
    .mine <- .ess[seq_along(.ess) %% .n == (.i %% .n)]
    .filter <- if (length(.mine) == 0L) "^$" else
      paste0("^(", paste(.mine, collapse = "|"), ")$")
  } else {
    .filter <- paste0("^(?!(", paste(.slowAll, collapse = "|"), ")$)")
  }
}
# Locally (and on CRAN) .filter stays NULL -> run everything.

test_check("nlmixr2est", stop_on_failure = FALSE, filter = .filter, perl = TRUE)
