# Fit the Panhard & Samson (2009) simulation design and summarize bias/RMSE.
#
#   Rscript design/iov-panhard/run.R <n> <nrep> <iovMethod> [outfile]
#
# e.g.  Rscript design/iov-panhard/run.R 24 1000 twoLevel results-24-twoLevel.rds
#
# Uses the INSTALLED package (library(), not load_all) so the timings and the
# optimization level are the real ones.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a)) b else a

# NLMIXR2EST_DEV=1 uses devtools::load_all() for in-development runs; the real
# reproduction runs against the installed package.
if (nzchar(Sys.getenv("NLMIXR2EST_DEV"))) {
  suppressMessages(devtools::load_all(".", quiet = TRUE, helpers = FALSE))
} else {
  suppressMessages(library(nlmixr2est))
}
source("design/iov-panhard/simulate.R")

.args <- commandArgs(trailingOnly = TRUE)
nSub <- as.integer(.args[1] %||% "24")
nRep <- as.integer(.args[2] %||% "20")
iovMethod <- .args[3] %||% "twoLevel"
outFile <- .args[4] %||% sprintf("design/iov-panhard/results-%d-%s.rds", nSub, iovMethod)

#' Pull the estimates the paper's table reports out of one fit
#'
#' Both `iovMethod` paths present the same way once `.uiFinalizeIov()` has run:
#' `$omega` is a list with the between-subject block in `$id` and the
#' inter-occasion variances in `$occ`.
panhardExtract <- function(fit, iovMethod) {
  .fx <- fit$fixef
  .om <- fit$omega$id
  .o <- fit$omega$occ
  .psi <- c(lV = .o["iov.lV", "iov.lV"],
            lKa = .o["iov.lKa", "iov.lKa"],
            lAUC = .o["iov.lAUC", "iov.lAUC"])
  c(mu.lV = .fx[["tlV"]], mu.lKa = .fx[["tlKa"]], mu.lAUC = .fx[["tlAUC"]],
    omega.lV = .om["eta.lV", "eta.lV"],
    omega.lKa = .om["eta.lKa", "eta.lKa"],
    omega.lAUC = .om["eta.lAUC", "eta.lAUC"],
    psi.lV = .psi[["lV"]], psi.lKa = .psi[["lKa"]], psi.lAUC = .psi[["lAUC"]],
    sigmaAdd = .fx[["add.sd"]], sigmaProp = .fx[["prop.sd"]])
}

panhardTrue <- c(
  mu.lV = panhardTruth$mu[["lV"]],
  mu.lKa = panhardTruth$mu[["lKa"]],
  mu.lAUC = panhardTruth$mu[["lAUC"]],
  omega.lV = panhardTruth$omega[["lV"]],
  omega.lKa = panhardTruth$omega[["lKa"]],
  omega.lAUC = panhardTruth$omega[["lAUC"]],
  psi.lV = panhardTruth$psi[["lV"]],
  psi.lKa = panhardTruth$psi[["lKa"]],
  psi.lAUC = panhardTruth$psi[["lAUC"]],
  sigmaAdd = panhardTruth$sigma,
  sigmaProp = panhardTruth$sigma)

#' Relative bias and RMSE over the usable replicates
panhardSummary <- function(res) {
  .ok <- stats::complete.cases(res)
  if (!any(.ok)) {
    return(data.frame(true = panhardTrue, biasPct = NA_real_, rmsePct = NA_real_))
  }
  .d <- sweep(res[.ok, , drop = FALSE], 2, panhardTrue, "-")
  data.frame(true = panhardTrue,
             biasPct = 100 * colMeans(.d) / panhardTrue,
             rmsePct = 100 * sqrt(colMeans(.d^2)) / abs(panhardTrue))
}

.ctl <- saemControl(nBurn = 200, nEm = 300, nmc = 3, seed = 99,
                    print = 0L, covMethod = "", calcTables = FALSE,
                    iovMethod = iovMethod)

.res <- matrix(NA_real_, nrow = nRep, ncol = length(panhardTrue),
               dimnames = list(NULL, names(panhardTrue)))
.done <- 0L
# resume: a run of 1000 replicates is hours, so pick up where a previous one
# stopped rather than starting over
if (file.exists(outFile)) {
  .prev <- try(readRDS(outFile), silent = TRUE)
  if (!inherits(.prev, "try-error") && identical(dim(.prev$estimates)[2], ncol(.res))) {
    .n <- min(nrow(.prev$estimates), nRep)
    .res[seq_len(.n), ] <- .prev$estimates[seq_len(.n), , drop = FALSE]
    .done <- .prev$done %||% .n
    cat(sprintf("resuming after %d replicates\n", .done))
  }
}
.t0 <- proc.time()
.save <- function(done) {
  saveRDS(list(n = nSub, nRep = nRep, iovMethod = iovMethod, done = done,
               estimates = .res, summary = panhardSummary(.res)), outFile)
}
for (.r in seq_len(nRep)) {
  if (.r <= .done) next
  .d <- panhardSim(nSub, seed = 1000L + .r)
  .f <- try(suppressWarnings(suppressMessages(
    nlmixr2(panhardModel(), .d, est = "saem", control = .ctl))), silent = TRUE)
  if (!inherits(.f, "try-error")) {
    .e <- try(panhardExtract(.f, iovMethod), silent = TRUE)
    if (!inherits(.e, "try-error")) .res[.r, ] <- .e
  }
  cat(sprintf("rep %d/%d  (%.1f s elapsed)\n", .r, nRep,
              (proc.time() - .t0)[["elapsed"]]))
  # checkpoint every 25 replicates so a killed run keeps its work
  if (.r %% 25L == 0L) .save(.r)
}

.summary <- panhardSummary(.res)
print(round(.summary, 2))
cat(sprintf("\n%d/%d replicates usable; %.1f s total\n",
            sum(stats::complete.cases(.res)), nRep,
            (proc.time() - .t0)[["elapsed"]]))
.save(nRep)
