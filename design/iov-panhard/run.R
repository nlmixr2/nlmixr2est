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
panhardExtract <- function(fit, iovMethod) {
  .fx <- fit$fixef
  .om <- if (is.list(fit$omega)) fit$omega$id else fit$omega
  .psi <- if (identical(iovMethod, "twoLevel")) {
    c(lV = .om["rx.iov.lV.1", "rx.iov.lV.1"],
      lKa = .om["rx.iov.lKa.1", "rx.iov.lKa.1"],
      lAUC = .om["rx.iov.lAUC.1", "rx.iov.lAUC.1"])
  } else {
    .o <- fit$omega$occ
    c(lV = .o["iov.lV", "iov.lV"],
      lKa = .o["iov.lKa", "iov.lKa"],
      lAUC = .o["iov.lAUC", "iov.lAUC"])
  }
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

.ctl <- saemControl(nBurn = 200, nEm = 300, nmc = 3, seed = 99,
                    print = 0L, covMethod = "", calcTables = FALSE,
                    iovMethod = iovMethod)

.res <- matrix(NA_real_, nrow = nRep, ncol = length(panhardTrue),
               dimnames = list(NULL, names(panhardTrue)))
.t0 <- proc.time()
for (.r in seq_len(nRep)) {
  .d <- panhardSim(nSub, seed = 1000L + .r)
  .f <- try(suppressWarnings(suppressMessages(
    nlmixr2(panhardModel(), .d, est = "saem", control = .ctl))), silent = TRUE)
  if (!inherits(.f, "try-error")) {
    .e <- try(panhardExtract(.f, iovMethod), silent = TRUE)
    if (!inherits(.e, "try-error")) .res[.r, ] <- .e
  }
  cat(sprintf("rep %d/%d  (%.1f s elapsed)\n", .r, nRep,
              (proc.time() - .t0)[["elapsed"]]))
}

.ok <- stats::complete.cases(.res)
.bias <- 100 * colMeans(sweep(.res[.ok, , drop = FALSE], 2, panhardTrue, "-")) / panhardTrue
.rmse <- 100 * sqrt(colMeans(sweep(.res[.ok, , drop = FALSE], 2, panhardTrue, "-")^2)) /
  abs(panhardTrue)
.summary <- data.frame(true = panhardTrue, biasPct = .bias, rmsePct = .rmse)
print(round(.summary, 2))
cat(sprintf("\n%d/%d replicates usable; %.1f s total\n",
            sum(.ok), nRep, (proc.time() - .t0)[["elapsed"]]))
saveRDS(list(n = nSub, nRep = nRep, iovMethod = iovMethod,
             estimates = .res, summary = .summary), outFile)
