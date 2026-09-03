# Tabulate a Panhard & Samson reproduction run against the paper's own numbers.
#
#   Rscript design/iov-panhard/report.R design/iov-panhard/results-24-twoLevel.rds ...
#
# Each argument is a results file written by run.R.  Prints one relative-bias
# and one RMSE table with a column per run.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a)) b else a

.files <- commandArgs(trailingOnly = TRUE)
if (length(.files) == 0L) {
  .files <- list.files("design/iov-panhard", pattern = "^results-.*[.]rds$",
                       full.names = TRUE)
}
if (length(.files) == 0L) stop("no results files given or found", call. = FALSE)

.res <- lapply(.files, readRDS)
.lab <- vapply(.res, function(x) sprintf("n=%d %s", x$n, x$iovMethod),
               character(1))

.col <- function(what) {
  .m <- do.call(cbind, lapply(.res, function(x) x$summary[[what]]))
  dimnames(.m) <- list(rownames(.res[[1]]$summary), .lab)
  round(.m, 2)
}

cat("True values:\n")
print(.res[[1]]$summary["true"])
cat("\nReplicates used:\n")
print(setNames(vapply(.res, function(x) sum(stats::complete.cases(x$estimates)),
                      integer(1)), .lab))
cat("\nRelative bias (%):\n")
print(.col("biasPct"))
cat("\nRelative RMSE (%):\n")
print(.col("rmsePct"))

# Panhard & Samson (2009) Table 1, the extended-SAEM columns, transcribed from
# arXiv:0803.4437v1.  Rows are named to match this harness's summary rows so the
# two tables can be printed side by side.
.paper <- data.frame(
  row.names = c("mu.lV", "mu.lKa", "mu.lAUC",
                "omega.lV", "omega.lKa", "omega.lAUC",
                "psi.lV", "psi.lKa", "psi.lAUC"),
  bias24 = c(0.01, 0.48, -0.08, -5.13, -3.99, -4.88, -8.67, -10.94, -5.37),
  bias40 = c(-0.06, 0.02, -0.11, -3.45, -3.23, -1.51, -5.93, -7.06, -4.92),
  rmse24 = c(3.9, 14.4, 1.0, 38.7, 42.4, 34.5, 69.4, 73.5, 43.6),
  rmse40 = c(2.91, 10.79, 0.79, 30.30, 33.49, 27.41, 58.78, 62.00, 33.31))

.side <- function(what, pcols) {
  .m <- .col(what)
  .r <- intersect(rownames(.paper), rownames(.m))
  .out <- .m[.r, , drop = FALSE]
  for (.i in seq_along(pcols)) {
    .out <- cbind(.out, .paper[.r, pcols[[.i]]])
    colnames(.out)[ncol(.out)] <- names(pcols)[.i]
  }
  .out
}
.want <- c("paper n=24" = "bias24", "paper n=40" = "bias40")
cat("\nRelative bias (%) vs Panhard & Samson Table 1 (extended SAEM):\n")
print(.side("biasPct", as.list(.want)))
.want <- c("paper n=24" = "rmse24", "paper n=40" = "rmse40")
cat("\nRelative RMSE (%) vs Panhard & Samson Table 1 (extended SAEM):\n")
print(.side("rmsePct", as.list(.want)))

cat("\nTwo rows of the paper's table are deliberately not reproduced:\n",
    "* beta_V, beta_ka, beta_AUC -- the arXiv HTML never states the beta\n",
    "  values used to generate the data, so this harness simulates under H0\n",
    "  (beta = 0), where a relative bias in beta is undefined.\n",
    "* sigma^2 -- the paper ties one sigma to g = 1 + f, while this harness\n",
    "  estimates add.sd and prop.sd separately (combined1).  Estimating two\n",
    "  parameters where the paper estimates one makes each individually\n",
    "  noisier, so the sigma rows are not comparable.\n", sep = "")
