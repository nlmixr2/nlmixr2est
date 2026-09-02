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

cat("\nPaper (Table 1, SAEM, n=24) for comparison -- see the note in\n",
    "simulate.R: the arXiv HTML does not state the beta values used to\n",
    "generate the data, so this harness runs under H0 and the beta rows of\n",
    "the paper's table are not reproduced here.\n", sep = "")
