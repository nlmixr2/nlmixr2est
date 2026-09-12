## Where to put vaeControl(covSelectColinearCut=) -- the abs(cor) at which two
## covariates are treated as near-interchangeable (R/vaeCovColinear.R).
##
## Clustering is a LABEL, never a constraint: it only lets the covariate M-step
## keep last iteration's incumbent when a mate does not beat it by a covariate's
## L0 cost, and report the mates that came close.  So the cost of a cut that is
## too low is a needlessly sticky selection, and the cost of one that is too high
## is that the chattering it exists to stop goes unaddressed.  What this script
## measures is where, on real designs, groups actually start merging.
##
## Run: Rscript tools/vaeColinearPreflight.R
##
## MEASURED 2026-09-11 over the 20 nlmixr2data datasets that yield more than
## one covariate group.  Number of datasets in which a cluster MERGES two groups:
##
##   cut    0.50  0.60  0.70  0.80  0.85  0.90  0.95  0.99
##   binds     6     5     3     2     2     2     2     2
##
## The count is flat from 0.80 all the way to 0.99 and only the two datasets with
## genuinely near-duplicate covariates (pump, rats) are in it; below 0.80 the
## count climbs as ordinarily-correlated covariates (warfarin's WT/AGE/SEX,
## wbcSim's V2I/V1I/CLI) start merging.  So anywhere in [0.80, 0.99] behaves the
## same on this corpus, and 0.90 sits in the middle of that plateau: high enough
## not to label merely-correlated covariates interchangeable, low enough to still
## catch the duplicates.  Re-run this before changing the default.

library(nlmixr2est)

.cuts <- c(0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.99)

## Every nlmixr2data dataset vaeCovariates() will take.  A dataset with one
## covariate can never merge two groups, so it is reported but uninformative.
.dataSets <- function() {
  .nm <- utils::data(package = "nlmixr2data")$results[, "Item"]
  .out <- list()
  for (.n in .nm) {
    ## a lazy-loaded dataset is not an export, so getFromNamespace() misses it
    .d <- try(get(.n, envir = asNamespace("nlmixr2data")), silent = TRUE)
    if (inherits(.d, "try-error") || !is.data.frame(.d)) next
    .out[[.n]] <- .d
  }
  .out
}

## vaeCovariates() builds exactly the design the search sees, so measuring
## through it measures what the M-step would measure.
.row <- function(name, d) {
  .res <- try(vaeCovariates(d, warn = FALSE), silent = TRUE)
  if (inherits(.res, "try-error") || nrow(.res) == 0L) return(NULL)
  .per <- lapply(.cuts, function(cut) vaeCovariates(d, warn = FALSE,
                                                    colinearCut = cut))
  data.frame(data = name, groups = length(unique(.res$group)), cut = .cuts,
             clusters = vapply(.per, function(r) length(unique(r$cluster)),
                               integer(1)),
             binds = vapply(.per, function(r) .vaeClusterBinds(r$cluster,
                                                               r$group),
                            logical(1)))
}

.dat <- .dataSets()
.all <- do.call(rbind, lapply(names(.dat), function(n) .row(n, .dat[[n]])))
if (is.null(.all)) {
  cat("no dataset produced covariates\n")
} else {
  .multi <- .all[.all$groups > 1L, ]
  cat("\n-- datasets with more than one covariate group --\n")
  print(.multi, row.names = FALSE)
  cat("\n-- how many datasets bind at each cut --\n")
  print(stats::aggregate(binds ~ cut, data = .multi, FUN = sum), row.names = FALSE)
}
