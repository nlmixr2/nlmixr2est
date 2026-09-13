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


## ---------------------------------------------------------------------------
## Part 2: where to put vaeControl(covSelectPhiJoin=) -- the abs(cor) between
## LATENT DIMENSIONS at which they join a correlated group for the
## cross-parameter covariate refinement.  Different quantity from part 1, which
## is about correlation between COVARIATE COLUMNS.
##
## Measured through the shipped gate rather than by recomputing cor() here: the
## fit reports $vae$phiPairOn, the sticky adjacency the refinement actually
## used, so sweeping the join threshold and recording where a pair stops joining
## brackets the correlation the running code saw.  That also covers
## covSelectPhiCor="suffStat", whose EMA is internal and has no R-side analogue.
##
## MEASURED 2026-09-12, one-compartment oral with a declared cl/v omega block
## (45 subjects), sweeping the join threshold and reading $vae$phiPairOn:
##
##   source     joins at 0.70  0.75  0.80  0.85  0.90  0.95
##   suffStat             yes   yes    no    no    no    no
##   mu                   yes   yes    no    no    no    no
##   resid                yes   yes   yes    no    no    no
##
## So on a WELL-IDENTIFIED fit the cl/v dims sit at 0.75-0.80 under the default
## "suffStat" (and one bracket higher, 0.80-0.85, under "resid" -- which is the
## measurement behind the advice to raise the join threshold if you select it).
## A join default of 0.8 would therefore group a model with nothing wrong with
## it, which is why it is 0.9.  Re-run this before changing either default.

.phiCuts <- c(0.70, 0.75, 0.80, 0.85, 0.90, 0.95)

## one-compartment oral, clearance and volume given a correlated omega block so
## the refinement gate can open at all
.phiModel <- function() {
  ini({
    tka <- log(1.5); tcl <- log(2.7); tv <- log(31)
    eta.ka ~ 0.3
    eta.cl + eta.v ~ c(0.09, 0.07, 0.09)
    add.sd <- 0.3
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.phiSim <- function(nid = 45L, seed = 11L) {
  set.seed(seed)
  wt <- stats::runif(nid, 50, 100)
  .sig <- matrix(c(0.09, 0.075, 0.075, 0.09), 2, 2)
  z <- matrix(stats::rnorm(nid * 2L), nid, 2L) %*% chol(.sig)
  cl <- 2.7 * exp(0.7 * log(wt / mean(wt)) + z[, 1])
  v <- 31 * exp(z[, 2])
  ka <- 1.5 * exp(stats::rnorm(nid, 0, 0.3))
  tms <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
  do.call(rbind, lapply(seq_len(nid), function(i) {
    ke <- cl[i] / v[i]
    f <- 320 / v[i] * ka[i] / (ka[i] - ke) * (exp(-ke * tms) - exp(-ka[i] * tms))
    rbind(data.frame(ID = i, TIME = 0, AMT = 320, EVID = 1, DV = 0, WT = wt[i]),
          data.frame(ID = i, TIME = tms, AMT = 0, EVID = 0,
                     DV = f + stats::rnorm(length(tms), 0, 0.25), WT = wt[i]))
  }))
}

.phiJoined <- function(d, src, cut) {
  f <- try(suppressMessages(suppressWarnings(
    nlmixr2(.phiModel, d, est = "vae",
            control = vaeControl(iters = 60L, itersBurnIn = 15L,
                                 calcTables = FALSE, covSelectPhiCor = src,
                                 covSelectPhiJoin = cut,
                                 covSelectPhiLeave = cut - 0.05))))
    , silent = TRUE)
  if (inherits(f, "try-error") || is.null(f$vae$phiPairOn)) return(NA)
  .a <- f$vae$phiPairOn
  ## the cl/v pair: the two dims sharing the declared omega block
  isTRUE(.a["eta.cl", "eta.v"] == 1L)
}

.d <- .phiSim()
.phi <- do.call(rbind, lapply(c("suffStat", "mu", "resid"), function(src) {
  data.frame(source = src, cut = .phiCuts,
             joined = vapply(.phiCuts, function(cut) .phiJoined(.d, src, cut),
                             logical(1)))
}))
cat("\n-- cl/v pair joined, by covSelectPhiCor source and join threshold --\n")
print(.phi, row.names = FALSE)
cat("\n-- highest threshold at which the pair still joins --\n")
print(stats::aggregate(cut ~ source, data = .phi[.phi$joined %in% TRUE, ],
                       FUN = max), row.names = FALSE)
