## Preflight for the correlation-aware VAE covariate M-step.
##
## Two questions, answered against real data rather than assumed:
##
##   1. STATIC (no fit).  For each covariate design, is any pair of covariate
##      GROUPS correlated enough to form a colinearity cluster?  Only a cluster
##      that MERGES TWO GROUPS can activate the hysteresis and near-tie passes --
##      two shapes of one covariate share a group and are excluded by
##      construction.  If no design binds, those mechanisms cannot change any
##      existing result, whatever the threshold.
##
##   2. FITTED.  For each model, how correlated are the latent dims?  This is
##      what sets covSelectPhiJoin, and it is measured, not guessed.  Also
##      reports whether the model declares correlated etas, since the cross-phi
##      refinement is gated on that.
##
## Run from the package root:  Rscript tools/vaeColinearPreflight.R
## Add --fits to include part 2 (slow: one short fit per model).
##
## Measured 2026-08-13 (nlmixr2est 7.0.3), which is where the defaults come from:
##
##   covariate designs, max cross-group |cor|
##     theo_sd 0.000 (one group, so no cross-group pair exists)
##     neonatal_wt 0.203, wideData(30) 0.509, warfarin 0.775
##   -> nothing binds at 0.80 or above, so a 0.9 covSelectColinearCut is inert on
##      every design here.  0.80 would be uncomfortably close to warfarin's 0.775.
##
##   latent dims, theo_sd 1-cmt, 3 etas
##     corMu 0.824, corResid 0.888
##   -> 0.80 would group an ordinary well-identified model, so it is too low.
##      0.90 clears the 'suffStat'/'mu' sources with room but sits only 0.012
##      above corResid, so the 'resid' source is NOT safe to default to.

suppressMessages(library(nlmixr2est))
if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  stop("nlmixr2data is needed for the preflight datasets")
}

.cuts <- c(0.80, 0.90, 0.95)

## ---- part 1: static covariate colinearity --------------------------------

#' Largest absolute correlation between columns of two DIFFERENT covariate
#' groups.  Within-group pairs are never examined, which is exactly the rule the
#' cluster builder uses -- so this is the quantity that decides whether a
#' cluster can bind at a given cut.
maxCrossGroupCor <- function(covMat, group) {
  if (is.null(covMat) || ncol(covMat) < 2L || length(unique(group)) < 2L) {
    return(list(max = 0, pair = NA_character_))
  }
  .sd <- apply(covMat, 2, stats::sd)
  .ok <- which(is.finite(.sd) & .sd > 0)
  if (length(.ok) < 2L) return(list(max = 0, pair = NA_character_))
  .r <- abs(stats::cor(covMat[, .ok, drop = FALSE]))
  .r[!is.finite(.r)] <- 0
  .g <- group[.ok]
  .best <- 0
  .pair <- NA_character_
  .nm <- colnames(covMat)[.ok]
  for (a in seq_along(.ok)) {
    for (b in seq_len(a - 1L)) {
      if (.g[a] == .g[b]) next
      if (.r[a, b] > .best) {
        .best <- .r[a, b]
        .pair <- paste0(.nm[b], " ~ ", .nm[a])
      }
    }
  }
  list(max = .best, pair = .pair)
}

staticRow <- function(label, data) {
  d <- as.data.frame(data)
  names(d) <- toupper(names(d))
  .cov <- tryCatch(
    nlmixr2est:::.vaeCovariateSearch(d, unique(d$ID),
                                     nlmixr2est:::.vaeResolveShapes(
                                       c("power", "lin", "log", "identity",
                                         "center", "hockey"))$rules,
                                     "median", NULL, 0.05),
    error = function(e) NULL)
  if (is.null(.cov) || ncol(.cov$covMat) == 0L) {
    return(data.frame(design = label, nCol = 0L, nGroup = 0L,
                      maxCrossGroup = NA_real_, pair = NA_character_,
                      binds0.80 = FALSE, binds0.90 = FALSE, binds0.95 = FALSE))
  }
  .m <- maxCrossGroupCor(.cov$covMat, .cov$covGroup)
  data.frame(design = label, nCol = ncol(.cov$covMat),
             nGroup = length(unique(.cov$covGroup)),
             maxCrossGroup = round(.m$max, 3), pair = .m$pair,
             binds0.80 = .m$max >= .cuts[1],
             binds0.90 = .m$max >= .cuts[2],
             binds0.95 = .m$max >= .cuts[3])
}

## The wide nuisance-covariate design from test-vae-l0learn-fit.R
wideData <- function(nCov = 30L, seed = 42L) {
  set.seed(seed)
  d <- nlmixr2data::theo_sd
  ids <- unique(d$ID)
  out <- do.call(rbind, lapply(seq_len(4L), function(rep) {
    x <- d
    x$ID <- x$ID + (rep - 1L) * length(ids)
    x
  }))
  nid <- length(unique(out$ID))
  cv <- as.data.frame(matrix(stats::rnorm(nid * nCov), nid, nCov))
  names(cv) <- paste0("C", seq_len(nCov))
  cv$ID <- unique(out$ID)
  merge(out, cv, by = "ID")
}

cat("\n== Part 1: covariate colinearity (static, no fit) ==\n\n")
static <- rbind(
  staticRow("theo_sd", nlmixr2data::theo_sd),
  staticRow("neonatal_wt", nlmixr2data::neonatal_wt),
  staticRow("warfarin", nlmixr2data::warfarin),
  staticRow("wideData(30)", wideData(30L)))
print(static, row.names = FALSE)
cat("\nA cluster can only activate hysteresis/near-tie reporting where 'binds' is",
    "TRUE.\n")

## ---- part 2: latent-dim correlation (needs fits) --------------------------

if (!("--fits" %in% commandArgs(trailingOnly = TRUE))) {
  cat("\n(part 2 skipped; re-run with --fits to measure phi correlations)\n")
  quit(save = "no")
}

offDiagMax <- function(m) {
  if (is.null(m) || !is.matrix(m) || ncol(m) < 2L) return(NA_real_)
  r <- abs(stats::cor(m))
  diag(r) <- 0
  r[!is.finite(r)] <- 0
  max(r)
}

#' `returnVae = TRUE` hands back the raw training list, which is the only place
#' the posterior means and the fitted covariate centers are both exposed -- the
#' fitted object keeps only their difference, as $etaMat.
phiRow <- function(label, model, data, ctl) {
  ctl$returnVae <- TRUE
  f <- tryCatch(
    suppressMessages(suppressWarnings(
      nlmixr2est::nlmixr2(model, data, est = "vae", control = ctl))),
    error = function(e) {cat("  [", label, "failed:", conditionMessage(e), "]\n"); NULL})
  if (is.null(f) || is.null(f$mu)) return(NULL)
  om <- f$omegaMat
  data.frame(model = label,
             nEta = ncol(f$mu),
             corMu = round(offDiagMax(f$mu), 3),
             corResid = round(offDiagMax(f$mu - f$zPopMat), 3),
             omegaOffDiag = !is.null(om) && ncol(om) > 1L &&
               any(abs(om[upper.tri(om)]) > 0),
             anySelected = if (is.null(f$selected)) NA else any(f$selected))
}

cat("\n== Part 2: latent-dim correlation (fitted) ==\n\n")
ctl <- nlmixr2est::vaeControl(iters = 60L, itersBurnIn = 15L, calcTables = FALSE)

theoModel <- function() {
  ini({
    tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
    eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

phi <- do.call(rbind, list(
  phiRow("theo_sd 1cmt (3 etas)", theoModel, nlmixr2data::theo_sd, ctl)))
if (!is.null(phi)) print(phi, row.names = FALSE)

cat("\ncorMu drives covSelectPhiCor='mu'; the default 'suffStat' source is an EMA\n",
    "of the same posterior means, so it tracks corMu closely.  Set\n",
    "covSelectPhiJoin above the largest value here that should NOT group.\n", sep = "")
