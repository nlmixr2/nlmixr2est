# Full FOCEI covariance (structural theta + residual sigma + Omega, diagonal or block)
# with a two-tier ladder.  Both tiers use the SAME direction-set model and the SAME
# R-matrix assembly, differing only in how the 3rd-order tensor is obtained:
#
#   1. fd2 (DEFAULT): analytic 2nd-order sensitivities plus Shi(2021) adaptive finite
#      differences of those EXACT 2nd-order sensitivities for the 3rd-order term.  It
#      reproduces the exact tier to ~1e-3 (it differences exact derivatives, not the
#      objective, so no catastrophic cancellation) AND avoids building the large
#      3rd-order augmented model, whose C compile is the dominant cost: for a
#      covariate model the 3rd-order model is ~80 KB and takes ~40 s to compile vs
#      ~4 s for the 2nd-order model, so fd2 is ~10x faster to build.
#   2. exact3 (FALLBACK): analytic exact 3rd-order sensitivities, used only when the
#      2nd-order finite-difference step fails (non-finite), so the exact tier is
#      always available but rarely paid for.
#
# Both tiers return the SAME full covariance (theta + sigma + Omega, diagonal or
# block), so the fallback never silently drops the Omega/residual blocks.

#' Omega blocks (connected components of the random-effect covariance)
#'
#' Each connected component is a block that can be set with one `ini()` formula.
#' @param Om the estimated Omega matrix
#' @param tol off-diagonal magnitude treated as a non-zero link
#' @return list of eta-index vectors, one per block
#' @author Hidde van de Beek
#' @noRd
.omegaBlocks <- function(Om, tol = 1e-10) {
  .n <- nrow(Om)
  .adj <- abs(Om) > tol
  diag(.adj) <- TRUE
  .comp <- integer(.n)
  .k <- 0L
  for (.i in seq_len(.n)) {
    if (.comp[.i] == 0L) {
      .k <- .k + 1L
      .q <- .i
      while (length(.q)) {
        .v <- .q[1]
        .q <- .q[-1]
        if (.comp[.v] == 0L) {
          .comp[.v] <- .k
          .q <- c(.q, which(.adj[.v, ] & .comp == 0L))
        }
      }
    }
  }
  unname(split(seq_len(.n), .comp))
}

#' Which Omega element pairs are fixed
#'
#' @param idf the UI iniDf
#' @param pairs matrix of eta-index pairs, each row c(a, b)
#' @return logical vector, TRUE where the pair is a fixed Omega element
#' @author Hidde van de Beek
#' @noRd
.omegaFixed <- function(idf, pairs) {
  if (nrow(pairs) == 0L) {
    return(logical(0))
  }
  vapply(seq_len(nrow(pairs)), function(.k) {
    .r <- which(!is.na(idf$neta1) &
                  ((idf$neta1 == pairs[.k, 1] & idf$neta2 == pairs[.k, 2]) |
                     (idf$neta1 == pairs[.k, 2] & idf$neta2 == pairs[.k, 1])))
    length(.r) > 0L && isTRUE(idf$fix[.r[1]])
  }, logical(1))
}

#' Diagonal Omega rows and their mu-referenced structural thetas
#'
#' @param ui an rxode2 UI object
#' @return list with `etaRows` (the diagonal Omega ini rows, ordered by eta
#'   index), `etaNames`, and `thetaForEta` (the paired structural theta name, or
#'   NA for a non-mu-referenced eta)
#' @author Hidde van de Beek
#' @noRd
.foceiEtaThetaMap <- function(ui) {
  .idf <- ui$iniDf
  .etaRows <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2, , drop = FALSE]
  .etaRows <- .etaRows[order(.etaRows$neta1), , drop = FALSE]
  .muRef <- ui$muRefDataFrame
  list(etaRows = .etaRows,
       etaNames = .etaRows$name,
       thetaForEta = .muRef$theta[match(.etaRows$name, .muRef$eta)])
}

#' Which of the named parameters are fixed in an ini block
#'
#' @param idf the UI iniDf
#' @param nm character vector of parameter names
#' @return logical vector, TRUE where the parameter is fixed
#' @author Hidde van de Beek
#' @noRd
.iniIsFixed <- function(idf, nm) {
  if (is.null(idf$fix)) {
    return(rep(FALSE, length(nm)))
  }
  .i <- match(nm, idf$name)
  !is.na(.i) & !is.na(idf$fix[.i]) & idf$fix[.i]
}

#' Name the Omega variance/covariance rows of the covariance matrix
#'
#' @param pairs matrix of eta-index pairs, each row c(a, b)
#' @param onm per-eta display name (the mu-referenced theta, or the eta name when
#'   not mu-referenced)
#' @return character vector, `om.<name>` for a variance and `cov.<name>.<name>`
#'   for a covariance
#' @author Hidde van de Beek
#' @noRd
.foceiOmegaCovNames <- function(pairs, onm) {
  apply(pairs, 1, function(.pr) {
    if (.pr[1] == .pr[2]) {
      paste0("om.", onm[.pr[1]])
    } else {
      paste0("cov.", onm[.pr[1]], ".", onm[.pr[2]])
    }
  })
}

#' Full covariance for a fitted nlmixr2 FOCEI object
#'
#' Returns the full observed-information covariance over the structural thetas,
#' residual sigma, and Omega (diagonal or block) variances and covariances, via a
#' two-tier ladder: fd2 (the default, 2nd-order model with Shi finite differences)
#' then exact3 (the 3rd-order fallback).  Each tier returns the same full set of
#' parameters; the analytic path is never partially applied.
#' @param fit a fitted nlmixr2 focei object
#' @param covFull FALSE for the structural-theta block only, TRUE for theta +
#'   sigma + Omega
#' @return list with `cov`, `se`, `params`, and `method` (`"analytic-fd2"` /
#'   `"analytic"`), or NULL if no tier succeeds
#' @author Hidde van de Beek
#' @noRd
.foceiCov <- function(fit, covFull = FALSE) {
  # each tier is guarded: an unexpected error (e.g. a singular Omega solve) drops
  # to the next.  fd2 is the default: same result to ~1e-3 as exact3 but ~10x
  # faster to build (it skips compiling the large 3rd-order augmented model).
  # exact3 is the exact fallback.
  .r <- tryCatch(.foceiCovAnalytic(fit, sens = "fd2", covFull = covFull),
                 error = function(e) NULL)
  if (!is.null(.r)) {
    .r$method <- "analytic-fd2"
    return(.r)
  }
  .r <- tryCatch(.foceiCovAnalytic(fit, sens = "exact3", covFull = covFull),
                 error = function(e) NULL)
  if (!is.null(.r)) {
    .r$method <- "analytic"
    return(.r)
  }
  NULL
}
