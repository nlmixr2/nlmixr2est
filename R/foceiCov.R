# Full FOCEI covariance (theta + residual sigma + Omega): analytic 2nd-order
# sensitivities with Shi(2021) finite differences for the 3rd-order term, else the
# caller falls back to the finite-difference covariance.

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
#' residual sigma, and Omega variances and covariances, from the analytic 2nd-order
#' sensitivities with Shi finite differences for the 3rd-order term.  The analytic
#' path is never partially applied: it returns the full parameter set or NULL.
#' @param fit a fitted nlmixr2 focei object
#' @param covFull FALSE for the structural-theta block only, TRUE for theta +
#'   sigma + Omega
#' @return list with `cov`, `se`, `params`, and `method` (`"analytic-fd2"`), or NULL
#' @author Hidde van de Beek
#' @noRd
.foceiCov <- function(fit, covFull = FALSE) {
  # Single analytic tier (fd2): 2nd-order sensitivities with Shi(2021) finite
  # differences for the 3rd-order term.  Any failure returns NULL and the caller
  # falls back to the finite-difference covariance.
  .r <- tryCatch(.foceiCovAnalytic(fit, covFull = covFull), error = function(e) NULL)
  if (!is.null(.r)) {
    .r$method <- "analytic-fd2"
    return(.r)
  }
  NULL
}
