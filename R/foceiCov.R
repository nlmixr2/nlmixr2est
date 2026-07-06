# Full FOCEI covariance (structural theta + residual sigma + Omega, diagonal or
# block) from the exact analytic observed information.  The augmented sensitivity
# ODE carries the 1st/2nd/3rd-order model sensitivities (rxode2 `.rxSens`); if it
# fails to solve, or the model is out of scope, the analytic path returns NULL and
# the caller falls back to the finite-difference covariance.

#' Omega blocks (connected components of the random-effect covariance), so each
#' block can be set with one `ini()` formula.  Block structure comes from the
#' DECLARED model (`idf` off-diagonal Omega rows), NOT the converged values -- a
#' declared 2x2 block whose off-diagonal converges to ~0 must stay one block so its
#' covariance parameter is retained.  Returns a list of eta-index vectors.
#' @noRd
.omegaBlocks <- function(Om, idf) {
  n <- nrow(Om); adj <- diag(TRUE, n)
  .off <- idf[!is.na(idf$neta1) & !is.na(idf$neta2) & idf$neta1 != idf$neta2, , drop = FALSE]
  for (k in seq_len(nrow(.off))) {
    a <- .off$neta1[k]; b <- .off$neta2[k]
    if (a >= 1L && a <= n && b >= 1L && b <= n) adj[a, b] <- adj[b, a] <- TRUE
  }
  comp <- integer(n); k <- 0L
  for (i in seq_len(n)) if (comp[i] == 0L) {
    k <- k + 1L; q <- i
    while (length(q)) { v <- q[1]; q <- q[-1]
      if (comp[v] == 0L) { comp[v] <- k; q <- c(q, which(adj[v, ] & comp == 0L)) } }
  }
  unname(split(seq_len(n), comp))
}

#' Free (non-fixed) Omega variance/covariance element pairs (each row `c(a, b)`,
#' `a >= b`): the lower triangle of every DECLARED Omega block, minus the fixed
#' elements.  Shared by both analytic-cov callers.
#' @noRd
.foceiOmegaPairs <- function(Om, ini) {
  blocks <- .omegaBlocks(Om, ini)
  pairs <- do.call(rbind, lapply(blocks, function(b)
    do.call(rbind, lapply(seq_along(b), function(a) cbind(b[a], b[seq_len(a)])))))
  pairs[!.omegaFixed(ini, pairs), , drop = FALSE]
}

#' Which rows of `pairs` (each `c(a, b)`) are FIXED Omega elements in `idf`.
#' @noRd
.omegaFixed <- function(idf, pairs) {
  if (nrow(pairs) == 0L) return(logical(0))
  vapply(seq_len(nrow(pairs)), function(k) {
    r <- which(!is.na(idf$neta1) & ((idf$neta1 == pairs[k, 1] & idf$neta2 == pairs[k, 2]) |
                                    (idf$neta1 == pairs[k, 2] & idf$neta2 == pairs[k, 1])))
    length(r) > 0L && isTRUE(idf$fix[r[1]])
  }, logical(1))
}

#' Diagonal Omega rows and their mu-referenced structural thetas
#'
#' @param ui an rxode2 UI object
#' @return list with `etaRows` (diagonal Omega ini rows, ordered by eta index),
#'   `etaNames`, and `thetaForEta` (the paired structural theta name, or NA for a
#'   non-mu-referenced eta)
#' @noRd
.foceiEtaThetaMap <- function(ui) {
  .idf <- ui$iniDf
  .etaRows <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2, , drop = FALSE]
  .etaRows <- .etaRows[order(.etaRows$neta1), , drop = FALSE]
  list(etaRows = .etaRows, etaNames = .etaRows$name,
       thetaForEta = ui$muRefDataFrame$theta[match(.etaRows$name, ui$muRefDataFrame$eta)])
}

#' Which of the named parameters are fixed in an ini block
#' @param idf the UI iniDf
#' @param nm character vector of parameter names
#' @return logical vector, TRUE where the parameter is fixed
#' @noRd
.iniIsFixed <- function(idf, nm) {
  if (is.null(idf$fix)) return(rep(FALSE, length(nm)))
  .i <- match(nm, idf$name)
  !is.na(.i) & !is.na(idf$fix[.i]) & idf$fix[.i]
}

#' Name the Omega variance/covariance rows of the covariance matrix
#' @param pairs matrix of eta-index pairs, each row c(a, b)
#' @param onm per-eta display name (the random-effect / eta name)
#' @return character vector, `om.<eta>` (variance) / `cov.<eta>.<eta>` (covariance)
#' @noRd
.foceiOmegaCovNames <- function(pairs, onm) {
  apply(pairs, 1, function(.pr) if (.pr[1] == .pr[2]) paste0("om.", onm[.pr[1]])
        else paste0("cov.", onm[.pr[1]], ".", onm[.pr[2]]))
}
