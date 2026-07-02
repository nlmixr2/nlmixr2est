# Full FOCEI/FOCE covariance (structural theta + residual sigma + Omega, diagonal or block)
# with a two-tier robustness ladder (Matthew's note that the augmented ODE is more
# likely to fail solving the higher its order):
#
#   1. analytic, EXACT 3rd-order sensitivities      (foceiCovAnalytic sens="exact3")
#   2. analytic, 2nd-order + Shi(2021) finite differences of the 2nd-order
#      sensitivities for the 3rd-order term         (foceiCovAnalytic sens="fd2")
#        -- a much lighter ODE (no O(neta^3) states) that solves where tier 1
#           does not, and reproduces it to ~1e-5 (it differences EXACT
#           derivatives, not the objective, so no catastrophic cancellation).
#
# Both tiers return the SAME full covariance (theta + sigma + Omega, diagonal or
# block), so the fallback never silently drops the Omega/residual blocks.

#' Omega blocks (connected components of the random-effect covariance), so each
#' block can be set with one `ini()` formula.  Returns a list of eta-index vectors.
#' @noRd
.omegaBlocks <- function(Om, tol = 1e-10) {
  n <- nrow(Om); adj <- abs(Om) > tol; diag(adj) <- TRUE
  comp <- integer(n); k <- 0L
  for (i in seq_len(n)) if (comp[i] == 0L) {
    k <- k + 1L; q <- i
    while (length(q)) { v <- q[1]; q <- q[-1]
      if (comp[v] == 0L) { comp[v] <- k; q <- c(q, which(adj[v, ] & comp == 0L)) } }
  }
  unname(split(seq_len(n), comp))
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

#' Full covariance for a fitted nlmixr2 FOCEI/FOCE object
#'
#' Returns the full observed-information covariance over the structural thetas,
#' residual sigma, and Omega (diagonal or block) variances and covariances, via a
#' two-tier ladder: analytic exact 3rd-order, then 2nd-order with Shi finite
#' differences.  Each tier returns the same full set of parameters; the analytic
#' path is never partially applied.
#'
#' @param fit a fitted nlmixr2 focei object.
#' @return list with `cov`, `se`, `params`, and `method` (`"analytic"` /
#'   `"analytic-fd2"`), or `NULL` if no tier succeeds.
#' @noRd
foceiCov <- function(fit, covFull = FALSE) {
  # each tier guarded: an unexpected error (e.g. singular Omega solve) drops to the next
  r <- tryCatch(foceiCovAnalytic(fit, sens = "exact3", covFull = covFull), error = function(e) NULL)
  if (!is.null(r)) { r$method <- "analytic"; return(r) }
  r <- tryCatch(foceiCovAnalytic(fit, sens = "fd2", covFull = covFull), error = function(e) NULL)
  if (!is.null(r)) { r$method <- "analytic-fd2"; return(r) }
  NULL
}
