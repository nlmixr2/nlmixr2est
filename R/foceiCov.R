# Orchestrator for the analytic / finite-difference covariance paths.
#
# By default returns the fast, exact, finite-difference-free analytic
# structural-theta covariance (foceiCovAnalytic).  Omega and residual SEs are
# OPTIONAL (they require the slower full finite-difference path in the
# non-Cholesky variance-covariance parameterization, foceiCovFD).  If the
# analytic path is out of scope or fails, falls back to the finite-difference
# path (mirrors ferx's all-or-nothing fallback).

#' Covariance for a fitted nlmixr2 FOCEI/FOCE object (analytic + FD fallback)
#'
#' @param fit a fitted nlmixr2 focei object.
#' @param omega logical; also return random-effect (Omega) SEs.  Requires the
#'   finite-difference path (slower).  Default `FALSE`.
#' @param residual logical; also return residual-error (sigma) SEs.  Requires
#'   the finite-difference path.  Default `FALSE`.
#' @param fallback logical; if the analytic path is out of scope / fails, fall
#'   back to the finite-difference path.  Default `TRUE`.
#' @return list with `cov`, `se`, `params`, and `method`
#'   (`"analytic"` / `"fd"` / `"fd-fallback"`), or `NULL`.
#' @export
foceiCov <- function(fit, omega = FALSE, residual = FALSE, fallback = TRUE) {
  if (isTRUE(omega) || isTRUE(residual)) {
    # full covariance (theta + Omega + residual), non-Cholesky variance scale
    r <- foceiCovFD(fit)
    if (is.null(r)) return(NULL)
    r$method <- "fd"
    return(r)
  }
  # fast path: analytic structural-theta block
  r <- foceiCovAnalytic(fit)
  if (!is.null(r)) { r$method <- "analytic"; return(r) }
  if (!isTRUE(fallback)) return(NULL)
  # fallback: finite-difference, restricted to the structural-theta block
  rf <- foceiCovFD(fit)
  if (is.null(rf)) return(NULL)
  muRef <- fit$finalUi$muRefDataFrame
  ini <- fit$finalUi$iniDf
  etaRows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2, , drop = FALSE]
  etaRows <- etaRows[order(etaRows$neta1), , drop = FALSE]
  thN <- muRef$theta[match(etaRows$name, muRef$eta)]
  thN <- intersect(thN, rf$params)
  list(cov = rf$cov[thN, thN, drop = FALSE], se = rf$se[thN], params = thN, method = "fd-fallback")
}
