# Orchestrator for the analytic / finite-difference covariance paths.
#
# By default returns the fast, exact, finite-difference-free analytic
# structural-theta covariance (foceiCovAnalytic).  Omega and residual SEs are
# OPTIONAL (they require the slower full finite-difference path in the
# non-Cholesky variance-covariance parameterization, foceiCovFD).
#
# Fallback ladder (mirrors ferx's all-or-nothing behaviour, and Matthew's note
# that the augmented ODE is more likely to fail solving):
#   1. analytic structural-theta block (needs 2nd/3rd-order sensitivities);
#   2. if that is out of scope / fails -> foceiCovFD (needs only 1st-order
#      sensitivities, so it survives where the big augmented solve does not);
#   3. if the sensitivity machinery does not load at all (e.g. no symengine, or
#      rxode2 cannot build the sensitivity model) -> the fit's own built-in
#      finite-difference covariance (`fit$cov`, the Gill-Hessian covMethod result
#      -- which the covMethod="r"/"s" scaling fix in this same change set leaves
#      correctly scaled).
# So a finite-difference covariance is always produced when the sensitivities are
# unavailable; the analytic path is never partially applied.

#' The fit's built-in finite-difference (Gill-Hessian) covariance, optionally
#' restricted to a set of parameters.  Used as the final fallback when the
#' sensitivity-based paths cannot run.
#' @noRd
.foceiCovBuiltin <- function(fit, keep = NULL) {
  cv <- fit$cov
  if (is.null(cv) || !is.matrix(cv) || nrow(cv) == 0L) return(NULL)
  if (!is.null(keep)) {
    keep <- intersect(keep, rownames(cv))
    if (length(keep) == 0L) return(NULL)
    cv <- cv[keep, keep, drop = FALSE]
  }
  list(cov = cv, se = setNames(sqrt(abs(diag(cv))), rownames(cv)),
       params = rownames(cv), method = "fd-builtin")
}

#' Covariance for a fitted nlmixr2 FOCEI/FOCE object (analytic + FD fallback)
#'
#' @param fit a fitted nlmixr2 focei object.
#' @param omega logical; also return random-effect (Omega) SEs.  Requires the
#'   finite-difference path (slower).  Default `FALSE`.
#' @param residual logical; also return residual-error (sigma) SEs.  Requires
#'   the finite-difference path.  Default `FALSE`.
#' @param fallback logical; if the analytic path is out of scope / fails, fall
#'   back to the finite-difference path (and, if the sensitivities do not load at
#'   all, to the fit's built-in finite-difference covariance).  Default `TRUE`.
#' @return list with `cov`, `se`, `params`, and `method` (`"analytic"` / `"fd"` /
#'   `"fd-fallback"` / `"fd-builtin"`), or `NULL`.
#' @export
foceiCov <- function(fit, omega = FALSE, residual = FALSE, fallback = TRUE) {
  # structural (mu-referenced) theta names, in eta order
  muRef <- fit$finalUi$muRefDataFrame
  ini <- fit$finalUi$iniDf
  etaRows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2, , drop = FALSE]
  etaRows <- etaRows[order(etaRows$neta1), , drop = FALSE]
  thN <- muRef$theta[match(etaRows$name, muRef$eta)]

  if (isTRUE(omega) || isTRUE(residual)) {
    # full covariance (theta + Omega + residual), non-Cholesky variance scale
    r <- foceiCovFD(fit)
    if (!is.null(r)) { r$method <- "fd"; return(r) }
    if (!isTRUE(fallback)) return(NULL)
    return(.foceiCovBuiltin(fit))                 # sens unavailable -> built-in FD
  }
  # fast path: analytic structural-theta block
  r <- foceiCovAnalytic(fit)
  if (!is.null(r)) { r$method <- "analytic"; return(r) }
  if (!isTRUE(fallback)) return(NULL)
  # fallback 1: finite-difference (1st-order sensitivities), structural-theta block
  rf <- foceiCovFD(fit)
  if (!is.null(rf)) {
    thNk <- intersect(thN, rf$params)
    return(list(cov = rf$cov[thNk, thNk, drop = FALSE], se = rf$se[thNk],
                params = thNk, method = "fd-fallback"))
  }
  # fallback 2: sensitivities did not load at all -> the fit's built-in FD covariance
  .foceiCovBuiltin(fit, thN)
}
