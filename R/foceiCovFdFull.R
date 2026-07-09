# covType="fd", covFull=TRUE: the full theta+sigma+Omega covariance by finite differences
# of the objective, computed in C++ (foceiCalcRFdFull) at the FD seam.  These two helpers
# are the R glue: the enumeration/names the C++ loop needs, and the install of the stashed
# natural cov as the fit's $cov (the FD counterpart to .foceiInstallAnalyticCov).

#' Parameter enumeration for the C++ FD-full covariance (foceiCalcRFdFull).
#'
#' From the live cov-step environment `e`, returns the free structural + residual theta
#' positions (0-based, into `op_focei.fullTheta`) and the free `Omega` lower-triangle
#' elements (1-based row/col), plus matching natural-scale parameter names (`om.<theta>` /
#' `cov.<theta>.<theta>` as in [foceiCovAnalytic]).  `NULL` if there is nothing to do.
#' @param e focei cov-step environment
#' @noRd
.foceiFdFullParams <- function(e) {
  ui <- get("ui", e)
  ini <- ui$iniDf
  isFix <- if (is.null(ini$fix)) rep(FALSE, nrow(ini)) else ini$fix
  isFix[is.na(isFix)] <- FALSE
  thRows <- which(!is.na(ini$ntheta) & !isFix)                  # structural + residual thetas
  if (length(thRows) == 0L) return(NULL)
  thPos <- as.integer(ini$ntheta[thRows] - 1L)                  # 0-based position in fullTheta
  thNames <- ini$name[thRows]
  Om <- get("omega", e)
  pairs <- .foceiOmegaPairs(Om, ini)                            # free Omega lower-triangle (a>=b)
  if (is.null(pairs) || nrow(pairs) == 0L) {
    return(list(thPos = thPos, omA = integer(0), omB = integer(0), names = thNames))
  }
  map <- .foceiEtaThetaMap(ui)
  onm <- map$etaNames                                          # Omega named by the eta
  omNames <- .foceiOmegaCovNames(pairs, onm)
  list(thPos = thPos, omA = as.integer(pairs[, 1]), omB = as.integer(pairs[, 2]),
       names = c(thNames, omNames))
}

#' Install the C++ FD-full covariance (`e[".fdFullCov"]`) as `fit$cov` when
#' `covType = "fd"` and `covFull = TRUE`.  No-op if it was not computed, or if it is not
#' finite / positive-definite (the native theta-only FD cov is then kept).  The FD
#' counterpart to [.foceiInstallAnalyticCov].
#' @param .ret focei fit environment
#' @noRd
.foceiInstallFdFullCov <- function(.ret) {
  if (!exists(".fdFullCov", envir = .ret, inherits = FALSE)) return(invisible())
  .cov <- get(".fdFullCov", envir = .ret)
  if (!is.matrix(.cov) || !all(is.finite(.cov))) return(invisible())
  .ev <- suppressWarnings(eigen(.cov, symmetric = TRUE, only.values = TRUE)$values)
  if (any(diag(.cov) <= 0) || !all(is.finite(.ev)) || min(.ev) <= 0) return(invisible())
  .ret$cov <- .cov
  invisible()
}
