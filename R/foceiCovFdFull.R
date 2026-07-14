# covType="fd", covFull=TRUE: the full theta+sigma+Omega covariance by finite differences,
# computed in C++ (foceiCalcRFdFull) at the FD seam.  C++ stashes the full Hessian inverse
# `.fdFullCov` (= Rinv_full) and, for the S-using methods, the full OPG cross-product
# `.fdFullS` (= Sfull).  These helpers are the R glue: the enumeration/names the C++ loop
# needs, and the install of the full cov as the fit's $cov per covMethod -- "r" -> Rinv_full,
# "s" -> solve(Sfull), "r,s" -> the true sandwich Rinv_full %*% Sfull %*% Rinv_full (the FD
# counterpart to .foceiInstallAnalyticCov).

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

#' Install the C++ FD-full covariance as `fit$cov` (and `fit$covR/covS/covRS`) when
#' `covFull = TRUE`.  Routes on the fit's final `covMethod` string: an "r,s"-family cov
#' installs the true full sandwich `Rinv_full %*% Sfull %*% Rinv_full`, an "s"-family cov
#' installs `solve(Sfull)`, and an "r"-family cov installs `Rinv_full` (the Hessian inverse).
#' No-op (native theta-only cov kept) if the stashed pieces are absent/non-finite, if the
#' assembled cov is not positive-definite, or for a non-FD covMethod (analytic/failed/
#' boundary).  The FD counterpart to [.foceiInstallAnalyticCov].
#' @param .ret focei fit environment
#' @noRd
.foceiInstallFdFullCov <- function(.ret) {
  if (!exists(".fdFullCov", envir = .ret, inherits = FALSE)) return(invisible())
  .cm <- if (exists("covMethod", envir = .ret, inherits = FALSE)) .ret$covMethod else ""
  if (length(.cm) != 1L || is.na(.cm)) .cm <- ""
  # native covMethod strings: r|r+||r| , s|s+||s| , and "<r>,<s>" for the sandwich.
  .type <- if (grepl("^(r\\+?|\\|r\\|),(s\\+?|\\|s\\|)$", .cm)) "r,s"
    else if (grepl("^(r\\+?|\\|r\\|)$", .cm)) "r"
    else if (grepl("^(s\\+?|\\|s\\|)$", .cm)) "s"
    else return(invisible())   # analytic / failed / "" / boundary -> keep the native cov
  .Rinv <- get(".fdFullCov", envir = .ret)
  if (!is.matrix(.Rinv) || !all(is.finite(.Rinv))) return(invisible())
  .S <- if (exists(".fdFullS", envir = .ret, inherits = FALSE)) get(".fdFullS", envir = .ret) else NULL
  if (.type != "r" && (!is.matrix(.S) || !all(is.finite(.S)))) return(invisible())
  .covS <- if (is.null(.S)) NULL else tryCatch(solve(.S), error = function(e) NULL)
  if (.type != "r" && is.null(.covS)) return(invisible())
  .covRS <- if (is.null(.S)) NULL else .Rinv %*% .S %*% .Rinv
  .cov <- switch(.type, "r" = .Rinv, "s" = .covS, "r,s" = .covRS)
  if (is.null(.cov) || !is.matrix(.cov) || !all(is.finite(.cov))) return(invisible())
  # keep the var-cov dimnames on the assembled products
  dimnames(.cov) <- dimnames(.Rinv)
  # PD guard: an indefinite assembled cov installs negative variances -> NaN SEs.  Reject and
  # keep the native theta-only cov rather than a plausible-looking wrong one.
  .ev <- suppressWarnings(eigen(.cov, symmetric = TRUE, only.values = TRUE)$values)
  if (any(diag(.cov) <= 0) || !all(is.finite(.ev)) || min(.ev) <= 0) return(invisible())
  .ret$cov <- .cov
  .ret$covR <- .Rinv
  if (!is.null(.covS)) {
    dimnames(.covS) <- dimnames(.Rinv)
    .ret$covS <- .covS
  }
  if (!is.null(.covRS)) {
    dimnames(.covRS) <- dimnames(.Rinv)
    .ret$covRS <- .covRS
  }
  .foceiCovCondition(.ret, .cov, .ev)
  invisible()
}
