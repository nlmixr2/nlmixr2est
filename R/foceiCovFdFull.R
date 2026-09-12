# covFull=TRUE finite-difference full theta+sigma+Omega covariance: C++
# (foceiCalcRFdFull) stashes the Hessian inverse `.fdFullCov` (Rinv_full) and the OPG
# cross-product `.fdFullS` (Sfull); these helpers enumerate the parameters and install
# the cov per covMethod (the FD counterpart to .foceiInstallAnalyticCov).

#' Parameter enumeration for the C++ FD-full covariance (foceiCalcRFdFull).
#'
#' From the cov-step environment `e`, the free structural+residual theta positions
#' (0-based) and free `Omega` lower-triangle elements (1-based), plus matching
#' `om.<theta>` / `cov.<theta>.<theta>` names.  `NULL` if there is nothing to do.
#' @param e focei cov-step environment
#' @noRd
.foceiFdFullParams <- function(e) {
  ui <- get("ui", e)
  # bounded-parameter transforms put thetas on an internal scale the theta-sized
  # Jacobian hook cannot correct for a full cov; bow out and keep the native cov.
  if (!is.null(ui$boundedTransforms) && length(ui$boundedTransforms) > 0L) return(NULL)
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
#' `covFull = TRUE`, routing on the fit's `covMethod`: "r,s" -> the sandwich
#' `Rinv %*% S %*% Rinv`, "s" -> `solve(S)`, "r" -> `Rinv`.  No-op (native cov kept)
#' if the pieces are absent/non-finite, the cov is not positive-definite, or covMethod
#' is not an FD method.  FD counterpart to [.foceiInstallAnalyticCov].
#' @param .ret focei fit environment
#' @noRd
.foceiInstallFdFullCov <- function(.ret) {
  if (!exists(".fdFullCov", envir = .ret, inherits = FALSE)) return(invisible(FALSE))
  # The fit env's covMethod records what the NATIVE theta-only step produced, not what was
  # asked for -- C++ downgrades it to "s" when that step's "r" fails.  The full R/S pieces
  # used below are computed independently of that step, so route on the REQUESTED control;
  # otherwise a requested "r,s" silently installs solve(S) with a usable .Rinv in hand.
  # foceiControl() stores covMethod as an integer code ("r,s"=1, "r"=2, "s"=3, ""=0) with
  # covType separating "r" from "analytic", so decode it rather than reading a string.
  .env <- if (exists("covMethod", envir = .ret, inherits = FALSE)) .ret$covMethod else ""
  if (length(.env) != 1L || !is.character(.env) || is.na(.env)) .env <- ""
  .cm <- .env
  .code <- tryCatch(rxode2::rxGetControl(.ret$ui, "covMethod", NA_integer_),
                    error = function(e) NA_integer_)
  .cty <- tryCatch(rxode2::rxGetControl(.ret$ui, "covType", "fd"), error = function(e) "fd")
  if (is.numeric(.code) && length(.code) == 1L && !is.na(.code) &&
      !identical(.cty, "analytic")) {
    .req <- switch(as.character(as.integer(.code)), "1" = "r,s", "2" = "r", "3" = "s", "")
    if (nzchar(.req)) .cm <- .req
  }
  .type <- .covFdType(.cm)
  if (!nzchar(.type)) return(invisible(FALSE))  # analytic / failed / "" / boundary -> keep native
  .Rinv <- get(".fdFullCov", envir = .ret)
  if (!is.matrix(.Rinv) || !all(is.finite(.Rinv))) return(invisible(FALSE))
  .S <- if (exists(".fdFullS", envir = .ret, inherits = FALSE)) get(".fdFullS", envir = .ret) else NULL
  if (.type != "r" && (!is.matrix(.S) || !all(is.finite(.S)))) return(invisible(FALSE))
  .covS <- if (is.null(.S)) NULL else tryCatch(solve(.S), error = function(e) NULL)
  if (.type != "r" && is.null(.covS)) return(invisible(FALSE))
  .covRS <- if (is.null(.S)) NULL else .Rinv %*% .S %*% .Rinv
  .cov <- switch(.type, "r" = .Rinv, "s" = .covS, "r,s" = .covRS)
  if (is.null(.cov) || !is.matrix(.cov) || !all(is.finite(.cov))) return(invisible(FALSE))
  dimnames(.cov) <- dimnames(.Rinv)
  # PD guard: reject an indefinite cov (negative variances -> NaN SEs), keep the native cov.
  .ev <- suppressWarnings(eigen(.cov, symmetric = TRUE, only.values = TRUE)$values)
  if (any(diag(.cov) <= 0) || !all(is.finite(.ev)) || min(.ev) <= 0) return(invisible(FALSE))
  # The theta-only covariance the native step produced -- and the r/s/sandwich pieces
  # behind it -- are about to be replaced.  Cache them first so setCov() can swap back
  # to the theta-only shape without recomputing anything (they are already in hand).
  .nat <- lapply(stats::setNames(c("covR", "covS", "covRS"), c("r", "s", "r,s")),
                 function(.n) {
                   if (exists(.n, envir = .ret, inherits = FALSE)) get(.n, envir = .ret) else NULL
                 })
  # covMethod="s"/"r" write only e["cov"] -- the chosen covariance is not always
  # mirrored into covR/covS/covRS -- so cache the installed native under its own type too
  .envType <- .covFdType(.env)
  if (nzchar(.envType) && is.null(.nat[[.envType]]) &&
        exists("cov", envir = .ret, inherits = FALSE)) {
    .nat[[.envType]] <- get("cov", envir = .ret)
  }
  .ret$cov <- .cov
  # Keep the reported covMethod consistent with what was installed: routing on the
  # requested control can install a sandwich where the env still says "s".  Only rewrite
  # the TYPE when it differs, so the env's "r+"/"|r|" decorations survive when they
  # agree; either way the name carries the " (full)" scope suffix.
  .ret$covMethod <- .covFullName(if (identical(.type, .envType)) .env else .type)
  .ret$covR <- .Rinv
  if (!is.null(.covS)) {
    dimnames(.covS) <- dimnames(.Rinv)
    .ret$covS <- .covS
  }
  if (!is.null(.covRS)) {
    dimnames(.covRS) <- dimnames(.Rinv)
    .ret$covRS <- .covRS
  }
  for (.n in names(.nat)) .covCacheAdd(.ret, .n, .nat[[.n]])
  .covCacheAdd(.ret, .covFullName("r"), .Rinv)
  .covCacheAdd(.ret, .covFullName("s"), .covS)
  .covCacheAdd(.ret, .covFullName("r,s"), .covRS)
  .covCacheDrop(.ret, .ret$covMethod)
  .covCacheDrop(.ret, .covFullName(.type))
  .foceiCovCondition(.ret, .cov, .ev)
  # Report the swap: the SEs the C++ step derived from the native theta-only
  # covariance describe a matrix that is no longer $cov, so the caller must
  # refresh the parameter table.
  invisible(TRUE)
}
