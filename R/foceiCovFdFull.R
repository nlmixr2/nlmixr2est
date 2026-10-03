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
  if (!is.null(ui$boundedTransforms) && length(ui$boundedTransforms) > 0L) {
    return(NULL)
  }
  ini <- ui$iniDf
  isFix <- if (is.null(ini$fix)) rep(FALSE, nrow(ini)) else ini$fix
  isFix[is.na(isFix)] <- FALSE
  thRows <- which(!is.na(ini$ntheta) & !isFix) # structural + residual thetas
  if (length(thRows) == 0L) {
    return(NULL)
  }
  thPos <- as.integer(ini$ntheta[thRows] - 1L) # 0-based position in fullTheta
  thNames <- ini$name[thRows]
  Om <- get("omega", e)
  pairs <- .foceiOmegaPairs(Om, ini) # free Omega lower-triangle (a>=b)
  if (is.null(pairs) || nrow(pairs) == 0L) {
    return(list(thPos = thPos, omA = integer(0), omB = integer(0), names = thNames))
  }
  map <- .foceiEtaThetaMap(ui)
  onm <- map$etaNames # Omega named by the eta
  omNames <- .foceiOmegaCovNames(pairs, onm)
  list(thPos = thPos, omA = as.integer(pairs[, 1]), omB = as.integer(pairs[, 2]), names = c(thNames, omNames))
}

#' The three FD-full covariances, each checked by `.covGuard()`
#'
#' The sandwich is usable only when `Rinv` and `solve(S)` both are: built from
#' an indefinite `Rinv` it is still positive semi-definite, so it can pass the
#' check on its own.
#' @param Rinv full Hessian inverse (`.fdFullCov`)
#' @param S full score cross-product (`.fdFullS`), or `NULL`
#' @return named list ("r", "s", "r,s") of `.covGuard()` results, the matrices
#'   dimnamed like `Rinv`
#' @noRd
.foceiFdFullShapes <- function(Rinv, S) {
  .dn <- dimnames(Rinv)
  .shapes <- list(
    r = Rinv,
    s = if (is.matrix(S)) tryCatch(solve(S), error = function(e) NULL),
    "r,s" = if (is.matrix(S)) Rinv %*% S %*% Rinv
  )
  for (.n in names(.shapes)) {
    if (is.matrix(.shapes[[.n]])) {
      dimnames(.shapes[[.n]]) <- .dn
    }
    .shapes[[.n]] <- .covGuard(.shapes[[.n]])
  }
  if (!.shapes$r$ok || !.shapes$s$ok) {
    .shapes[["r,s"]] <- list(ok = FALSE, reason = paste0("needs a positive-definite ", if (.shapes$r$ok) "S" else "R"))
  }
  .shapes
}

#' Install the C++ FD-full covariance as `fit$cov` (and `fit$covR/covS/covRS`) when
#' `covFull = TRUE`, routing on the requested covMethod: "r,s" -> the sandwich
#' `Rinv %*% S %*% Rinv`, "s" -> `solve(S)`, "r" -> `Rinv`.  The native cov is
#' kept -- with a warning when the requested shape is not usable, silently when
#' covMethod is not an FD method or the pieces were not computed.  Every usable
#' shape, native or full, is cached for `setCov()`; an unusable one is neither
#' stored nor cached.  FD counterpart to [.foceiInstallAnalyticCov].
#' @param .ret focei fit environment
#' @return invisibly TRUE when the full covariance was installed
#' @noRd
.foceiInstallFdFullCov <- function(.ret) {
  if (!exists(".fdFullCov", envir = .ret, inherits = FALSE)) {
    return(invisible(FALSE))
  }
  # The fit env's covMethod records what the NATIVE theta-only step produced, not what was
  # asked for -- C++ downgrades it to "s" when that step's "r" fails.  The full R/S pieces
  # used below are computed independently of that step, so route on the REQUESTED control;
  # otherwise a requested "r,s" silently installs solve(S) with a usable .Rinv in hand.
  .env <- if (.covIsName(.ret$covMethod)) .ret$covMethod else ""
  .cty <- tryCatch(rxode2::rxGetControl(.ret$ui, "covType", "fd"), error = function(e) "fd")
  .req <- tryCatch(rxode2::rxGetControl(.ret$ui, "covMethod", NA_integer_), error = function(e) NA_integer_)
  .req <- if (identical(.cty, "analytic")) "" else .covMethodFromSlot(.req)
  .type <- .covFdType(if (nzchar(.req)) .req else .env)
  .S <- get0(".fdFullS", envir = .ret, inherits = FALSE)
  if (!nzchar(.type) || (.type != "r" && is.null(.S))) {
    return(invisible(FALSE))
  } # analytic / failed / "" / boundary, or no S computed -> keep native
  .full <- .foceiFdFullShapes(get(".fdFullCov", envir = .ret), .S)
  # the native theta-only pieces, cached so setCov() can swap to that shape without
  # recomputing anything (they are already in hand)
  .nat <- stats::setNames(mget(c("covR", "covS", "covRS"), envir = .ret, ifnotfound = list(NULL)), c("r", "s", "r,s"))
  .installed <- .full[[.type]]$ok
  if (.installed) {
    # covMethod="s"/"r" write only e["cov"] -- the chosen covariance is not always
    # mirrored into covR/covS/covRS -- so cache the native about to be replaced too
    .envType <- .covFdType(.env)
    if (nzchar(.envType) && is.null(.nat[[.envType]])) {
      .nat[[.envType]] <- .ret$cov
    }
    # Keep the reported covMethod consistent with what was installed: routing on the
    # requested control can install a sandwich where the env still says "s".  Only
    # rewrite the TYPE when it differs, so the env's "r+"/"|r|" decorations survive when
    # they agree; either way the name carries the " (full)" scope suffix.
    .covInstall(
      .ret,
      .full[[.type]]$cov,
      .covFullName(if (identical(.type, .envType)) .env else .type),
      stash = FALSE,
      refresh = "none"
    )
    for (.n in names(.full)) {
      if (.full[[.n]]$ok) assign(c(r = "covR", s = "covS", "r,s" = "covRS")[[.n]], .full[[.n]]$cov, envir = .ret)
    }
  } else {
    .covRejectWarn(.ret, .covFullName(.type), .full[[.type]]$reason)
  }
  for (.n in names(.nat)) {
    if (.covGuard(.nat[[.n]])$ok) .covCacheAdd(.ret, .n, .nat[[.n]])
  }
  for (.n in names(.full)) {
    if (.full[[.n]]$ok) .covCacheAdd(.ret, .covFullName(.n), .full[[.n]]$cov)
  }
  .covCacheDrop(.ret, .ret$covMethod)
  .covCacheDrop(.ret, .covFullName(.type))
  # Report the swap: the SEs the C++ step derived from the native theta-only
  # covariance describe a matrix that is no longer $cov, so the caller must
  # refresh the parameter table.
  invisible(.installed)
}
