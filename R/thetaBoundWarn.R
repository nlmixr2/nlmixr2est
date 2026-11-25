#' Check the `thetaBoundWarn` control option
#' @param thetaBoundWarn relative tolerance; `NA`/`0` is off
#' @return invisibly the option
#' @noRd
.thetaBoundWarnCheck <- function(thetaBoundWarn) {
  checkmate::assertNumber(thetaBoundWarn, lower = 0, upper = 1, na.ok = TRUE, .var.name = "thetaBoundWarn")
}

#' The theta-bound warning tolerance of a fit
#'
#' @param fit the finalized fit
#' @return the relative tolerance below which an estimate counts as at its
#'   bound; `NA` when the warning is disabled
#' @noRd
.thetaBoundWarnTol <- function(fit) {
  .ctl <- tryCatch(fit$control, error = function(e) NULL)
  .v <- if (is.list(.ctl)) .ctl[["thetaBoundWarn"]] else NULL
  if (is.null(.v)) {
    return(1e-4)
  }
  .v <- suppressWarnings(as.double(.v))
  if (length(.v) != 1L || !is.finite(.v) || .v <= 0) {
    return(NA_real_)
  }
  .v
}

#' Estimated thetas that ended at (or within tolerance of) a finite bound
#'
#' The tolerance is relative to the bound range when both bounds are finite,
#' otherwise to the larger of `|bound|` and the initial distance from it (or
#' 1 when the initial estimate is unknown).  The IOV rewrite's own
#' (`nlmixr2iov`) thetas are left out: their bound is not the user's.
#' @param iniDf the fit's `iniDf` (final estimates, bounds on the original scale)
#' @param tol relative tolerance; `NA` names nothing
#' @param ini0 named initial theta estimates, or `NULL`
#' @return the names of the thetas at a bound
#' @noRd
.thetaAtBound <- function(iniDf, tol, ini0 = NULL) {
  .need <- c("ntheta", "name", "est", "lower", "upper", "fix")
  if (is.na(tol) || !is.data.frame(iniDf) || !all(.need %in% names(iniDf))) {
    return(character(0))
  }
  .keep <- !is.na(iniDf$ntheta) & !iniDf$fix & is.finite(iniDf$est)
  if ("backTransform" %in% names(iniDf)) {
    .keep <- .keep & !(!is.na(iniDf$backTransform) & grepl("^nlmixr2iov", iniDf$backTransform))
  }
  .d <- iniDf[.keep, , drop = FALSE]
  if (nrow(.d) == 0L) {
    return(character(0))
  }
  .i0 <- if (is.null(ini0)) rep(NA_real_, nrow(.d)) else unname(ini0[.d$name])
  .near <- function(bound, other) {
    .scale <- ifelse(
      is.finite(other),
      other - bound,
      pmax(abs(bound), ifelse(is.finite(.i0), abs(.i0 - bound), 1))
    )
    is.finite(bound) & abs(.d$est - bound) <= tol * abs(.scale)
  }
  .d$name[.near(.d$lower, .d$upper) | .near(.d$upper, .d$lower)]
}

#' Warn when an estimated theta ended at its bound
#'
#' An estimate stopped by its `lower`/`upper` bound rather than by the
#' likelihood gives unreliable estimates and standard errors (#1505).  Emits
#' one concise warning (collected into `$runInfo`) naming those thetas.
#' @param fit the finalized fit
#' @param env the estimation environment; its `ui` holds the initial estimates
#' @return the names of the thetas at a bound, invisibly
#' @noRd
.thetaBoundWarn <- function(fit, env) {
  .tol <- .thetaBoundWarnTol(fit)
  if (is.na(.tol)) {
    return(invisible(character(0)))
  }
  .ui0 <- if (is.environment(env)) get0("ui", envir = env, inherits = FALSE)
  .ini0 <- tryCatch(
    {
      .d0 <- .ui0$iniDf
      .d0 <- .d0[!is.na(.d0$ntheta), , drop = FALSE]
      setNames(.d0$est, .d0$name)
    },
    error = function(e) NULL
  )
  .at <- .thetaAtBound(tryCatch(fit$ui$iniDf, error = function(e) NULL), .tol, .ini0)
  if (length(.at) > 0L) {
    warning(.runInfoNameMsg("estimate at its bound: ", paste0("\"", .at, "\"")), call. = FALSE)
  }
  invisible(.at)
}
