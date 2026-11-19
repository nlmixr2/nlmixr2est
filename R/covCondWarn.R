#' Check the `condWarn` control option
#' @param condWarn correlation condition number that triggers the warning;
#'   `NA`/`Inf` is off
#' @return invisibly `TRUE`; an invalid value is an error
#' @noRd
.condWarnCheck <- function(condWarn) {
  checkmate::assertNumber(condWarn, lower = 1, na.ok = TRUE, .var.name = "condWarn")
  invisible(TRUE)
}

#' Condition number of a covariance's correlation matrix
#'
#' Rows with a missing or non-positive variance are dropped first.
#' @param cov covariance matrix
#' @return single number, `NA` when it cannot be computed
#' @noRd
.covCondCor <- function(cov) {
  if (!is.matrix(cov) || nrow(cov) == 0L) {
    return(NA_real_)
  }
  .d <- diag(cov)
  .keep <- which(is.finite(.d) & .d > 0)
  if (length(.keep) == 0L) {
    return(NA_real_)
  }
  .c <- cov[.keep, .keep, drop = FALSE]
  .e <- try(eigen(stats::cov2cor(.c), symmetric = TRUE, only.values = TRUE), silent = TRUE)
  if (inherits(.e, "try-error") || any(!is.finite(.e$values))) {
    return(NA_real_)
  }
  .a <- abs(.e$values)
  max(.a) / min(.a)
}

#' Warning text for a large condition number
#' @param cond the correlation-matrix condition number
#' @param condWarn the threshold; `NA`/`Inf` gives none
#' @return the message, or `NULL` when there is nothing to warn about
#' @noRd
.covCondMsg <- function(cond, condWarn) {
  if (length(cond) != 1L || is.na(cond) || length(condWarn) != 1L || !is.finite(condWarn) || cond <= condWarn) {
    return(NULL)
  }
  sprintf(
    "covariance condition number is large (%s > %s)",
    format(cond, digits = 5L),
    format(condWarn, digits = 5L)
  )
}

#' The `condWarn` threshold a fit was estimated with
#'
#' The fit's own control, then its `foceiControl`; `1000` when neither has one.
#' @param fit the fit
#' @param env the fit environment
#' @return single number
#' @noRd
.covCondWarnOf <- function(fit, env) {
  .ctl <- tryCatch(fit$control, error = function(e) NULL)
  .v <- if (is.list(.ctl)) .ctl$condWarn else NULL
  if (is.null(.v)) {
    .v <- tryCatch(env$foceiControl$condWarn, error = function(e) NULL)
  }
  if (is.null(.v)) 1000 else .v
}

#' Warn when the fit's covariance is ill-conditioned
#'
#' Uses the `conditionNumberCor` the fit reports, or computes it from the
#' installed covariance when the fit has none.
#' @param fit the finalized fit (or its environment)
#' @return invisibly the warning text, or `NULL` when none was given
#' @noRd
.covCondWarn <- function(fit) {
  .env <- if (is.environment(fit)) fit else tryCatch(fit$env, error = function(e) NULL)
  if (!is.environment(.env) || !exists("cov", envir = .env, inherits = FALSE)) {
    return(invisible(NULL))
  }
  .cw <- .covCondWarnOf(fit, .env)
  if (length(.cw) != 1L || !is.finite(.cw)) {
    return(invisible(NULL))
  }
  .cond <- if (exists("conditionNumberCor", envir = .env, inherits = FALSE)) {
    suppressWarnings(as.double(get("conditionNumberCor", envir = .env, inherits = FALSE)))
  } else {
    NA_real_
  }
  if (length(.cond) != 1L || !is.finite(.cond)) {
    .cond <- .covCondCor(get("cov", envir = .env, inherits = FALSE))
  }
  .msg <- .covCondMsg(.cond, .cw)
  if (is.null(.msg)) {
    return(invisible(NULL))
  }
  .covWarn(.msg)
}
