#' Check the `corWarn` control option
#' @param corWarn absolute correlation that triggers the warning; `NA`/`Inf` is off
#' @return invisibly `TRUE`; an invalid value is an error
#' @noRd
.corWarnCheck <- function(corWarn) {
  checkmate::assertNumber(corWarn, lower = 0, na.ok = TRUE, .var.name = "corWarn")
  invisible(TRUE)
}

# the most pairs a correlation warning lists
.corWarnMaxPairs <- 3L

#' Pairs of estimates whose correlation reaches a threshold
#'
#' Only the upper triangle is scanned; rows with a missing or non-positive
#' variance and the names in `skip` are left out.
#' @param cov covariance matrix with dimnames
#' @param corWarn absolute correlation threshold; `NA`/`Inf` gives none
#' @param skip names to leave out (fixed parameters)
#' @return data frame with `p1`, `p2` and `r`, strongest first
#' @noRd
.covCorPairs <- function(cov, corWarn, skip = character(0)) {
  .none <- data.frame(p1 = character(0), p2 = character(0), r = numeric(0))
  if (!is.matrix(cov) || nrow(cov) < 2L || !is.finite(corWarn)) {
    return(.none)
  }
  .n <- rownames(cov)
  if (is.null(.n)) {
    return(.none)
  }
  .d <- diag(cov)
  .keep <- which(is.finite(.d) & .d > 0 & !(.n %in% skip))
  if (length(.keep) < 2L) {
    return(.none)
  }
  .c <- cov[.keep, .keep, drop = FALSE]
  .r <- .c / sqrt(outer(diag(.c), diag(.c)))
  .w <- which(upper.tri(.r) & is.finite(.r) & abs(.r) >= corWarn, arr.ind = TRUE)
  if (nrow(.w) == 0L) {
    return(.none)
  }
  .ret <- data.frame(
    p1 = rownames(.c)[.w[, 1]],
    p2 = rownames(.c)[.w[, 2]],
    r = .r[.w]
  )
  .ret[order(-abs(.ret$r)), , drop = FALSE]
}

#' Warning text naming highly correlated estimates
#' @param pairs result of `.covCorPairs()`
#' @param max most pairs to list
#' @return the message, or `NULL` when there are no pairs
#' @noRd
.covCorMsg <- function(pairs, max = .corWarnMaxPairs) {
  if (nrow(pairs) == 0L) {
    return(NULL)
  }
  .p <- utils::head(pairs, max)
  .txt <- paste0(.p$p1, "/", .p$p2, " (", sprintf("%.2f", .p$r), ")", collapse = ", ")
  if (nrow(pairs) > max) {
    .txt <- paste0(.txt, ", ...")
  }
  paste0("highly correlated estimates: ", .txt)
}

#' The `corWarn` threshold a fit was estimated with
#'
#' The fit's own control, then its `foceiControl`; `0.95` when neither has one.
#' @param fit the fit
#' @param env the fit environment
#' @return single number
#' @noRd
.covCorWarnOf <- function(fit, env) {
  .ctl <- tryCatch(fit$control, error = function(e) NULL)
  .v <- if (is.list(.ctl)) .ctl$corWarn else NULL
  if (is.null(.v)) {
    .v <- tryCatch(env$foceiControl$corWarn, error = function(e) NULL)
  }
  if (is.null(.v)) 0.95 else .v
}

#' Warn when the fit's covariance shows highly correlated estimates
#' @param fit the finalized fit
#' @return invisibly the (clipped) warning text, or `NULL` when none was given
#' @noRd
.covCorWarn <- function(fit) {
  .env <- if (is.environment(fit)) fit else tryCatch(fit$env, error = function(e) NULL)
  if (!is.environment(.env) || !exists("cov", envir = .env, inherits = FALSE)) {
    return(invisible(NULL))
  }
  .ini <- tryCatch(fit$iniDf, error = function(e) NULL)
  # a full covariance names an omega "om.<eta>"
  .skip <- if (is.data.frame(.ini)) .ini$name[.ini$fix] else character(0)
  .skip <- c(.skip, paste0("om.", .skip))
  .pairs <- .covCorPairs(get("cov", envir = .env, inherits = FALSE), .covCorWarnOf(fit, .env), .skip)
  .msg <- .covCorMsg(.pairs)
  if (is.null(.msg)) {
    return(invisible(NULL))
  }
  .covWarn(.msg)
}
