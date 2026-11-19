#' The ETA shrinkage warning threshold of a fit
#'
#' @param fit the finalized fit
#' @return the threshold (percent); `NA` when the warning is disabled
#' @noRd
.shrinkEtaWarnLimit <- function(fit) {
  .ctl <- tryCatch(fit$control, error = function(e) NULL)
  .v <- if (is.list(.ctl)) .ctl[["shrinkEtaWarn"]] else NULL
  if (is.null(.v)) {
    return(30)
  }
  .v <- suppressWarnings(as.double(.v))
  if (length(.v) != 1L || !is.finite(.v)) {
    return(NA_real_)
  }
  .v
}

#' Etas whose SD shrinkage is at or above a threshold
#'
#' @param shrink a fit's shrinkage table, or a list of them (IOV fits)
#' @param limit the threshold in percent; `NA`/`Inf` returns nothing
#' @return character vector of the offending eta names
#' @noRd
.shrinkEtaHigh <- function(shrink, limit) {
  if (length(limit) != 1L || !is.finite(limit)) {
    return(character(0))
  }
  if (is.list(shrink) && !is.data.frame(shrink)) {
    return(unique(unlist(lapply(shrink, .shrinkEtaHigh, limit = limit))))
  }
  if (!is.data.frame(shrink) || !("sd shrinkage (%)" %in% row.names(shrink))) {
    return(character(0))
  }
  # the last column is the residual (IWRES) shrinkage
  .nm <- setdiff(names(shrink), "IWRES")
  .v <- vapply(
    .nm,
    function(n) {
      suppressWarnings(as.double(shrink["sd shrinkage (%)", n]))
    },
    double(1)
  )
  .nm[is.finite(.v) & .v >= limit]
}

#' Warn when a fit's ETA shrinkage is high
#'
#' Emits one concise warning (collected into `$runInfo`) naming the etas whose
#' SD shrinkage is at or above the fit's `shrinkEtaWarn` control.
#' @param fit the finalized fit
#' @return the offending eta names, invisibly
#' @noRd
.shrinkEtaWarn <- function(fit) {
  .limit <- .shrinkEtaWarnLimit(fit)
  if (is.na(.limit)) {
    return(invisible(character(0)))
  }
  .shrink <- tryCatch(fit$shrink, error = function(e) NULL)
  if (is.null(.shrink)) {
    return(invisible(character(0)))
  }
  .high <- .shrinkEtaHigh(.shrink, .limit)
  if (length(.high) > 0L) {
    warning(.shrinkEtaWarnMsg(.high, .limit), call. = FALSE)
  }
  invisible(.high)
}

#' The high-shrinkage warning, kept to one `$runInfo` line
#'
#' @param etas the offending eta names
#' @param limit the threshold in percent
#' @return the message
#' @noRd
.shrinkEtaWarnMsg <- function(etas, limit) {
  .pre <- sprintf("high eta shrinkage (>=%s%%): ", format(limit))
  .n <- length(etas)
  for (.k in rev(seq_len(.n))) {
    .msg <- paste0(.pre, paste(etas[seq_len(.k)], collapse = ", "))
    if (.k < .n) {
      .msg <- paste0(.msg, sprintf(", +%d more", .n - .k))
    }
    if (nchar(.msg) < 75L || .k == 1L) {
      return(.msg)
    }
  }
}
