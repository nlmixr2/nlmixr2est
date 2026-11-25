#' Start polling the control's stop file and time limit
#'
#' Only the outermost estimation polls; a nested one keeps its file and
#' time limit.
#' @param control the estimation's control
#' @return `TRUE` when this call started polling
#' @noRd
.nmStopFileStart <- function(control) {
  if (!is.list(control) || .nmStopActive()) {
    return(FALSE)
  }
  .file <- control$stopFile
  .file <- if (checkmate::testString(.file, min.chars = 1L)) {
    normalizePath(.file, mustWork = FALSE)
  } else {
    ""
  }
  .maxTime <- control$maxTime
  .maxTime <- if (checkmate::testNumber(.maxTime, lower = 0, finite = TRUE)) {
    as.double(.maxTime)
  } else {
    NA_real_
  }
  if (.file == "" && is.na(.maxTime)) {
    return(FALSE)
  }
  .nmStopFileSet(.file)
  .nmStopTimeSet(.maxTime)
  TRUE
}

#' Stop polling the stop file and time limit, noting a stop in `$runInfo`
#' @return `TRUE` when a stop was requested
#' @noRd
.nmStopFileEnd <- function() {
  .stopped <- .nmStopFileStopped()
  .msg <- .nmStopMessage()
  .nmStopFileSet("")
  .nmStopTimeSet(NA_real_)
  if (.stopped) {
    warning(.msg, call. = FALSE)
  }
  invisible(.stopped)
}

#' Run an R outer optimizer, ending it at a stop request
#'
#' @param opt the optimizer, called with `...`
#' @param ... passed to `opt`
#' @return `opt`'s result, or `NULL` when a stop request ended it
#' @noRd
.foceiOuterOpt <- function(opt, ...) {
  tryCatch(opt(...), error = function(e) {
    if (.nmStopFileRequested()) {
      return(NULL)
    }
    stop(e)
  })
}

#' Run an nlm-family optimizer, ending it at a stop request
#'
#' A stop request ends the optimizer with an error (or a failed
#' evaluation); its result is then the lowest evaluation, marked
#' `nlmStopped`.
#' @param expr the optimizer call
#' @param par name of the estimate in its result
#' @param objective name of the objective in its result
#' @return the optimizer's result
#' @noRd
.nlmStopOpt <- function(expr, par = "par", objective = "fval") {
  .nlmStopOn(TRUE)
  on.exit(.nlmStopOn(FALSE))
  .ret <- tryCatch(expr, error = function(e) {
    if (.nmStopFileRequested() && !is.null(.nlmStopBest())) {
      return(list())
    }
    stop(e)
  })
  if (!.nmStopFileRequested()) {
    return(.ret)
  }
  .nmStopFileConsume()
  .best <- .nlmStopBest()
  if (!is.null(.best)) {
    .ret[[par]] <- .best$par
    .ret[[objective]] <- .best$f
  }
  .ret$hessian <- NULL
  .ret$nlmStopped <- TRUE
  .ret
}
