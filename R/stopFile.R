#' Start polling the control's stop file
#'
#' Only the outermost estimation polls; a nested one keeps its file.
#' @param control the estimation's control
#' @return `TRUE` when this call started polling
#' @noRd
.nmStopFileStart <- function(control) {
  if (!is.list(control)) {
    return(FALSE)
  }
  .file <- control$stopFile
  if (!checkmate::testString(.file, min.chars = 1L) || .nmStopFileGet() != "") {
    return(FALSE)
  }
  .nmStopFileSet(normalizePath(.file, mustWork = FALSE))
  TRUE
}

#' Stop polling the stop file, noting a stop request in `$runInfo`
#' @return `TRUE` when a stop was requested
#' @noRd
.nmStopFileEnd <- function() {
  .stopped <- .nmStopFileRequested()
  .nmStopFileSet("")
  if (.stopped) {
    warning("stopped early by user request", call. = FALSE)
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
