.preProcessHooks <- new.env(parent=emptyenv())

#' This adds a pre-processing hook to nlmixr2est
#'
#' This pre-processing hook is run before the estimation process begins.  It is
#' useful for modifying the user interface, the estimation object, the data, or
#' the control object before the estimation process begins.  The function must
#' take four arguments: ui, est, data, and control.  The function must return a
#' list with elements 'ui', 'est', 'data', and/or 'control'.  If the element is
#' not returned, the original object is used.  If the element is returned, the
#' original object is replaced with the new object.
#'
#' @param name Character vector representing the name of the hook
#' @param fun The function to run
#' @return The function that was added (invisibly)
#' @export
#' @author Matthew L. Fidler
#' @family preProcessHooks
#' @keywords internal
preProcessHooksAdd <- function(name, fun) {
  checkmate::assertCharacter(name, len=1,
                             pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                             min.chars=1)
  checkmate::assertFunction(fun, args=c("ui", "est", "data", "control"))
  if (exists(name, envir=.preProcessHooks)) {
    stop("nlmixr2est preprocessing hook '", name, "' already exists",
         call.=FALSE)
  }
  assign(name, fun, envir=.preProcessHooks)
  invisible(fun)
}
#' Remove the hook from nlmixr2est
#'
#' This removes the hook from nlmixr2est.  If the hook does not exist, a warning
#' is issued.
#'
#' @param name Character vector representing the name of the hook
#' @return logical indicating if the hook was removed (invisibly)
#' @export
#' @author Matthew L. Fidler
#' @family preProcessHooks
#' @keywords internal
preProcessHooksRm <- function(name) {
  checkmate::assertCharacter(name, len=1,
                             pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                             min.chars=1)
  if (exists(name, envir=.preProcessHooks)) {
    rm(list=name, envir=.preProcessHooks)
    invisible(TRUE)
  } else {
    warning("nlmixr2est preprocessing hook '", name, "' does not exist",
            call.=FALSE)
    invisible(FALSE)
  }
}
#' Return a list of all pre-processing hooks
#'
#'
#' @param name when specified, the name of the hook, otherwise `NULL`
#'   to list all hooks
#'
#' @return a charcter vector listing all pre-processing hooks or the
#'   function for the hook
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
#' @family preProcessHooks
preProcessHooks <- function(name=NULL) {
  if (is.null(name)) {
    ls(envir=.preProcessHooks, all.names=TRUE)
  } else {
    checkmate::assertCharacter(name, len=1,
                               pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                               min.chars=1)
    if (exists(name, envir=.preProcessHooks)) {
      get(name, envir=.preProcessHooks)
    } else {
      stop("nlmixr2est preprocessing hook '", name, "' does not exist",
           call.=FALSE)
    }
  }
}

.preProcessHooksRun <- function(env, est) {
  .ui <- env$ui
  .est <- est
  if (is.null(.est) && is.character(nlmixr2global$nlmixr2pipeEst)) {
    .est <- est <- nlmixr2global$nlmixr2pipeEst
  }
  .data <- env$data
  if (is.null(.data) &&
        inherits(nlmixr2global$nlmixr2pipeData, "data.frame")) {
    .data <- env$data <- nlmixr2global$nlmixr2pipeData
  }
  .control <- env$control
  if (is.null(.control) && inherits(nlmixr2global$nlmixr2pipeControl, "list")) {
    .control <- env$control <- nlmixr2global$nlmixr2pipeControl
  }
  for (name in preProcessHooks()) {
    .fun <- get(name, envir=.preProcessHooks)
    .ret <- .fun(.ui, .est, .data, .control)
    if (is.null(.ret) || length(.ret) == 0) {
    } else if (checkmate::testList(.ret, max.len=4, min.len=1) &&
          !is.null(names(.ret)) &&
          checkmate::checkSubset(names(.ret), c("ui", "est", "data", "control"))) {
      if (!is.null(.ret$ui)) {
        .ui <- .ret$ui
      }
      if (!is.null(.ret$est)) {
        .est <- .ret$est
      }
      if (!is.null(.ret$data)) {
        .data <- .ret$data
      }
      if (!is.null(.ret$control)) {
        .control <- .ret$control
      }
    } else {
      stop("nlmixr2est preprocessing hook '", name, "' must return a list with elements 'ui', 'est', 'data', and/or 'control'",
             call.=FALSE)
    }
  }
  env$ui <- .ui
  env$data <- .data
  env$control <- .control
  .est
}
