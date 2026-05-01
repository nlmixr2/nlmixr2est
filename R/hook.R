.postFinalObjectHooks <- new.env(parent=emptyenv())

#' This adds a pre-final table processing hook to nlmixr2est
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
postFinalObjectHooksAdd <- function(name, fun) {
  checkmate::assertCharacter(name, len=1,
                             pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                             min.chars=1)
  checkmate::assertFunction(fun, args="ret")
  if (exists(name, envir=.preProcessHooks)) {
    stop("nlmixr2est pre-table hook '", name, "' already exists",
         call.=FALSE)
  }
  assign(name, fun, envir=.postFinalObjectHooks)
  invisible(fun)
}



#' Remove the hook from the pre-final table nlmixr2est
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
postFinalObjectHooksRm <- function(name) {
  checkmate::assertCharacter(name, len=1,
                             pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                             min.chars=1)
  if (exists(name, envir=.postFinalObjectHooks)) {
    rm(list=name, envir=.postFinalObjectHooks)
    invisible(TRUE)
  } else {
    warning("nlmixr2est pre-table hook '", name, "' does not exist",
            call.=FALSE)
    invisible(FALSE)
  }
}

#' Return a list of all pre-final table processing hooks
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
postFinalObjectHooks <- function(name=NULL) {
  if (is.null(name)) {
    ls(envir=.postFinalObjectHooks, all.names=TRUE)
  } else {
    checkmate::assertCharacter(name, len=1,
                               pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                               min.chars=1)
    if (exists(name, envir=.postFinalObjectHooks)) {
      get(name, envir=.postFinalObjectHooks)
    } else {
      stop("nlmixr2est pre-final table hook '", name, "' does not exist",
           call.=FALSE)
    }
  }
}

.postFinalObjectHooksRun <- function(ret) {
  # This is for updating the return object from the estimation process.
  # This is used for modifying the return object after the estimation process is complete, but
  # before it is returned to the user.
  .ret <- ret
  for (name in postFinalObjectHooks()) {
    .fun <- get(name, envir=.postFinalObjectHooks)
    .ret <- .fun(.ret)
  }
  .ret
}



.preFinalParTableHooks <- new.env(parent=emptyenv())
#' This adds a pre-final table processing hook to nlmixr2est
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
preFinalParTableHooksAdd <- function(name, fun) {
  checkmate::assertCharacter(name, len=1,
                             pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                             min.chars=1)
  checkmate::assertFunction(fun, args="env")
  if (exists(name, envir=.preProcessHooks)) {
    stop("nlmixr2est pre-table hook '", name, "' already exists",
         call.=FALSE)
  }
  assign(name, fun, envir=.preFinalParTableHooks)
  invisible(fun)
}


#' This adds a pre-final table processing hook to nlmixr2est
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
preFinalParTableHooksAdd <- function(name, fun) {
  checkmate::assertCharacter(name, len=1,
                             pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                             min.chars=1)
  checkmate::assertFunction(fun, args="env")
  if (exists(name, envir=.preProcessHooks)) {
    stop("nlmixr2est pre-table hook '", name, "' already exists",
         call.=FALSE)
  }
  assign(name, fun, envir=.preFinalParTableHooks)
  invisible(fun)
}

#' Remove the hook from the pre-final table nlmixr2est
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
preFinalParTableHooksRm <- function(name) {
  checkmate::assertCharacter(name, len=1,
                             pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                             min.chars=1)
  if (exists(name, envir=.preFinalParTableHooks)) {
    rm(list=name, envir=.preFinalParTableHooks)
    invisible(TRUE)
  } else {
    warning("nlmixr2est pre-table hook '", name, "' does not exist",
            call.=FALSE)
    invisible(FALSE)
  }
}

#' Return a list of all pre-final table processing hooks
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
preFinalParTableHooks <- function(name=NULL) {
  if (is.null(name)) {
    ls(envir=.preFinalParTableHooks, all.names=TRUE)
  } else {
    checkmate::assertCharacter(name, len=1,
                               pattern="^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                               min.chars=1)
    if (exists(name, envir=.preFinalParTableHooks)) {
      get(name, envir=.preFinalParTableHooks)
    } else {
      stop("nlmixr2est pre-final table hook '", name, "' does not exist",
           call.=FALSE)
    }
  }
}

.preFinalParTableHooksRun <- function(env) {
  # This is for updating:
  # - `cov` before running
  # - `theta` before running
  # - `thetaNames` before running
  # - `thetaDf` before running
  # - Possibly `etaObf`  before running
  # - Possibly updating `ui` before running
  # These should check to make sure that the objects exist before
  # modifying them.
  for (name in preFinalParTableHooks()) {
    .fun <- get(name, envir=.preFinalParTableHooks)
    .ret <- .fun(env)
  }
  invisible(NULL)
}

.preProcessHooks <- new.env(parent=emptyenv())

.orderPreProcessHookNames <- function(names) {
  .last <- ".preProcessBoundedTransform"
  if (.last %in% names) {
    c(names[names != .last], .last)
  } else {
    names
  }
}

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
    .orderPreProcessHookNames(ls(envir=.preProcessHooks, all.names=TRUE))
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
  nlmixr2global$preProcessHookWarnings <- character(0)
  .ret <- .collectWarn({
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
  }, lst=TRUE)
  nlmixr2global$preProcessHookWarnings <- .ret[[2]]
  .ret[[1]]
}
