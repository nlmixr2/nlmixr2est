.postFinalObjectHooks <- new.env(parent=emptyenv())

#' Register a post-final-object hook run after estimation completes
#'
#' `fun` must take one argument (`ret`) and return the finalized return
#' object (or `NULL`/nothing to leave it unchanged).
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



#' Remove a post-final-object hook from nlmixr2est
#'
#' Warns if the hook does not exist.
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

#' List all post-final-object processing hooks
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
  # Runs registered post-final-object hooks to update the return object
  # after estimation completes, before it is returned to the user.
  .ret <- ret
  for (name in postFinalObjectHooks()) {
    .fun <- get(name, envir=.postFinalObjectHooks)
    .tmp <- .fun(.ret)
    if (is.null(.tmp) || length(.tmp) == 0) {
    } else {
      .ret <- .tmp
    }
  }
  .ret
}



.preFinalParTableHooks <- new.env(parent=emptyenv())
#' Register a pre-final parameter table hook run before final table assembly
#'
#' `fun` must take one argument (`env`) and update it in place (e.g. `cov`,
#' `theta`, `thetaNames`, `thetaDf`) before the final tables are built.
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
  if (exists(name, envir=.preFinalParTableHooks)) {
    stop("nlmixr2est pre-table hook '", name, "' already exists",
         call.=FALSE)
  }
  assign(name, fun, envir=.preFinalParTableHooks)
  invisible(fun)
}

#' Remove a pre-final parameter table hook from nlmixr2est
#'
#' Warns if the hook does not exist.
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

#' List all pre-final parameter table processing hooks
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
  # Updates cov/theta/thetaNames/thetaDf/etaObf/ui before final table
  # assembly; hooks should check objects exist before modifying them.
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

#' Register a pre-processing hook run before estimation begins
#'
#' `fun` must take four arguments (`ui`, `est`, `data`, `control`) and return
#' a list with any of `'ui'`, `'est'`, `'data'`, `'control'` to override; a
#' non-returned element is left unchanged.
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
#' Remove a pre-processing hook from nlmixr2est
#'
#' Warns if the hook does not exist.
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
#' List all pre-processing hooks
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
  # Keep per-call copies of the original-model restore info; nested nlmixr2()
  # calls during estimation (setOfv/addCwres/...) reset the globals before the
  # post-estimation restore in .nlmixrEstUpdatesOrigModel() runs (issue #741)
  env$nlmixrPureInputUi <- nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi
  env$uiUnfix <- nlmixr2global$nlmixr2EstEnv$uiUnfix
  .ret[[1]]
}
