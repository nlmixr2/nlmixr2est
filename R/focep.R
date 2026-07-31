#' Control options for the focep (FOCE+) estimation method
#'
#' This is the first order conditional estimation without eta/residual
#' interaction, but keeping the live conditional residual variance R (the
#' `foce = "foce+"` option of [foceiControl()]).  It is the `foce` method with
#' `foce = "foce+"` forced; use `foce` (est = "foce") for the NONMEM-matching
#' frozen-R behavior.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `focei` instead
#' @param foce FOCE residual-variance mode; for `focepControl()` this is
#'   always `"foce+"` and cannot be changed -- use `foceControl()` for
#'   `"nonmem"`
#' @return focepControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' focepControl()
focepControl <- function(sigdig=3,
                         ...,
                         interaction=FALSE,
                         foce="foce+") {
  .control <- foceiControl(sigdig=sigdig, ...,
                           interaction=FALSE, foce="foce+")
  class(.control) <- "focepControl"
  .control
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.focepControl <- function(control, env) {
  assign("focepControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.focep <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- focepControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("focepControl", .ctl)
  if (inherits(.ctl, "foceiControl") ||
        inherits(.ctl, "foControl") ||
        inherits(.ctl, "foiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to focepControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(focepControl, .ctl)
  } else if (!inherits(.ctl, "focepControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- focepControl()
  } else {
    .ctl <- do.call(focepControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.focep <- function(x, ...) {
  .env <- x[[1]]
  if (exists("focepControl", .env)) {
    .control <- get("focepControl", .env)
    if (inherits(.control, "focepControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "focepControl")) return(.control)
  }
  stop("cannot find focep related control object", call.=FALSE)
}

.focepControlToFoceiControl <- function(env, assign=TRUE) {
  .focepControl <- env$focepControl
  .ui <- env$ui
  .n <- names(.focepControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.focepControl$interaction)
                                     }
                                     .focepControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.focep <- function(x, ...) {
  .env <- x[[1]]
  .focepControlToFoceiControl(.env, assign=FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.focep <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'focep'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="focepControl")
  .focepControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$focepControl <- .control
  env$est <- "focep"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="focep")
}
attr(nlmixr2Est.focep, "iov") <- TRUE
attr(nlmixr2Est.focep, "covPresent") <- TRUE
attr(nlmixr2Est.focep, "unbounded") <- .foUnbounded
