#' Control options for the foce estimation method
#'
#' This is the first order option without the interaction between
#' residuals and etas.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `focei` instead
#' @return foceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' foceControl()
foceControl <- function(sigdig=3,
                        ...,
                        interaction=FALSE) {
  .control <- foceiControl(sigdig=sigdig, ...,
                           interaction=FALSE)
  class(.control) <- "foceControl"
  .control
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.foceControl <- function(control, env) {
  ## eval(rxode2::rxUiDeparse(control, "control"))
  assign("foceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.foce <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- foceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("foceControl", .ctl)
  if (inherits(.ctl, "foceiControl") ||
        inherits(.ctl, "foControl") ||
        inherits(.ctl, "foiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to foceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(foceControl, .ctl)
  } else if (!inherits(.ctl, "foceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- foceControl()
  } else {
    .ctl <- do.call(foceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.foce <- function(x, ...) {
  .env <- x[[1]]
  if (exists("foceControl", .env)) {
    .control <- get("foceControl", .env)
    if (inherits(.control, "foceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "foceControl")) return(.control)
  }
  stop("cannot find foce related control object", call.=FALSE)
}

.foceControlToFoceiControl <- function(env, assign=TRUE) {
  .foceControl <- env$foceControl
  .ui <- env$ui
  .n <- names(.foceControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.foceControl$interaction)
                                     }
                                     .foceControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.foce <- function(x, ...) {
  .env <- x[[1]]
  .foceControlToFoceiControl(.env, assign=FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.foce <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'foce'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .uiApplyIov(env)
  .foceiFamilyControl(env, ..., type="foceControl")
  .foceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$foceControl <- .control
  env$est <- "foce"
  .ui <- env$ui
  .uiFinalizeIov(.foceiFamilyReturn(env, .ui, ..., est="foce"))
}
attr(nlmixr2Est.foce, "covPresent") <- TRUE
