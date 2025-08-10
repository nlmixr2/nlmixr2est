#' Control options for the posthoc estimation method
#'
#' This option is for simply getting the maximum a-prior (MAP) also
#' called the posthoc estimates
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiConrol()`
#' @param maxOuterIterations ignored, posthoc always sets this to 0.
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`, though you can set it to be `TRUE` as well.
#' @return posthocControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' posthocControl()
posthocControl <- function(sigdig=3,
                           ...,
                           interaction=FALSE,
                           maxOuterIterations=NULL) {
  .control <- foceiControl(sigdig=sigdig, ...,
                           maxOuterIterations=0L,
                           interaction=interaction)
  class(.control) <- "posthocControl"
  .control
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.posthocControl <- function(control, env) {
  assign("posthocControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.posthoc <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- posthocControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("posthocControl", .ctl)
  if (!inherits(.ctl, "posthocControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- posthocControl()
  } else {
    .ctl <- do.call(posthocControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.posthoc <- function(x, ...) {
  .env <- x[[1]]
  if (exists("posthocControl", .env)) {
    .control <- get("posthocControl", .env)
    if (inherits(.control, "posthocControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "posthocControl")) return(.control)
  }
  stop("cannot find posthoc related control object", call.=FALSE)
}

.posthocControlToFoceiControl <- function(env, assign=TRUE) {
  .posthocControl <- env$posthocControl
  .ui <- env$ui
  .n <- names(.posthocControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "maxOuterIterations") {
                                       return(0L)
                                     }
                                     if (n == "interaction") {
                                       return(.posthocControl$interaction)
                                     }
                                     .posthocControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.posthoc <- function(x, ...) {
  .env <- x[[1]]
  .posthocControlToFoceiControl(.env, assign=FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.posthoc <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiRandomOnIdOnly(.ui,
                                   " for the estimation routine 'posthoc'",
                                   .var.name=.ui$modelName)
  .control <- env$control
  env$posthocControl <- .control
  .foceiFamilyControl(env, ..., type="posthocControl")
  .posthocControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  rxode2::rxAssignControlValue(.ui, "maxOuterIterations", 0L)
  .ret <- .foceiFamilyReturn(env, .ui, ..., est="posthoc")
  .ret
}
attr(nlmixr2Est.posthoc, "covPresent") <- TRUE
