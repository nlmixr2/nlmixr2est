#' Control options for the mufocei estimation method
#'
#' This is the mu-referenced-FOCEI-family closed-form-regression (`"lin"`)
#' variant of FOCEI -- see `foceiControl(muModel=)` for what this changes.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `mufoceiControl()`
#'   this is always `"lin"` and cannot be changed -- use `irlsfoceiControl()`
#'   for the IRLS variant.
#' @return mufoceiControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mufoceiControl()
mufoceiControl <- function(sigdig=3,
                           ...,
                           muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., muModel="lin")
  class(.control) <- "mufoceiControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mufoceiControl <- function(control, env) {
  assign("mufoceiControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mufocei <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- mufoceiControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("mufoceiControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to mufoceiControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(mufoceiControl, .ctl)
  } else if (!inherits(.ctl, "mufoceiControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- mufoceiControl()
  } else {
    .ctl <- do.call(mufoceiControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mufocei <- function(x, ...) {
  .env <- x[[1]]
  if (exists("mufoceiControl", .env)) {
    .control <- get("mufoceiControl", .env)
    if (inherits(.control, "mufoceiControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "mufoceiControl")) return(.control)
  }
  stop("cannot find mufocei related control object", call.=FALSE)
}

.mufoceiControlToFoceiControl <- function(env, assign=TRUE) {
  .mufoceiControl <- env$mufoceiControl
  .n <- names(.mufoceiControl)
  .foceiControl <- setNames(lapply(.n, function(n) .mufoceiControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mufocei <- function(x, ...) {
  .env <- x[[1]]
  .mufoceiControlToFoceiControl(.env, assign=FALSE)
}
