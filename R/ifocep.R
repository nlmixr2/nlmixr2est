#' Control options for the ifocep (IRLS FOCE+) estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' FOCE+ (no interaction, live conditional residual variance R); it is
#' `ifoce` with `foce = "foce+"` forced (see `foceiControl(foce=)`).
#'
#' @inheritSection mfoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `ifocei` instead
#' @param muModel Selects the regression variant; for `ifocepControl()`
#'   this is always `"irls"` and cannot be changed -- use `mfocepControl()`
#'   for the closed-form OLS variant.
#' @param foce FOCE residual-variance mode; for `ifocepControl()` this is
#'   always `"foce+"` and cannot be changed -- use `ifoceControl()` for
#'   `"nonmem"`
#' @return ifocepControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' ifocepControl()
ifocepControl <- function(sigdig=4,
                             ...,
                             interaction=FALSE,
                             muModel=c("irls", "lin", "none"),
                             foce="foce+") {
  .control <- foceiControl(sigdig=sigdig, ..., interaction=FALSE,
                           muModel="irls", foce="foce+")
  class(.control) <- "ifocepControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.ifocepControl <- function(control, env) {
  assign("ifocepControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifocep <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- ifocepControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("ifocepControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to ifocepControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(ifocepControl, .ctl)
  } else if (!inherits(.ctl, "ifocepControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- ifocepControl()
  } else {
    .ctl <- do.call(ifocepControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.ifocep <- function(x, ...) {
  .env <- x[[1]]
  if (exists("ifocepControl", .env)) {
    .control <- get("ifocepControl", .env)
    if (inherits(.control, "ifocepControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "ifocepControl")) return(.control)
  }
  stop("cannot find ifocep related control object", call.=FALSE)
}

.ifocepControlToFoceiControl <- function(env, assign=TRUE) {
  .ifocepControl <- env$ifocepControl
  .n <- names(.ifocepControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.ifocepControl$interaction)
                                     }
                                     .ifocepControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.ifocep <- function(x, ...) {
  .env <- x[[1]]
  .ifocepControlToFoceiControl(.env, assign=FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.ifocep <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'ifocep'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="ifocepControl")
  .ifocepControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$ifocepControl <- .control
  env$est <- "ifocep"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="ifocep")
}
attr(nlmixr2Est.ifocep, "iov") <- TRUE
attr(nlmixr2Est.ifocep, "covPresent") <- TRUE
attr(nlmixr2Est.ifocep, "unbounded") <- .foUnbounded
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.ifocep, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
