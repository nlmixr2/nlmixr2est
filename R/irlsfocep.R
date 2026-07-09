#' Control options for the irlsfocep (IRLS FOCE+) estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' FOCE+ (no interaction, live conditional residual variance R); it is
#' `irlsfoce` with `foce = "foce+"` forced (see `foceiControl(foce=)`).
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `irlsfocei` instead
#' @param muModel Selects the regression variant; for `irlsfocepControl()`
#'   this is always `"irls"` and cannot be changed -- use `mufocepControl()`
#'   for the closed-form OLS variant.
#' @param foce FOCE residual-variance mode; for `irlsfocepControl()` this is
#'   always `"foce+"` and cannot be changed -- use `irlsfoceControl()` for
#'   `"nonmem"`
#' @return irlsfocepControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' irlsfocepControl()
irlsfocepControl <- function(sigdig=3,
                             ...,
                             interaction=FALSE,
                             muModel=c("irls", "lin", "none"),
                             foce="foce+") {
  .control <- foceiControl(sigdig=sigdig, ..., interaction=FALSE,
                           muModel="irls", foce="foce+")
  class(.control) <- "irlsfocepControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.irlsfocepControl <- function(control, env) {
  assign("irlsfocepControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.irlsfocep <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- irlsfocepControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("irlsfocepControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to irlsfocepControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(irlsfocepControl, .ctl)
  } else if (!inherits(.ctl, "irlsfocepControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- irlsfocepControl()
  } else {
    .ctl <- do.call(irlsfocepControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.irlsfocep <- function(x, ...) {
  .env <- x[[1]]
  if (exists("irlsfocepControl", .env)) {
    .control <- get("irlsfocepControl", .env)
    if (inherits(.control, "irlsfocepControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "irlsfocepControl")) return(.control)
  }
  stop("cannot find irlsfocep related control object", call.=FALSE)
}

.irlsfocepControlToFoceiControl <- function(env, assign=TRUE) {
  .irlsfocepControl <- env$irlsfocepControl
  .n <- names(.irlsfocepControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.irlsfocepControl$interaction)
                                     }
                                     .irlsfocepControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.irlsfocep <- function(x, ...) {
  .env <- x[[1]]
  .irlsfocepControlToFoceiControl(.env, assign=FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.irlsfocep <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'irlsfocep'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="irlsfocepControl")
  .irlsfocepControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$irlsfocepControl <- .control
  env$est <- "irlsfocep"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="irlsfocep")
}
attr(nlmixr2Est.irlsfocep, "iov") <- TRUE
attr(nlmixr2Est.irlsfocep, "covPresent") <- TRUE
attr(nlmixr2Est.irlsfocep, "unbounded") <- .foUnbounded
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.irlsfocep, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
