#' Control options for the mufocep (mu-referenced FOCE+) estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' FOCE+ (no interaction, live conditional residual variance R); it is
#' `mufoce` with `foce = "foce+"` forced (see `foceiControl(foce=)`).
#'
#' @inheritSection mufoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `mufocei` instead
#' @param muModel Selects the regression variant; for `mufocepControl()`
#'   this is always `"lin"` and cannot be changed -- use `irlsfocepControl()`
#'   for the IRLS variant.
#' @param foce FOCE residual-variance mode; for `mufocepControl()` this is
#'   always `"foce+"` and cannot be changed -- use `mufoceControl()` for
#'   `"nonmem"`
#' @return mufocepControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mufocepControl()
mufocepControl <- function(sigdig=3,
                           ...,
                           interaction=FALSE,
                           muModel=c("lin", "irls", "none"),
                           foce="foce+") {
  .control <- foceiControl(sigdig=sigdig, ..., interaction=FALSE,
                           muModel="lin", foce="foce+")
  class(.control) <- "mufocepControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mufocepControl <- function(control, env) {
  assign("mufocepControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mufocep <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- mufocepControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("mufocepControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to mufocepControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(mufocepControl, .ctl)
  } else if (!inherits(.ctl, "mufocepControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- mufocepControl()
  } else {
    .ctl <- do.call(mufocepControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mufocep <- function(x, ...) {
  .env <- x[[1]]
  if (exists("mufocepControl", .env)) {
    .control <- get("mufocepControl", .env)
    if (inherits(.control, "mufocepControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "mufocepControl")) return(.control)
  }
  stop("cannot find mufocep related control object", call.=FALSE)
}

.mufocepControlToFoceiControl <- function(env, assign=TRUE) {
  .mufocepControl <- env$mufocepControl
  .n <- names(.mufocepControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.mufocepControl$interaction)
                                     }
                                     .mufocepControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mufocep <- function(x, ...) {
  .env <- x[[1]]
  .mufocepControlToFoceiControl(.env, assign=FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.mufocep <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'mufocep'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="mufocepControl")
  .mufocepControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$mufocepControl <- .control
  env$est <- "mufocep"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="mufocep")
}
attr(nlmixr2Est.mufocep, "iov") <- TRUE
attr(nlmixr2Est.mufocep, "covPresent") <- TRUE
attr(nlmixr2Est.mufocep, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.mufocep, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
