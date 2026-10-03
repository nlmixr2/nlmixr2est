#' Control options for the mfocep (mu-referenced FOCE+) estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' FOCE+ (no interaction, live conditional residual variance R); it is
#' `mfoce` with `foce = "foce+"` forced (see `foceiControl(foce=)`).
#'
#' @inheritSection mfoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `mfocei` instead
#' @param muModel Selects the regression variant; for `mfocepControl()`
#'   this is always `"lin"` and cannot be changed -- use `ifocepControl()`
#'   for the IRLS variant.
#' @param foce FOCE residual-variance mode; for `mfocepControl()` this is
#'   always `"foce+"` and cannot be changed -- use `mfoceControl()` for
#'   `"nonmem"`
#' @return mfocepControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mfocepControl()
mfocepControl <- function(sigdig = 3, ..., interaction = FALSE, muModel = c("lin", "irls", "none"), foce = "foce+") {
  .control <- foceiControl(sigdig = sigdig, ..., interaction = FALSE, muModel = "lin", foce = "foce+")
  class(.control) <- "mfocepControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mfocepControl <- function(control, env) assign("mfocepControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfocep <- function(control) .foceiFamilyValidCtl(control, "mfocepControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mfocep <- function(x, ...) .nmObjGetControlByClass(x, "mfocepControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mfocep <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "mfocepControl", assign = FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.mfocep <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'mfocep'", .var.name = .ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type = "mfocepControl")
  .foceiFamilyControlToFoceiControl(env, "mfocepControl")
  on.exit({
    if (exists("control", envir = .ui)) {
      rm("control", envir = .ui)
    }
  })
  env$mfocepControl <- .control
  env$est <- "mfocep"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est = "mfocep")
}
attr(nlmixr2Est.mfocep, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.mfocep, "iov") <- TRUE
attr(nlmixr2Est.mfocep, "covPresent") <- TRUE
attr(nlmixr2Est.mfocep, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.mfocep, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
