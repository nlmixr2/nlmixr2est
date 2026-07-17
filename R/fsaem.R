#' Control options for the fsaem (fast SAEM) estimation method
#'
#' `fsaem` is the fast-SAEM (f-SAEM) of Karimi, Lavielle and Moulines
#' (2020).  It is sugar over `saem` with `saemControl(fast=TRUE)` forced:
#' the MCMC simulation step samples the individual random effects from an
#' independent Metropolis-Hastings proposal centered at each subject's
#' conditional MAP estimate, accelerating the early SAEM iterations.  All
#' other options are the `saemControl()` options; see there (in particular
#' `fastKernel`, `fastCov`, `fastIter` and `fastLik`) for the fast-specific
#' tuning knobs.
#'
#' @inheritParams saemControl
#' @param ... Parameters used in the default `saemControl()`
#' @param fast Always `TRUE` for `fsaemControl()` and cannot be changed --
#'   use `saemControl(fast=FALSE)` (or `est="saem"`) for standard SAEM.
#' @return fsaemControl object (a `saemControl` with `fast=TRUE`)
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' fsaemControl()
fsaemControl <- function(..., fast=TRUE) {
  .control <- saemControl(..., fast=TRUE)
  class(.control) <- "fsaemControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.fsaemControl <- function(control, env) {
  assign("fsaemControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.fsaem <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- fsaemControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) {
    .ctl <- do.call("fsaemControl", .ctl)
  }
  if (inherits(.ctl, "saemControl") && !inherits(.ctl, "fsaemControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to fsaemControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(fsaemControl, .ctl)
  } else if (!inherits(.ctl, "fsaemControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- fsaemControl()
  } else {
    class(.ctl) <- NULL
    .ctl <- do.call(fsaemControl, .ctl)
  }
  .ctl
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.fsaem <- function(env, ...) {
  .ui <- env$ui
  # fsaem can fit a general log-likelihood endpoint (ll() ~ expr) via the FOCEi
  # inner (the fast kernel supplies the observation likelihood); only require a
  # (transformably) normal model for the ordinary continuous case.
  if (!.fsaemGeneralLik(.ui)) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'fsaem'", .var.name=.ui$modelName)
  }
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'fsaem'",
                             .var.name=.ui$modelName)
  rxode2::assertRxUiMixedOnly(.ui, .noRandomEffectMsg("fsaem"), .var.name=.ui$modelName)
  rxode2::warnRxBounded(.ui, " which are ignored in 'fsaem'", .var.name=.ui$modelName)
  if (length(.ui$mixProbs) > 0) {
    message("mixture SAEM computation scales with the number of sub-populations")
  }
  .control <- env$control
  if (is.null(.control)) {
    .control <- fsaemControl()
  }
  if (!inherits(.control, "fsaemControl")) {
    .control <- do.call(fsaemControl, .control)
  }
  env$fsaemControl <- .control
  # fsaem is sugar for saemControl(fast=TRUE): run through the standard SAEM
  # machinery with a plain saemControl carrying fast=TRUE.
  .saem <- .control
  class(.saem) <- "saemControl"
  env$control <- .saem
  .saemFamilyControl(env, ...)
  on.exit({
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  .saemFamilyFit(env, ...)
}
attr(nlmixr2Est.fsaem, "covPresent") <- TRUE
attr(nlmixr2Est.fsaem, "unbounded") <- TRUE
attr(nlmixr2Est.fsaem, "mu") <- TRUE
attr(nlmixr2Est.fsaem, "iov") <- TRUE
