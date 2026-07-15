# est="imp" -- NONMEM METHOD=IMP: importance-sampling EM WITHOUT the per-iteration
# MAP search.  Shares the impmap kernel (src/imp.cpp impOuter); the only difference
# is the E-step proposal, which for "imp" is centered at each subject's running
# conditional mean with covariance gamma*Omega (over-dispersed vs the posterior, a
# robust importance-sampling proposal) instead of the MAP mode + inner Hessian.
# Selected by mapIter = 0 / est = "imp" (impOuter reads op_focei.isImp).

#' Control for the imp (importance-sampling EM without MAP search) method
#'
#' A convenience wrapper around [impmapControl()] with `mapIter = 0`, i.e. the
#' importance-sampling proposal is centered at the running conditional mean rather
#' than re-optimized to the MAP mode each iteration (NONMEM METHOD=IMP).  See
#' [impmapControl()] for the full parameter list.
#'
#' @inheritParams impmapControl
#' @param ... Parameters passed to [impmapControl()].
#' @return An `impmapControl` object with `mapIter = 0`.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' impControl()
impControl <- function(...) {
  .ctl <- impmapControl(...)
  .ctl$mapIter <- 0L
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.imp <- function(control) {
  .ctl <- getValidNlmixrCtl.impmap(control)
  .ctl$mapIter <- 0L
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.imp <- function(x, ...) {
  nmObjGetControl.impmap(x, ...)
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.imp <- function(x, ...) {
  nmObjGetFoceiControl.impmap(x, ...)
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.imp <- function(env, ...) {
  .ui <- env$ui
  # General (dnorm/ll) likelihoods flow through the shared FOCEI inner problem, so
  # only require transformable normality when the rxode2 build has no llik support.
  if (!rxode2hasLlik()) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'imp'", .var.name=.ui$modelName)
  }
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'imp'",
                             .var.name=.ui$modelName)
  .foceiFamilyControl(env, ..., type="impmapControl")
  # imp = the impmap kernel with no per-iteration MAP re-centering.
  env$control$mapIter <- 0L
  .control <- env$control
  on.exit({
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  env$impmapControl <- .control
  env$est <- "imp"
  .ui <- env$ui
  .impmapFamilyFit(env, .ui, ...)
}
attr(nlmixr2Est.imp, "covPresent") <- TRUE
attr(nlmixr2Est.imp, "unbounded") <- .foUnbounded
attr(nlmixr2Est.imp, "iov") <- TRUE
attr(nlmixr2Est.imp, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
