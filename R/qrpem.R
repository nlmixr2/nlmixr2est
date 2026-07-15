# est="qrpem" -- QRPEM (quasi-random parametric EM, Leary & Dunlavey PAGE
# 2012): the impmap importance-sampling EM with Sobol quasi-random importance
# samples and the SIR-accelerated non-mu/residual M-step.  Sugar over the
# impmap kernel: qrpemControl() = impmapControl(qr=TRUE, sir=TRUE, ...).

#' Control for the qrpem (quasi-random parametric EM) estimation method
#'
#' A convenience wrapper around [impmapControl()] defaulting `qr=TRUE` (Sobol
#' quasi-random importance samples) and `sir=TRUE` (SIR-accelerated non-mu /
#' residual-error M-step); explicitly supplied arguments win.  See
#' [impmapControl()] for the full parameter list.
#'
#' Note this is not know to be the same as the QRPEM implementation in
#' Phoenix NLME since the details of their method are not public.
#' However, this matches the QRPEM method of using quasi-random
#' parametric EM and SIR accelerated parameter convergence
#' described in Leary & Dunlavey (2012) PAGE 2012, 19(1): 1-6.
#'
#' @inheritParams impmapControl
#' @param ... Parameters passed to [impmapControl()].
#' @return An `impmapControl` object with the QRPEM defaults.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' qrpemControl()
qrpemControl <- function(..., qr=TRUE, sir=TRUE) {
  impmapControl(..., qr=qr, sir=sir)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.qrpem <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) return(qrpemControl())
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) {
    return(do.call("qrpemControl", .ctl))
  }
  # an explicit impmapControl()/qrpemControl() object passes through the
  # impmap validation unchanged (its qr/sir choices are respected)
  getValidNlmixrCtl.impmap(control)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.qrpem <- function(x, ...) {
  nmObjGetControl.impmap(x, ...)
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.qrpem <- function(x, ...) {
  nmObjGetFoceiControl.impmap(x, ...)
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.qrpem <- function(env, ...) {
  .ui <- env$ui
  # General (dnorm/ll) likelihoods flow through the shared FOCEI inner problem,
  # so only require transformable normality when rxode2 has no llik support.
  if (!rxode2hasLlik()) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'qrpem'", .var.name=.ui$modelName)
  }
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'qrpem'",
                             .var.name=.ui$modelName)
  .foceiFamilyControl(env, ..., type="impmapControl")
  .control <- env$control
  on.exit({
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  env$impmapControl <- .control
  env$est <- "qrpem"
  .ui <- env$ui
  .impmapFamilyFit(env, .ui, ...)
}
attr(nlmixr2Est.qrpem, "covPresent") <- TRUE
attr(nlmixr2Est.qrpem, "unbounded") <- .foUnbounded
attr(nlmixr2Est.qrpem, "iov") <- TRUE
attr(nlmixr2Est.qrpem, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
