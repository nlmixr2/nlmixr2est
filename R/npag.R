# est="npag" -- Nonparametric Adaptive Grid (Yamada 2021).  A deterministic,
# convex nonparametric population method: the parameter distribution is a set of
# discrete support points with weights, estimated by an adaptive grid over the
# support-point locations wrapped around Burke's interior-point solve for the
# weights.  Like the imp family it reuses the FOCEI inner likelihood machinery
# (src/inner.cpp) to fill the per-subject conditional-likelihood matrix; the
# estimation loop runs in C++ (npagOuter, src/npag.cpp), dispatched from
# foceiFitCpp_ via op_focei.isNpag.
#
# M0: the control and dispatch scaffold routes est="npag" to the C++ driver,
# which currently errors "not yet implemented".  Later milestones build the
# adaptive-grid cycle.  The control layers on impmapControl so the shared FOCEI
# family plumbing (inner model, mu index maps, covariance, tables) is reused; a
# dedicated .npFamilyFit with the nonparametric knobs is introduced in a later
# milestone.

#' Control for the npag (nonparametric adaptive grid) method
#'
#' A wrapper around [impmapControl()] that reuses the shared FOCEI family
#' plumbing for the nonparametric adaptive grid engine.  The nonparametric
#' support-point knobs are added in a later milestone.
#'
#' Note: the npag objective is the nonparametric marginal log-likelihood and uses
#' a different constant convention than NONMEM/FOCEI, so its `-2LL` is NOT
#' comparable to nlmixr2's FOCEI/SAEM/FOCE `-2LL`.  Compare npag runs to each
#' other or to Pmetrics NPAG.
#'
#' @inheritParams impmapControl
#' @param points Initial Sobol grid size (support points).
#' @param cycles Maximum adaptive-grid cycles.
#' @param gammaOptimize Optimize the residual-error magnitude (gamma) each cycle
#'   (only valid for normal endpoints; censoring and transform-both-sides are
#'   supported).
#' @param ... Parameters passed to [impmapControl()].
#' @return An `impmapControl` object tagged for the npag engine.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' npagControl()
npagControl <- function(points = 2028L, cycles = 100L, gammaOptimize = TRUE, ...) {
  .ctl <- impmapControl(...)
  .ctl$est <- "npag"
  .ctl$points <- as.integer(points)
  .ctl$cycles <- as.integer(cycles)
  .ctl$gammaOptimize <- isTRUE(gammaOptimize)
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.npag <- function(control) {
  .ctl <- .npValidCtl(control, "npag")
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.npag <- function(x, ...) {
  nmObjGetControl.impmap(x, ...)
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.npag <- function(x, ...) {
  nmObjGetFoceiControl.impmap(x, ...)
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.npag <- function(env, ...) {
  .npEstCore(env, "npag", ...)
}
attr(nlmixr2Est.npag, "covPresent") <- TRUE
attr(nlmixr2Est.npag, "unbounded") <- .foUnbounded
attr(nlmixr2Est.npag, "iov") <- TRUE
attr(nlmixr2Est.npag, "mu") <- .npMuAttr

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mnpag <- function(env, ...) {
  .npEstCore(env, "mnpag", muModel="lin", ...)
}
attr(nlmixr2Est.mnpag, "covPresent") <- TRUE
attr(nlmixr2Est.mnpag, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mnpag, "iov") <- TRUE
# sugar variant is mu-referenced by definition
attr(nlmixr2Est.mnpag, "mu") <- function(control) TRUE

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.inpag <- function(env, ...) {
  .npEstCore(env, "inpag", muModel="irls", ...)
}
attr(nlmixr2Est.inpag, "covPresent") <- TRUE
attr(nlmixr2Est.inpag, "unbounded") <- .foUnbounded
attr(nlmixr2Est.inpag, "iov") <- TRUE
attr(nlmixr2Est.inpag, "mu") <- function(control) TRUE
