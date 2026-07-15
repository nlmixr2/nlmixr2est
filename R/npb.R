# est="npb" -- Nonparametric Bayes (Tatarinova 2013).  A truncated stick-breaking
# Dirichlet-process mixture sampled by a blocked Metropolis-within-Gibbs sampler;
# it yields posterior distributions / credible intervals of the mixing
# distribution that npag cannot give without bootstrap.  It shares the
# conditional-likelihood primitive and support-point representation with npag and
# reuses the FOCEI inner machinery (src/inner.cpp); the sampler runs in C++
# (npbOuter, src/npb.cpp), dispatched from foceiFitCpp_ via op_focei.isNpb.
#
# M0: control and dispatch scaffold only; the C++ driver errors "not yet
# implemented".  Later milestones build the Gibbs sampler.

#' Control for the npb (nonparametric Bayes) method
#'
#' A wrapper around [impmapControl()] that reuses the shared FOCEI family
#' plumbing for the nonparametric Bayes engine.  The stick-breaking sampler knobs
#' are added in a later milestone.
#'
#' @inheritParams impmapControl
#' @param ... Parameters passed to [impmapControl()].
#' @return An `impmapControl` object tagged for the npb engine.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' npbControl()
npbControl <- function(...) {
  .ctl <- impmapControl(...)
  .ctl$est <- "npb"
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.npb <- function(control) {
  .ctl <- getValidNlmixrCtl.impmap(control)
  .ctl$est <- "npb"
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.npb <- function(x, ...) {
  nmObjGetControl.impmap(x, ...)
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.npb <- function(x, ...) {
  nmObjGetFoceiControl.impmap(x, ...)
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.npb <- function(env, ...) {
  .npEstCore(env, "npb", ...)
}
attr(nlmixr2Est.npb, "covPresent") <- TRUE
attr(nlmixr2Est.npb, "unbounded") <- .foUnbounded
attr(nlmixr2Est.npb, "iov") <- TRUE
attr(nlmixr2Est.npb, "mu") <- .npMuAttr

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mnpb <- function(env, ...) {
  .npEstCore(env, "mnpb", muModel="lin", ...)
}
attr(nlmixr2Est.mnpb, "covPresent") <- TRUE
attr(nlmixr2Est.mnpb, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mnpb, "iov") <- TRUE
attr(nlmixr2Est.mnpb, "mu") <- function(control) TRUE

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.inpb <- function(env, ...) {
  .npEstCore(env, "inpb", muModel="irls", ...)
}
attr(nlmixr2Est.inpb, "covPresent") <- TRUE
attr(nlmixr2Est.inpb, "unbounded") <- .foUnbounded
attr(nlmixr2Est.inpb, "iov") <- TRUE
attr(nlmixr2Est.inpb, "mu") <- function(control) TRUE
