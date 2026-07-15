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
#' Note: the npb objective is the nonparametric marginal log-likelihood and uses
#' a different constant convention than NONMEM/FOCEI, so its `-2LL` is NOT
#' comparable to nlmixr2's FOCEI/SAEM/FOCE `-2LL`.  Compare npb runs to each other
#' or to Pmetrics NPAG.
#'
#' @inheritParams impmapControl
#' @param points Stick-breaking truncation level K (number of support points).
#' @param alpha Dirichlet-process concentration parameter.
#' @param burnin Number of burn-in Gibbs sweeps.
#' @param nsamp Number of post-burn-in Gibbs samples collected.
#' @param propSd Standard deviation of the Gaussian random-walk MH proposal for
#'   the support-point locations (eta space).
#' @param seed Random seed for the sampler.
#' @param cycles Unused for npb (kept for control compatibility).
#' @param gammaOptimize Unused for npb (kept for control compatibility).
#' @param ... Parameters passed to [impmapControl()].
#' @return An `impmapControl` object tagged for the npb engine.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' npbControl()
npbControl <- function(points = 50L, alpha = 1.0, burnin = 500L, nsamp = 500L,
                       propSd = 0.2, seed = 42L, cycles = 100L,
                       gammaOptimize = FALSE, ...) {
  .ctl <- impmapControl(...)
  .ctl$est <- "npb"
  .ctl$points <- as.integer(points)
  .ctl$cycles <- as.integer(cycles)
  .ctl$gammaOptimize <- isTRUE(gammaOptimize)
  .ctl$alpha <- as.numeric(alpha)
  .ctl$burnin <- as.integer(burnin)
  .ctl$nsamp <- as.integer(nsamp)
  .ctl$propSd <- as.numeric(propSd)
  .ctl$seed <- as.integer(seed)
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.npb <- function(control) {
  .ctl <- .npValidCtl(control, "npb")
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
