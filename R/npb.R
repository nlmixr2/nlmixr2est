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
#' @param points Stick-breaking truncation level K (number of support points).
#' @param alpha Dirichlet-process concentration parameter.
#' @param burnin Number of burn-in Gibbs sweeps.
#' @param nsamp Number of post-burn-in Gibbs samples collected.
#' @param nchains Number of independent chains (Gelman-Rubin R-hat convergence is
#'   reported when \code{nchains > 1}).
#' @param propSd Standard deviation of the Gaussian random-walk MH proposal for
#'   the support-point locations (eta space).
#' @param seed Random seed for the sampler.
#' @param residOptimize How to estimate the residual-error thetas (every endpoint's
#'   `add`/`prop`/`lnorm`, each transform `lambda`, each `ar`) and any non-mu
#'   structural "regressor" theta, with the sampled mixing distribution held fixed,
#'   using the bounded `bobyqa` on the EXTENDED LEAST SQUARES objective (see
#'   [npagControl()]; the `log(r)` term keeps the residual from collapsing to zero and
#'   the moment warm-start gives the saem-style SD).  \code{"alternate"} (default)
#'   re-fits them during burn-in and then holds them fixed for the sampling phase (so
#'   every collected draw shares the converged residual scale); \code{"final"} holds
#'   them at their initial values through sampling and fits once at the converged draw;
#'   \code{"none"} holds them at their initial values throughout.  Fixed residual
#'   parameters are always held.  Unlike npag, npb does not optimize the assay-error
#'   multiplier (gamma); the residual thetas are fit directly.
#' @param cycles Unused for npb (kept for control compatibility).
#' @param gammaOptimize Unused for npb (kept for control compatibility).
#' @param muExpand When `TRUE`, mu-expand non-mu structural fixed-effect thetas (a
#'   theta with no eta) into grid-estimable pseudo-etas before the fit; `FALSE`
#'   (default) leaves them to the residual step.
#' @param cores Number of threads used for the parallel per-subject conditional-
#'   likelihood solves in the Gibbs sweeps.  `NULL` (default) uses the current
#'   `rxode2` thread count (`rxode2::getRxThreads()`); an integer sets the thread
#'   count for the fit (restored afterwards).  With a fixed `seed` the fit is
#'   bit-for-bit identical regardless of the thread count.
#' @param rhoend Final trust-region radius (`rhoend`) of the inner bounded
#'   `bobyqa` that fits the residual-error thetas.  A fixed default of `1e-4`,
#'   matching the optimizer convergence tolerance `10^(-sigdig)` at
#'   `sigdig = 4` (npb has no `sigdig`, so this is not derived from it).
#' @param gamma,df Declared only so they are REJECTED rather than partially
#'   matched.  `gamma` is a prefix of `gammaOptimize`, so without an explicit
#'   formal R bound `gamma = 2` to it and silently turned the assay-error
#'   optimisation off; `df` is an importance-sampling proposal control a
#'   nonparametric engine never builds.  Passing either is an error.
#' @param ... Parameters passed to [impmapControl()], for the shared FOCEI-family
#'   scaffolding only.  The importance-sampling controls are **rejected** rather
#'   than accepted -- see [npagControl()] -- as are `cycles` and `gammaOptimize`,
#'   which this signature carries but `npb` does not use.  `impSeed` is rejected
#'   pointing at `seed`.
#' @return An `impmapControl` object tagged for the npb engine.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' npbControl()
npbControl <- function(points = 50L, alpha = 1.0, burnin = 500L, nsamp = 500L,
                       nchains = 1L, propSd = 0.2, seed = 42L,
                       residOptimize = c("alternate", "final", "none"),
                       cycles = 100L,
                       gammaOptimize = FALSE, muExpand = FALSE, cores = NULL,
                       rhoend = 1e-4, gamma, df, ...) {
  # `gamma` is declared ONLY to be rejected: it is a prefix of the real formal
  # gammaOptimize, so without it R partial-matches and npbControl(gamma = 2)
  # silently sets gammaOptimize = isTRUE(2).  An exact match beats a partial one,
  # which catches the name even through a forwarding wrapper.
  # validation list is SEPARATE from what is forwarded to impmapControl()
  .npbChk <- list(...)
  .npbExp <- .npCallNames(sys.call())
  if (!missing(gamma)) .npbExp <- union(.npbExp, "gamma")
  if (!missing(df)) .npbExp <- union(.npbExp, "df")
  # cycles/gammaOptimize are FORMALS here, documented unused for npb, so they
  # never reach ... ; fold them in so the shared validator sees them
  for (.n in intersect(names(as.list(match.call())[-1L]),
                       c("cycles", "gammaOptimize"))) {
    .npbChk[[.n]] <- get(.n)
  }
  .npAssertImpCtl(.npbChk, "npb", explicit = .npbExp)
  .ctl <- impmapControl(...)
  .npAssertBuilt(.ctl, "npb")
  .ctl$est <- "npb"
  checkmate::assertNumeric(rhoend, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  .ctl$rhoend <- as.numeric(rhoend)
  .ctl$points <- as.integer(points)
  .ctl$cycles <- as.integer(cycles)
  # NA -> use the default rxode2 thread count (resolved in .npEstCore)
  .ctl$npCores <- .npAssertCores(cores)
  .ctl$gammaOptimize <- isTRUE(gammaOptimize)
  .ctl$residOptimize <- match.arg(residOptimize)
  .ctl$muExpand <- isTRUE(muExpand)
  .ctl$alpha <- as.numeric(alpha)
  .ctl$burnin <- as.integer(burnin)
  .ctl$nsamp <- as.integer(nsamp)
  .ctl$nchains <- as.integer(nchains)
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
