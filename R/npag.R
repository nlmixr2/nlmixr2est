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
#' Note on residual error with a flexible support distribution: the residual
#' parameters are estimated against the nonparametric objective (see
#' \code{residOptimize}) with the support-point distribution held fixed.  Because
#' that distribution is flexible, it can absorb variability a parametric model
#' (FOCEI/SAEM) would attribute to residual error -- especially the additive term
#' of a combined additive+proportional model at low concentrations.  As a result
#' the additive coefficient of a combined error model may be estimated smaller
#' (sometimes toward zero) than the corresponding parametric fit, while the
#' proportional term and per-endpoint magnitudes are recovered well.  This is an
#' expected property of nonparametric estimation, not a convergence failure; use
#' \code{residOptimize = "none"} to hold the residual parameters at their initial
#' values if a fixed error model is desired.
#'
#' @inheritParams impmapControl
#' @param points Initial Sobol grid size (support points).
#' @param cycles Maximum adaptive-grid cycles.
#' @param gammaOptimize Use a global assay-error multiplier (gamma) as a per-cycle
#'   warm start for the overall residual magnitude, folded into the variance-scale
#'   coefficients (`add`/`prop`/`lnorm`).  The per-endpoint values, the add/prop
#'   ratio, and the transform/autocorrelation parameters come from
#'   \code{residOptimize}.  Only valid for normal endpoints; censoring and
#'   transform-both-sides are supported.
#' @param residOptimize How to estimate the residual-error thetas (every endpoint's
#'   `add`/`prop`/`lnorm`, each transform `lambda`, each `ar`) with the support
#'   points and weights held fixed, using Nelder-Mead on the nonparametric -2LL.
#'   \code{"alternate"} (default) optimizes them every cycle (block-coordinate
#'   ascent); \code{"final"} optimizes once at the converged support; \code{"none"}
#'   holds them at their initial values (only the gamma warm-start adjusts the
#'   overall magnitude).  Fixed residual parameters are always held.  The
#'   optimization uses the bounded \code{minqa::bobyqa}, honoring the ini-block
#'   lower/upper bounds of the residual parameters (e.g. keeping an additive SD
#'   non-negative); an unbounded optimizer could wander into an invalid region.
#' @param muExpand how to estimate non-mu structural fixed-effect parameters (a
#'   theta with no eta, e.g. `ke <- exp(tke)`; npag's grid otherwise covers only
#'   mu-referenced and residual/likelihood parameters).  `FALSE` (default) optimizes
#'   them directly as "regressors" in the residual step -- the bounded `bobyqa` moves
#'   them alongside the residual parameters, re-solving the ODE each candidate (they
#'   feed the states, so the ODE freeze is turned off) -- which identifies them well
#'   (e.g. recovering a clearance and a mixture proportion from poor starts).  `TRUE`
#'   instead uses the saem-style mu-expansion: inject a pseudo-eta
#'   (`ke <- exp(tke + eta.tke)`), grid-estimate, and recover it as a fixed effect at
#'   finalization (support-mean folded into the theta, injected random effect
#'   collapsed).  The regressor default usually identifies these parameters more
#'   sharply than the grid.
#' @param ... Parameters passed to [impmapControl()].
#' @return An `impmapControl` object tagged for the npag engine.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' npagControl()
npagControl <- function(points = 2028L, cycles = 100L, gammaOptimize = TRUE,
                        residOptimize = c("alternate", "final", "none"),
                        muExpand = FALSE, ...) {
  .ctl <- impmapControl(...)
  .ctl$est <- "npag"
  .ctl$points <- as.integer(points)
  .ctl$cycles <- as.integer(cycles)
  .ctl$gammaOptimize <- isTRUE(gammaOptimize)
  .ctl$residOptimize <- match.arg(residOptimize)
  .ctl$muExpand <- isTRUE(muExpand)
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
