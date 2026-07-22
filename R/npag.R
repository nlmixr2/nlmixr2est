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
#' @param points Initial Sobol grid size (support points).  `NULL` (default) picks
#'   it automatically from the number of support-point dimensions (etas):
#'   `max(2028, 512 * n_eta)` -- a fixed grid (Pmetrics uses 2028) covers a
#'   low-dimensional model but grows sparse and can collapse in high dimensions, so
#'   the auto size floors at 2028 and scales up per added eta.  Supply an integer to
#'   override.
#' @param cycles Maximum adaptive-grid cycles.
#' @param gammaOptimize Use a global assay-error multiplier (gamma) as a per-cycle
#'   warm start for the overall residual magnitude, folded into the variance-scale
#'   coefficients (`add`/`prop`/`lnorm`).  The per-endpoint values, the add/prop
#'   ratio, and the transform/autocorrelation parameters come from
#'   \code{residOptimize}.  Only valid for normal endpoints; censoring and
#'   transform-both-sides are supported.
#' @param residOptimize How to estimate the residual-error thetas (every endpoint's
#'   `add`/`prop`/`lnorm`, each transform `lambda`, each `ar`) with the support points
#'   and weights held fixed, using bounded \code{minqa::bobyqa} on the EXTENDED LEAST
#'   SQUARES objective `sum_obs((f-dv)^2/r + log(r))` at the posterior-mean etas.  The
#'   `log(r)` term keeps the residual from drifting to zero on a flexible support (which
#'   the marginal likelihood would reward), giving the saem/focei residual; each
#'   variance scale is warm-started from the per-endpoint moment (additive SD from
#'   `sqrt(mean(err^2))`, proportional from `sqrt(mean((err/f)^2))`, on the transform-
#'   both-sides scale).  \code{"alternate"} (default) optimizes every cycle (block-
#'   coordinate ascent); \code{"final"} optimizes once at the converged support;
#'   \code{"none"} holds them at their initial values.  Fixed residual parameters are
#'   always held.  After the residual thetas converge, a final adaptive-grid pass
#'   re-optimizes the support with them held constant so the support remains the
#'   nonparametric MLE (D(F) ~ 0) for the fitted residual.
#' @param muExpand how to estimate non-mu structural fixed-effect parameters (a
#'   theta with no eta, e.g. `ke <- exp(tke)`; npag's grid otherwise covers only
#'   mu-referenced and residual/likelihood parameters).  `FALSE` (default) optimizes
#'   them directly as "regressors" in the residual step -- `bobyqa` moves them alongside
#'   the residual parameters, re-deriving the posterior-mean etas each candidate (so the
#'   eta grid cannot stale-absorb the structural shift) -- which identifies them well
#'   (e.g. recovering a clearance from a poor start).  Not available for mix() models
#'   (the ELS step is not mixture-aware; component parameters are held).  `TRUE`
#'   instead uses the saem-style mu-expansion: inject a pseudo-eta
#'   (`ke <- exp(tke + eta.tke)`), grid-estimate, and recover it as a fixed effect at
#'   finalization (support-mean folded into the theta, injected random effect
#'   collapsed).  The regressor default usually identifies these parameters more
#'   sharply than the grid.
#' @param gridWidth support-point box half-width, in initial-eta SDs, for the
#'   `gridBounds="auto"` grid (default 4).  A narrower box focuses the initial Sobol
#'   grid on the plausible region -- useful for high-dimensional models where a wide
#'   box wastes points on near-zero-density support (which can collapse the fit).
#' @param gridBounds how to set the initial support-point box: `"auto"` (default)
#'   uses `+/- gridWidth * initial eta SD`; `"ini"` uses each mu-referenced
#'   parameter's ini-block lower/upper bounds where finite (else auto); `"both"`
#'   uses the ini bounds when present and auto otherwise.  For a high-dimensional
#'   model, bounded ini estimates + `"ini"` keep the grid in range.
#' @param dfScan Size of the Sobol scan used for the D(F) global-optimality
#'   certificate: `-1` (default) auto-sizes it to `max(2048, 2 * points)`, `0`
#'   skips the certificate (`npagDF` is `NA`), and a positive value sets an
#'   explicit scan size.  The scan does not affect the fit, only the reported
#'   certificate; a smaller scan is faster.
#' @param cores Number of threads used for the parallel per-subject conditional-
#'   likelihood solves.  `NULL` (default) uses the current `rxode2` thread count
#'   (`rxode2::getRxThreads()`); an integer sets the thread count for the fit
#'   (restored afterwards).  Results are independent of the thread count.
#' @param rhoend Final trust-region radius (`rhoend`) of the inner bounded
#'   `bobyqa` that fits the residual-error thetas each cycle.  A fixed default of
#'   `1e-4`, matching the optimizer convergence tolerance `10^(-sigdig)` at the
#'   default `sigdig = 4` (npag has no `sigdig`, so this is not derived from it).
#' @param ... Parameters passed to [impmapControl()].
#' @return An `impmapControl` object tagged for the npag engine.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' npagControl()
npagControl <- function(points = NULL, cycles = 100L, gammaOptimize = TRUE,
                        residOptimize = c("alternate", "final", "none"),
                        muExpand = FALSE, gridWidth = 4,
                        gridBounds = c("auto", "ini", "both"), dfScan = -1L,
                        cores = NULL, rhoend = 1e-4, ...) {
  .ctl <- impmapControl(...)
  .ctl$est <- "npag"
  checkmate::assertNumeric(rhoend, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  .ctl$rhoend <- as.numeric(rhoend)
  # NULL -> auto (scaled with the number of dimensions in .npEstCore); NA is the
  # sentinel that survives the control round-trip.
  .ctl$points <- if (is.null(points)) NA_integer_ else as.integer(points)
  .ctl$cycles <- as.integer(cycles)
  .ctl$dfScan <- .npAssertDfScan(dfScan)
  # NA -> use the default rxode2 thread count (resolved in .npEstCore)
  .ctl$npCores <- .npAssertCores(cores)
  .ctl$gammaOptimize <- isTRUE(gammaOptimize)
  .ctl$residOptimize <- match.arg(residOptimize)
  .ctl$muExpand <- isTRUE(muExpand)
  .ctl$gridWidth <- as.numeric(gridWidth)
  .ctl$gridBounds <- match.arg(gridBounds)
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
