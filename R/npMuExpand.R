# saem-style mu-expansion for the nonparametric engines.
#
# npag estimates only the mu-referenced (grid) parameters and the residual/
# likelihood (err-tagged) parameters.  A non-mu-referenced STRUCTURAL fixed-effect
# theta (err==NA, no eta, not fixed, not a mixture proportion -- e.g. `ke <-
# exp(tke)`) has nowhere to be estimated: it is not a grid dimension and it feeds
# the ODE, so the frozen residual step cannot touch it.  saem handles every
# parameter uniformly by making it mu-referenced; we do the same here -- inject a
# pseudo-eta (`param = theta + eta`) so the theta becomes a grid dimension estimated
# by the support-point mean-shift.  Because the structural parameter now rides the
# grid (which solves the ODE anyway) and not the residual step, the residual step
# keeps optimizing only err-tagged parameters and its ODE freeze stays valid.
#
# The injected eta starts with a modest variance so the initial Sobol box is wide
# enough to locate the theta; for a true fixed effect the nonparametric grid
# collapses that dimension toward a point mass (the "fix eta to 0" limit).

#' Substitute `theta` with `(theta + eta)` inside a single model expression.
#' @noRd
.npSubstThetaAddEta <- function(expr, theta, eta) {
  do.call("substitute",
          list(expr, setNames(list(bquote(.(as.name(theta)) + .(as.name(eta)))), theta)))
}

#' Mu-expand every non-mu structural fixed-effect theta of a nonparametric-engine
#' model so it becomes grid-estimable.  Returns the (possibly rewritten) ui.
#' @param ui rxode2 ui
#' @param initVar initial variance for each injected pseudo-eta (box half-width is
#'   gridWidth * sqrt(initVar) on the transformed scale)
#' @noRd
.npMuExpand <- function(ui, initVar = 0.25) {
  nlmixr2global$npMuExpandPairs <- NULL   # clear any stale pairs from a prior fit
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta), , drop = FALSE]
  .mr <- ui$muRefDataFrame
  .isFix <- !is.na(.th$fix) & .th$fix
  .isMu <- .th$name %in% .mr$theta
  .isErr <- !is.na(.th$err)
  # a mixture proportion (ui$mixProbs) is not a candidate: it is estimated by the
  # proportion EM, not the grid, and must not be injected/collapsed.
  .isMix <- .th$name %in% ui$mixProbs
  .cand <- .th$name[!.isFix & !.isMu & !.isMix & !.isErr]
  if (length(.cand) == 0L) return(ui)
  .allNames <- .iniDf$name
  .injEta <- character(0)   # injected pseudo-eta names
  .injTh <- character(0)    # the theta each pairs to
  for (.tn in .cand) {
    # a unique pseudo-eta name for this theta
    .etaN <- paste0("eta.", .tn)
    .k <- 1L
    while (.etaN %in% .allNames) { .etaN <- paste0("eta.", .tn, ".", .k); .k <- .k + 1L }
    .allNames <- c(.allNames, .etaN)
    # rewrite every model expression that uses the theta: theta -> (theta + eta)
    .le <- ui$lstExpr
    .used <- FALSE
    for (i in seq_along(.le)) {
      if (.tn %in% all.vars(.le[[i]])) {
        .new <- .npSubstThetaAddEta(.le[[i]], .tn, .etaN)
        ui <- eval(bquote(rxode2::model(ui, .(.new))))
        .used <- TRUE
      }
    }
    if (.used) {
      # inject the pseudo-eta with a FIXED omega: it represents a fixed effect
      # (eta -> 0), so it must NOT be a free parameter in the omega objective -- it is
      # only a support (grid) dimension for locating the theta.  Excluding it from the
      # free omega parameterization also sidesteps the npag+mixture multi-free-eta
      # omega setup limitation (same pattern IOV uses).  initVar sets the grid box
      # half-width (gridWidth * sqrt(initVar)); finalization folds its support-mean
      # into the theta and collapses its reported BSV to ~0.
      ui <- eval(bquote(rxode2::ini(ui, .(as.name(.etaN)) ~ fix(.(initVar)))))
      .injEta <- c(.injEta, .etaN); .injTh <- c(.injTh, .tn)
    }
  }
  # record the injected (eta, theta) pairs so finalization can recover each theta as
  # a fixed effect: fold the injected eta's support-mean into its theta and collapse
  # that eta's random effect (the "fix eta to 0" step).  A per-fit package-global
  # survives the control rebuild between here (a preprocessing hook) and .npFamilyFit
  # (control-stored values do not).
  nlmixr2global$npMuExpandPairs <- list(eta = .injEta, theta = .injTh)
  ui
}

# Which est strings are nonparametric-engine families (and their mu/irls sugar).
.npEstFamily <- c("npag", "mnpag", "inpag", "npb", "mnpb", "inpb")

#' Preprocessing hook: mu-expand non-mu structural fixed-effect thetas for the
#' nonparametric engines, before the rest of the pipeline builds on the ui.  Doing
#' it here (rather than mutating the ui mid-setup) keeps the injected pseudo-etas
#' consistent through covariate/mu processing and model compilation.  Off when
#' control$muExpand is FALSE.
#' @inheritParams nlmixr2
#' @return list(ui=) when the model was expanded, else NULL
#' @export
#' @author Matthew L. Fidler
.nlmixr0preProcessNpMuExpand <- function(ui, est, data, control) {
  if (is.null(est) || !(est %in% .npEstFamily)) return(NULL)
  # always clear the injected-pair record so a prior muExpand fit cannot leak into a
  # muExpand=FALSE fit (which returns before .npMuExpand clears it itself).
  nlmixr2global$npMuExpandPairs <- NULL
  if (!is.null(control) && isFALSE(control$muExpand)) return(NULL)
  .ui2 <- .npMuExpand(ui)
  if (identical(.ui2, ui)) return(NULL)
  list(ui = .ui2)
}

preProcessHooksAdd(".nlmixr0preProcessNpMuExpand", .nlmixr0preProcessNpMuExpand)
