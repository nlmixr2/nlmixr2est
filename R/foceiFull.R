# Convenience "f*" estimation methods: each is exactly its base Laplace/AGQ
# method with foceiControl(fast = TRUE, innerHessian = "conditional") as the
# default -- the full analytic conditional inner Hessian rather than the
# Gauss-Newton FOCEI one.  They are thin delegates, like the "*f" fast methods
# in R/foceiFast.R: the control validator forces both options through the base
# control constructor and the estimator dispatches to the base method.  The fit
# reports "Full Laplace"/"Full AGQ" because the label is driven by the
# conditional-curvature flag (src/inner.cpp), not by the est= name.
#
# The conditional curvature is only defined for Gaussian endpoints, so these
# methods refuse a generalized-likelihood model -- use the base
# laplace/agq method for those.

#' Force the full conditional inner Hessian through a base control constructor
#'
#' @param control the `getValidNlmixrControl` wrapper list
#' @param ctlFun the base control constructor (e.g. `laplaceControl`)
#' @return a base control object with `fast = TRUE` and
#'   `innerHessian = "conditional"`
#' @noRd
.foceiFullCtl <- function(control, ctlFun) {
  .ctl <- control[[1]]
  .l <- if (is.null(.ctl)) list() else unclass(.ctl)
  # a defaulted (not user-chosen) outer optimizer re-defaults under fast=TRUE
  # (nlminb -> lbfgsb3c); an explicit outerOpt is kept as given
  if (isTRUE(.l$outerOptDefault)) {
    .l$outerOpt <- NULL
    .l$outerOptTxt <- NULL
    .l$outerOptFun <- NULL
    .l$outerOptDefault <- NULL
  }
  .l$fast <- TRUE
  .l$innerHessian <- "conditional"
  do.call(ctlFun, .l)
}

## ---- flaplace / mflaplace / iflaplace -------------------------------------

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.flaplace <- function(control) .foceiFullCtl(control, laplaceControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.flaplace <- function(env, ...) nlmixr2Est.laplace(env, ...)
attr(nlmixr2Est.flaplace, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.flaplace, "covPresent") <- TRUE
attr(nlmixr2Est.flaplace, "unbounded") <- .foUnbounded
attr(nlmixr2Est.flaplace, "iov") <- TRUE

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mflaplace <- function(control) .foceiFullCtl(control, mlaplaceControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mflaplace <- function(env, ...) nlmixr2Est.mlaplace(env, ...)
attr(nlmixr2Est.mflaplace, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.mflaplace, "covPresent") <- TRUE
attr(nlmixr2Est.mflaplace, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mflaplace, "iov") <- TRUE
attr(nlmixr2Est.mflaplace, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.iflaplace <- function(control) .foceiFullCtl(control, ilaplaceControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.iflaplace <- function(env, ...) nlmixr2Est.ilaplace(env, ...)
attr(nlmixr2Est.iflaplace, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.iflaplace, "covPresent") <- TRUE
attr(nlmixr2Est.iflaplace, "unbounded") <- .foUnbounded
attr(nlmixr2Est.iflaplace, "iov") <- TRUE
attr(nlmixr2Est.iflaplace, "mu") <- .foceiFastMuAttr

## ---- fagq / mfagq / ifagq -------------------------------------------------

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.fagq <- function(control) .foceiFullCtl(control, agqControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.fagq <- function(env, ...) nlmixr2Est.agq(env, ...)
attr(nlmixr2Est.fagq, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.fagq, "covPresent") <- TRUE
attr(nlmixr2Est.fagq, "unbounded") <- .foUnbounded
attr(nlmixr2Est.fagq, "iov") <- TRUE

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfagq <- function(control) .foceiFullCtl(control, magqControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mfagq <- function(env, ...) nlmixr2Est.magq(env, ...)
attr(nlmixr2Est.mfagq, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.mfagq, "covPresent") <- TRUE
attr(nlmixr2Est.mfagq, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mfagq, "iov") <- TRUE
attr(nlmixr2Est.mfagq, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifagq <- function(control) .foceiFullCtl(control, iagqControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.ifagq <- function(env, ...) nlmixr2Est.iagq(env, ...)
attr(nlmixr2Est.ifagq, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.ifagq, "covPresent") <- TRUE
attr(nlmixr2Est.ifagq, "unbounded") <- .foUnbounded
attr(nlmixr2Est.ifagq, "iov") <- TRUE
attr(nlmixr2Est.ifagq, "mu") <- .foceiFastMuAttr
