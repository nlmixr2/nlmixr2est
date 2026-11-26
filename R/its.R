# est="its" -- iterative two-stage EM.  The impmap kernel (src/imp.cpp) with a
# deterministic E-step: each subject's MAP mode and inverse inner Hessian are
# the conditional mean and covariance, so no importance samples are drawn.

#' Control for the iterative two-stage (ITS) estimation method
#'
#' Each iteration finds every subject's MAP etas and conditional covariance
#' (the inverse Hessian of the individual objective), then updates the
#' mu-referenced thetas and Omega in closed form (Omega is the empirical
#' covariance of the etas plus the mean conditional covariance); the remaining
#' thetas and the residual error take a Newton step on the conditional
#' likelihood at the MAP etas.  The reported objective is a FOCE(i) evaluation
#' at the final estimates.
#'
#' ITS is fast and is often used as a warm start for another method, e.g.
#' `nlmixr2(fitIts, est = "focei")`.
#'
#' @param maxIter Maximum number of ITS iterations.
#' @param tol Relative tolerance on the change of the parameters (floored at
#'   1) and of the objective between iterations.
#' @param interaction Use the FOCEi (`TRUE`) or FOCE (`FALSE`) inner Hessian
#'   and final objective.
#' @param covMethod Covariance method computed at the final estimates, as in
#'   [foceiControl()].
#' @inheritParams impmapControl
#' @param ... Other parameters passed to [impmapControl()] and from there to
#'   [foceiControl()] (inner-problem and ODE controls).
#' @return An `impmapControl` object marked for ITS.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' itsControl()
itsControl <- function(
  maxIter = 50L,
  tol = 1e-4,
  interaction = TRUE,
  covMethod = c("r,s", "analytic", "r", "s", "sa", ""),
  sigdig = 3,
  ...
) {
  .dots <- list(...)
  # a round-tripped control carries the impmapControl names; they stand in
  # for the arguments that were not given
  if (missing(maxIter) && !is.null(.dots[["nIter"]])) {
    maxIter <- .dots[["nIter"]]
  }
  if (missing(tol) && !is.null(.dots[["ctol"]])) {
    tol <- .dots[["ctol"]]
  }
  if (missing(interaction) && !is.null(.dots[["interaction"]])) {
    interaction <- .dots[["interaction"]]
  }
  if (missing(covMethod) && !is.null(.dots[["covMethod"]])) {
    covMethod <- .dots[["covMethod"]]
  }
  checkmate::assertIntegerish(maxIter, lower = 0, len = 1, any.missing = FALSE, .var.name = "maxIter")
  checkmate::assertNumber(tol, lower = 0, finite = TRUE, .var.name = "tol")
  if (is.character(covMethod)) {
    covMethod <- .covMethodArg(covMethod, match.arg(covMethod))
  }
  .its <- list(
    nIter = as.integer(maxIter),
    ctol = as.double(tol),
    interaction = interaction,
    covMethod = covMethod,
    isample = 1L,
    mapIter = 1L,
    nBurn = 0L,
    gamma = 1.0,
    gammaMethod = "global",
    gammaRule = "floor",
    proposal = "normal",
    df = 0,
    auto = FALSE,
    iaccept = 0,
    nConvWindow = 1L,
    qr = FALSE,
    sir = FALSE,
    sirSample = 1L,
    its = TRUE
  )
  .dots[c(names(.its), "sigdig", "gammaMethodUser")] <- NULL
  do.call(impmapControl, c(list(sigdig = sigdig), .dots, .its))
}

#' Restore impmap's sampling defaults on an ITS control used by another method
#'
#' @param ctl validated impmapControl
#' @param est target estimation method
#' @return `ctl`, unchanged unless it is an ITS control and `est` is not "its"
#' @noRd
.itsControlRestoreSampling <- function(ctl, est) {
  if (!isTRUE(ctl$its) || identical(est, "its")) {
    return(ctl)
  }
  .def <- if (identical(est, "qrpem")) qrpemControl() else impmapControl()
  for (.n in c(
    "qr",
    "sir",
    "nIter",
    "isample",
    "gamma",
    "gammaMethod",
    "gammaRule",
    "proposal",
    "df",
    "auto",
    "iaccept",
    "nConvWindow",
    "ctol",
    "sirSample"
  )) {
    ctl[[.n]] <- .def[[.n]]
  }
  ctl$gammaMethodUser <- NULL
  ctl$its <- NULL
  .minfo(paste0("`est=\"", est, "\"` restores the importance-sampling defaults of an ITS control"))
  ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.its <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) {
    return(itsControl())
  }
  if (
    inherits(.ctl, "foceiControl") ||
      inherits(.ctl, "impmapControl") ||
      (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
  ) {
    # every ITS-defining field is forced by itsControl(), so a control from
    # another method (or a completed fit) is adopted for its shared options
    .lst <- unclass(.ctl)
    if (!is.null(attr(.ctl, "class")) && !isTRUE(.lst[["its"]])) {
      # only the foceiControl options are shared with another method
      .lst <- .lst[intersect(names(.lst), setdiff(names(formals(foceiControl)), "covMethod"))]
    }
    .lst$est <- NULL
    return(do.call(itsControl, .lst))
  }
  .minfo("invalid control for `est=\"its\"`, using default")
  itsControl()
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.its <- function(x, ...) {
  nmObjGetControl.impmap(x, ...)
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.its <- function(x, ...) {
  nmObjGetFoceiControl.impmap(x, ...)
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.its <- function(env, ...) {
  .ui <- env$ui
  if (!rxode2hasLlik()) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'its'", .var.name = .ui$modelName)
  }
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'its'", .var.name = .ui$modelName)
  rxode2::assertRxUiMixedOnly(.ui, .noRandomEffectMsg("its"), .var.name = .ui$modelName)
  .foceiFamilyControl(env, ..., type = "impmapControl")
  .control <- env$control
  on.exit(
    {
      if (is.environment(.ui) && exists("control", envir = .ui, inherits = FALSE)) {
        rm("control", envir = .ui)
      }
    },
    add = TRUE
  )
  env$impmapControl <- .control
  env$est <- "its"
  .ui <- env$ui
  .impmapFamilyFit(env, .ui, ...)
}
attr(nlmixr2Est.its, "nlmixr2LikPenalty") <- "full"
attr(nlmixr2Est.its, "covPresent") <- TRUE
attr(nlmixr2Est.its, "unbounded") <- .foUnbounded
attr(nlmixr2Est.its, "iov") <- TRUE
attr(nlmixr2Est.its, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
attr(nlmixr2Est.its, "nlmixr2Priors") <- "general"
