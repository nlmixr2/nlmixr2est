# est="impmap" -- importance-sampling EM.
#
# impmap centers a multivariate-normal proposal at each subject's MAP mode
# (produced by the mu-referenced FOCEI inner problem, muModel="lin"), draws
# importance samples, and updates the population parameters by an EM step.
# The numerical kernel lives in C++ (src/imp.cpp); this file is orchestration
# only: control construction, dispatch, and post-fit assembly.

# Importance-sampling / EM control names -- stripped when down-converting to a
# plain foceiControl for the MAP inner problem / output.
.impmapIsControlNames <- c("isample", "nIter", "mapIter", "gamma",
                           "iscaleMin", "iscaleMax", "iaccept",
                           "ctol", "nConvWindow", "impSeed")

#' Control options for the impmap (importance-sampling EM) estimation method
#'
#' A NONMEM-style Monte Carlo importance-sampling EM built on the mu-referenced
#' FOCEI MAP.  The proposal density for each subject is centered at the MAP mode
#' (`muModel="lin"`); mu-referenced population parameters are updated by the EM
#' gradient, while non-mu parameters (including the residual error) are updated
#' by finite differences.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param isample Number of importance samples drawn per subject per iteration
#'   (NONMEM ISAMPLE).
#' @param nIter Maximum number of importance-sampling EM iterations.
#' @param mapIter Number of MAP re-centering iterations per EM step; `> 0`
#'   re-centers the proposal at the MAP mode each iteration.
#' @param gamma Initial proposal-variance inflation factor (NONMEM ISCALE); the
#'   proposal covariance is `gamma` times the inverse of the inner information
#'   matrix at the mode.
#' @param iscaleMin,iscaleMax Lower/upper bounds for the adapted `gamma`
#'   (NONMEM ISCALE_MIN / ISCALE_MAX).
#' @param iaccept Target importance-sampling acceptance ratio used to adapt
#'   `gamma` (NONMEM IACCEPT).
#' @param ctol Convergence tolerance on the windowed objective-function change;
#'   `NULL` derives it from `sigdig`.
#' @param nConvWindow Length of the trailing iteration window used to average
#'   the objective-function change for convergence (NONMEM-style CTYPE).
#' @param muModel Mu-referencing variant for the MAP inner problem; for
#'   `impmapControl()` this is always `"lin"` and cannot be changed.
#' @param impSeed Base seed for the per-subject thread-safe (threefry) RNG
#'   streams; results are reproducible and independent of the thread count.
#' @return impmapControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' impmapControl()
impmapControl <- function(sigdig=3,
                          ...,
                          isample=300L,
                          nIter=100L,
                          mapIter=1L,
                          gamma=1.0,
                          iscaleMin=0.1,
                          iscaleMax=10.0,
                          iaccept=0.4,
                          ctol=NULL,
                          nConvWindow=10L,
                          impSeed=42L,
                          muModel=c("lin", "none")) {
  muModel <- match.arg(muModel)
  .control <- foceiControl(sigdig=sigdig, ..., muModel="lin")
  .control$isample <- as.integer(isample)
  .control$nIter <- as.integer(nIter)
  .control$mapIter <- as.integer(mapIter)
  .control$gamma <- as.double(gamma)
  .control$iscaleMin <- as.double(iscaleMin)
  .control$iscaleMax <- as.double(iscaleMax)
  .control$iaccept <- as.double(iaccept)
  .control$ctol <- if (is.null(ctol)) NULL else as.double(ctol)
  .control$nConvWindow <- as.integer(nConvWindow)
  .control$impSeed <- as.integer(impSeed)
  class(.control) <- "impmapControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.impmapControl <- function(control, env) {
  assign("impmapControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.impmap <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- impmapControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("impmapControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to impmapControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(impmapControl, .ctl)
  } else if (!inherits(.ctl, "impmapControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- impmapControl()
  } else {
    .ctl <- do.call(impmapControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.impmap <- function(x, ...) {
  .env <- x[[1]]
  if (exists("impmapControl", .env)) {
    .control <- get("impmapControl", .env)
    if (inherits(.control, "impmapControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "impmapControl")) return(.control)
  }
  stop("cannot find impmap related control object", call.=FALSE)
}

.impmapControlToFoceiControl <- function(env, assign=TRUE) {
  .impmapControl <- env$impmapControl
  .n <- setdiff(names(.impmapControl), .impmapIsControlNames)
  .foceiControl <- setNames(lapply(.n, function(n) .impmapControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.impmap <- function(x, ...) {
  .env <- x[[1]]
  .impmapControlToFoceiControl(.env, assign=FALSE)
}

#' Install the impmap control into the ui
#'
#' @param env Environment with ui in it
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.impmapFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- impmapControl()
  }
  if (!inherits(.control, "impmapControl")) {
    .control <- do.call(nlmixr2est::impmapControl, .control)
  }
  assign("control", .control, envir=.ui)
}

#' Fit the impmap family of models
#'
#' @param env Environment from nlmixr2Est
#' @param ui rxode2 ui object
#' @param ... Other arguments
#' @return fit environment
#' @author Matthew L. Fidler
#' @noRd
.impmapFamilyFit <- function(env, ui, ...) {
  # With est="impmap", foceiFitCpp_ runs impOuter() (src/imp.cpp) in place of
  # foceiOuter().  The FOCEI outer optimizer is turned off (maxOuterIterations=0)
  # because impOuter drives its own EM iteration; covariance is off for now.
  .control <- ui$control
  .control$maxOuterIterations <- 0L
  .control$covMethod <- 0L
  # 0-based index maps for the SIMPLE mu-referenced intercepts (theta = population
  # mean of an eta, no covariates): impOuter's M-step shifts each such theta by
  # the mean conditional eta.  Covariate mu-groups are excluded here because they
  # are handled by the regression update (updateMuGroups) instead.
  .env <- ui$foceiOptEnv  # builds foceiMuGroupTheta (the covariate-group thetas)
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta), ]
  .thNames <- .th[order(.th$ntheta), "name"]
  .etaRows <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, ]
  .etaNames <- .etaRows[order(.etaRows$neta1), "name"]
  .mr <- ui$muRefDataFrame
  .muThetaIdx <- as.integer(match(.mr$theta, .thNames) - 1L)
  .muEtaIdx <- as.integer(match(.mr$eta, .etaNames) - 1L)
  .covGroupTheta <- rxode2::rxGetControl(ui, "foceiMuGroupTheta", integer(0))
  .keep <- !is.na(.muThetaIdx) & !is.na(.muEtaIdx) &
    !(.muThetaIdx %in% .covGroupTheta)
  .control$impMuThetaIdx <- .muThetaIdx[.keep]
  .control$impMuEtaIdx <- .muEtaIdx[.keep]
  # 0-based theta indices of the estimated non-mu thetas (structural + residual
  # error) with sensitivity outputs in the sensitivity model; the M-step Newton
  # update maps its output columns back to these thetas.
  .control$impThetaSensIdx <- as.integer(.impmapEstTheta(ui)$all - 1L)
  assign("control", .control, envir=ui)
  .foceiFamilyReturn(env, ui, ..., est="impmap")
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.impmap <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'impmap'", .var.name=.ui$modelName)
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'impmap'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="impmapControl")
  on.exit({
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  env$impmapControl <- .control
  env$est <- "impmap"
  .ui <- env$ui
  .impmapFamilyFit(env, .ui, ...)
}
attr(nlmixr2Est.impmap, "covPresent") <- TRUE
attr(nlmixr2Est.impmap, "unbounded") <- .foUnbounded
attr(nlmixr2Est.impmap, "iov") <- TRUE
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook, R/mu2.R),
# gated on muModel/muRefCovAlg, exactly as the mufocei family does.
attr(nlmixr2Est.impmap, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
