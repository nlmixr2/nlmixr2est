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
                           "ctol", "nConvWindow", "impSeed", "impCov",
                           "qr", "qrShift", "qrRefresh", "sir", "sirSample",
                           # internal M-step index maps added in .impmapFamilyFit;
                           # not foceiControl() arguments, so they must be dropped
                           # when down-converting (e.g. .setOfvFo's do.call(foceiControl))
                           "impMuThetaIdx", "impMuEtaIdx", "impThetaSensIdx",
                           "impOmegaFixedEta")

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
#' @param iaccept Minimum importance-sampling effective-sample fraction
#'   (NONMEM IACCEPT).  The proposal scale `gamma` is kept at its efficient
#'   starting value while the achieved fraction stays at or above `iaccept`, and
#'   is inflated (toward `iscaleMax`) only when it drops below this floor.
#' @param ctol Convergence tolerance on the windowed objective-function change;
#'   `NULL` derives it from `sigdig`.
#' @param nConvWindow Length of the trailing iteration window used to average
#'   the objective-function change for convergence (NONMEM-style CTYPE).
#' @param muModel Mu-referencing variant for the MAP inner problem; for
#'   `impmapControl()` this is always `"lin"` and cannot be changed.
#' @param impSeed Base seed for the per-subject thread-safe (threefry) RNG
#'   streams; results are reproducible and independent of the thread count.
#' @param covMethod Covariance method.  `"imp"` (default) computes the
#'   Monte-Carlo importance-sampling observed-information covariance for the
#'   estimated thetas and Omega parameters (a finite-difference Hessian of the
#'   importance-sampling objective over fixed common-random-number samples),
#'   stashed as `$impCov` / `$impSe` and installed as the fit covariance; the
#'   theta standard errors match the Hessian-based FOCEI covariance, though the
#'   variance of a tightly-determined random effect (an Omega diagonal) can be
#'   over-estimated because the fixed samples barely span its prior variation.
#'   `"analytic"`, `"r,s"`, `"r"`, `"s"` instead compute the FOCEI covariance
#'   post-fit at the converged estimates (see [foceiControl()]); `""` skips the
#'   covariance step.
#' @param qr When `TRUE`, draw quasi-random (Sobol low-discrepancy) importance
#'   samples instead of pseudo-random Gaussian samples (QRPEM, Leary &
#'   Dunlavey PAGE 2012); the E-step integrals converge at O(1/N) instead of
#'   O(1/sqrt(N)).
#' @param qrShift Only used with `qr=TRUE`.  When `TRUE` each (iteration,
#'   subject) applies a random Cranley-Patterson shift to the Sobol points
#'   (seeded, thread-count independent); `FALSE` reuses one fixed Sobol point
#'   set everywhere (fully deterministic E-step, no RNG in the draw).
#' @param qrRefresh Only used with `qr=TRUE` and `qrShift=TRUE`.  When `TRUE`
#'   the shift is redrawn each iteration so residual quasi-random error
#'   averages out over the EM; `FALSE` draws one shift per subject at the fit
#'   start, making each EM iteration a deterministic map (smoothest objective
#'   trace).
#' @param sir When `TRUE`, accelerate the non-mu / residual-error M-step by
#'   SIR (sampling-importance-resampling): the theta-sensitivity Newton step
#'   uses `sirSample` equal-weight resampled points per subject instead of all
#'   `isample` weighted samples.
#' @param sirSample Number of SIR resampled points per subject; `NULL` uses
#'   `max(25, ceiling(isample/10))`.  Must be at most `isample`.
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
                          covMethod=c("imp", "analytic", "r,s", "r", "s", ""),
                          qr=FALSE,
                          qrShift=TRUE,
                          qrRefresh=TRUE,
                          sir=FALSE,
                          sirSample=NULL,
                          muModel=c("lin", "none")) {
  muModel <- match.arg(muModel)
  checkmate::assertLogical(qr, any.missing=FALSE, len=1, .var.name="qr")
  checkmate::assertLogical(qrShift, any.missing=FALSE, len=1, .var.name="qrShift")
  checkmate::assertLogical(qrRefresh, any.missing=FALSE, len=1, .var.name="qrRefresh")
  checkmate::assertLogical(sir, any.missing=FALSE, len=1, .var.name="sir")
  checkmate::assertIntegerish(isample, any.missing=FALSE, len=1, lower=1, .var.name="isample")
  checkmate::assertIntegerish(impSeed, any.missing=FALSE, len=1, .var.name="impSeed")
  .isample <- as.integer(isample)
  if (is.null(sirSample)) {
    .sirSample <- max(25L, as.integer(ceiling(.isample / 10)))
  } else {
    checkmate::assertIntegerish(sirSample, any.missing=FALSE, len=1, lower=1,
                                .var.name="sirSample")
    .sirSample <- as.integer(sirSample)
  }
  if (.sirSample > .isample) {
    stop("'sirSample' (", .sirSample, ") cannot exceed 'isample' (", .isample, ")",
         call.=FALSE)
  }
  # covMethod="imp" drives the Monte-Carlo importance-sampling covariance in the
  # C++ kernel (op_focei.impCov); every other token is the post-fit FOCEI
  # covariance, so hand foceiControl a valid token (the estimation pass forces
  # covMethod=0L regardless -- see .impmapFamilyFit).
  .dots <- list(...)
  .impCov <- isTRUE(.dots$impCov)   # may already be set on a round-tripped control
  .dots$impCov <- NULL              # internal field; do not forward to foceiControl
  if (is.character(covMethod)) {
    if (length(covMethod) == 1L && !nzchar(covMethod)) {
      covMethod <- ""
    } else {
      covMethod <- match.arg(covMethod)
    }
    .impCov <- identical(covMethod, "imp")
    .foceiCovMethod <- if (.impCov) "analytic" else covMethod
  } else {
    # round-trip: covMethod is already a resolved foceiControl integer slot;
    # keep the impCov flag from the incoming control (read above)
    .foceiCovMethod <- covMethod
  }
  .control <- do.call(foceiControl,
                      c(list(sigdig=sigdig), .dots,
                        list(covMethod=.foceiCovMethod, muModel="lin")))
  .control$impCov <- .impCov
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
  .control$qr <- qr
  .control$qrShift <- qrShift
  .control$qrRefresh <- qrRefresh
  .control$sir <- sir
  .control$sirSample <- .sirSample
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
  # np* internals (npBoxLower/npPoints/npResidFreeze ...) and the npag/npb user
  # knobs (points/cycles/gammaOptimize/... -- but NOT seed, a real foceiControl arg)
  # are nonparametric-engine control fields that foceiControl does not accept; drop
  # them so a downstream do.call(foceiControl, .) (e.g. .setOfvFo, general-likelihood
  # tables) does not error with "unused argument".
  .npKnobs <- c("points", "cycles", "gammaOptimize", "residOptimize", "muExpand",
                "gridWidth", "gridBounds",
                "alpha", "burnin", "nsamp", "nchains", "propSd", "est")
  .n <- .n[!grepl("^np[A-Z]", .n) & !(.n %in% .npKnobs)]
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
  # because impOuter drives its own EM iteration; the in-fit covariance is off
  # (it would bail on muModel="lin" anyway) -- the covariance is computed
  # post-fit on the full model at the converged estimates (.foceiRecomputeMuCov).
  .control <- ui$control
  .covMethodUser <- .control$covMethod  # restored on the fit env control below
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
  # 0-based eta indices whose Omega diagonal is FIXED; the EM Omega update
  # restores their rows/columns to the starting value so fix()ed variances hold.
  .etaOrd <- .etaRows[order(.etaRows$neta1), ]
  .control$impOmegaFixedEta <- as.integer(which(isTRUE(.etaOrd$fix) | .etaOrd$fix) - 1L)
  assign("control", .control, envir=ui)
  # Seed the importance-sampling RNG from the control (impSeed) right before the
  # fit, mirroring saem's set.seed(seed).  The E-step draws through rxode2's
  # threefry engine (getRxSeed1), so a fixed rxseed makes the fit reproducible
  # and independent of the thread count / ambient RNG state; rxWithSeed restores
  # the prior seed state afterward so it does not leak globally.
  .impSeed <- if (is.null(.control$impSeed)) 42L else as.integer(.control$impSeed)
  # est is "impmap" or "imp" (the no-MAP-search variant); pass it through so the
  # C++ kernel (impOuter) selects the proposal accordingly.
  .est <- if (exists("est", envir=env)) get("est", envir=env) else "impmap"
  .fit <- rxode2::rxWithSeed(.impSeed, rxseed=.impSeed,
                             code=.foceiFamilyReturn(env, ui, ..., est=.est))
  # The MC covariance (impCov=TRUE) is published with theta row/column names but
  # the Omega parameters come out unnamed on this path; fill them in (defensively,
  # only when the counts line up) so vcov()/$cov and the correlation are labelled.
  .impmapNameCov(.fit, ui)
  .impRestoreCovMethod(.fit, .covMethodUser)
  .fit
}

#' Restore the user's covMethod on the fit env's stored control
#'
#' The estimation pass forces covMethod=0L (the in-fit C++ step would bail on
#' muModel="lin"), and that runtime control is what gets stored on the fit env;
#' the post-fit recompute (.foceiRecomputeMuCov) reads the covMethod from there,
#' so put the user's choice back.
#' @noRd
.impRestoreCovMethod <- function(fit, covMethod) {
  .fenv <- tryCatch(fit$env, error = function(e) NULL)
  if (is.environment(.fenv) &&
        exists("impmapControl", envir = .fenv, inherits = FALSE)) {
    .ic <- get("impmapControl", envir = .fenv)
    .ic$covMethod <- covMethod
    assign("impmapControl", .ic, envir = .fenv)
  }
  invisible(fit)
}

#' Fill the Omega row/column names on the impmap covariance
#' @param fit impmap fit
#' @param ui rxode2 ui
#' @return Nothing, called for side effects
#' @noRd
.impmapNameCov <- function(fit, ui) {
  .fenv <- tryCatch(fit$env, error=function(e) NULL)
  if (is.null(.fenv) || is.null(.fenv$cov) || !is.matrix(.fenv$cov)) return(invisible())
  tryCatch({
    .dn <- dimnames(.fenv$cov)[[1]]
    if (is.null(.dn)) return(invisible())
    .empty <- which(is.na(.dn) | .dn == "")
    if (length(.empty) == 0L) return(invisible())
    .etaN <- .foceiEtaThetaMap(ui)$etaNames
    .op <- .foceiOmegaPairs(.fenv$omega, ui$iniDf)
    .omN <- .foceiOmegaCovNames(.op, .etaN)
    if (length(.omN) == length(.empty)) {
      .dn[.empty] <- .omN
      dimnames(.fenv$cov) <- list(.dn, .dn)
      if (!is.null(.fenv$fullCor) && is.matrix(.fenv$fullCor)) {
        dimnames(.fenv$fullCor) <- list(.dn, .dn)
      }
    }
  }, error=function(e) NULL)
  invisible()
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.impmap <- function(env, ...) {
  .ui <- env$ui
  # General (dnorm/ll) likelihoods flow through the shared FOCEI inner problem
  # (impEvalJointLik = likInner0), so only require transformable normality when
  # the rxode2 build has no llik support -- mirrors nlmixr2Est.focei.
  if (!rxode2hasLlik()) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'impmap'", .var.name=.ui$modelName)
  }
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
# gated on muModel/muRefCovAlg, exactly as the mfocei family does.
attr(nlmixr2Est.impmap, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
