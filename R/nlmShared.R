# Is `x` a usable Shi (2021) epsilon (a single finite value > 0)?
.isPositiveErr <- function(x) {
  checkmate::testNumber(x, finite = TRUE) && x > 0
}

# Assert a finite-difference epsilon (Shi 2021 `shiErr`/`hessErr`, Gill 1983
# `hessEps`) is strictly positive; 0 makes the searched step 0.
.assertPositiveEps <- function(x, .var.name = checkmate::vname(x), null.ok = FALSE) {
  checkmate::assertNumber(x, finite = TRUE, null.ok = null.ok, .var.name = .var.name)
  if (!is.null(x) && x <= 0) {
    stop("'", .var.name, "' must be > 0", call. = FALSE)
  }
  invisible(x)
}

#' Setup a nonlinear system for optimization
#'
#' @param par A named vector of initial estimates to setup the
#'   nonlinear model solving environment. The names of the parameter
#'   should match the names of the model to run (not `THETA[#]` as
#'   required in the `modelInfo` argument)
#'
#' @param ui rxode2 ui model
#'
#' @param data rxode2 compatible data for solving/setting up
#'
#' @param modelInfo A list with `predOnly` (predictions-only model in terms of
#'   `THETA[#]`/`DV`), `eventTheta` (0/1 per THETA flagging event-related
#'   parameters that need Shi2021 finite differences, same length as `par`),
#'   and `thetaGrad` (needed when solveType != 1; gives value/gradient per
#'   THETA). See `ui$nlmSensModel` or `ui$nlmRxModel` for examples.
#' @param control control structure; required: `rxControl`, `stickyRecalcN`,
#'   `maxOdeRecalc`, `odeRecalcFactor`. Optional: `solveType`, `eventType`,
#'   `shi21maxFD`, `shiErr`, `optimHessType`, `shi21maxHess`, `hessErr`,
#'   `useColor`, `printNcol`, `print`, `normType`, `scaleType`, `scaleCmin`,
#'   `scaleCmax`, `scaleTo`, `scaleC`, `gradTo` (default 0 if missing).
#' @param lower lower bounds, will be scaled if present
#' @param upper upper bounds, will be scaled if present
#' @return nlm solve environment; key fields: `$par.ini`, `$lower`, `$upper`
#'   (all scaled), and `$.ctl` (control structure).
#'
#' @details
#'
#' No rxode2 solving should occur between setup calls; prints the solving header if `print != 0`.
#'
#' @author Matthew Fidler
#' @keywords internal
#'
#' @export
.nlmSetupEnv <- function(par, ui, data, modelInfo, control, lower = NULL, upper = NULL) {
  .ctl <- control
  if (!any(names(.ctl) == "gradTo")) {
    .ctl$gradTo <- 0.0
  }
  if (!any(names(.ctl) == "solveType")) {
    if (any(names(modelInfo) == "thetaGrad")) {
      .ctl$solveType <- 2L
    } else {
      .ctl$solveType <- 1L
    }
  }
  if (!any(names(.ctl) == "eventType")) {
    .ctl$eventType <- 1L
  }
  if (!any(names(.ctl) == "optimHessType")) {
    .ctl$optimHessType <- 1L
  }
  if (!any(names(.ctl) == "shi21maxFD")) {
    .ctl$shi21maxFD <- 20L
  }
  if (!any(names(.ctl) == "shi21maxHess")) {
    .ctl$shi21maxHess <- 20L
  }
  # a non-positive Shi (2021) epsilon degenerates the step search to a 0 step
  if (!.isPositiveErr(.ctl$shiErr)) {
    .ctl$shiErr <- (.Machine$double.eps)^(1 / 3)
  }
  if (!.isPositiveErr(.ctl$hessErr)) {
    .ctl$hessErr <- (.Machine$double.eps)^(1 / 3)
  }
  # nlmSetup (nlm.cpp) reads control$iterPrintControl; external callers that
  # hand-build a control (e.g. babelmixr2 fmeMcmc) may omit it, so synthesize
  # one from the scalar print args instead of erroring in C.
  if (!inherits(.ctl$iterPrintControl, "iterPrintControl")) {
    .ctl$iterPrintControl <-
      .absorbIterPrintControl(
        print = if (is.null(.ctl$print)) 1L else .ctl$print,
        printNcol = .ctl$printNcol,
        useColor = .ctl$useColor
      )
  }

  .env <- new.env(parent = emptyenv())
  .env$rxControl <- .ctl$rxControl
  .env$thetaNames <- names(par)
  .f <- modelInfo
  if (any(names(.f) == "thetaGrad")) {
    .env$predOnly <- .f$predOnly
    nlmixr2global$nlmEnv$model <- .env$thetaGrad <- .f$thetaGrad
  } else {
    nlmixr2global$nlmEnv$model <- .env$predOnly <- .f$predOnly
  }
  ## Delay differential equation models need a dense-output solver so delay()
  ## history is recorded and interpolated (also by the forward-sensitivity /
  ## jump-sensitivity states).  nlm.cpp calls rxSolve_ directly, bypassing
  ## rxSolve()'s hasDelay enforcement, so replicate it here: dense dop853 (which
  ## needs no analytic Jacobian) unless an in-engine sensitivity method (>=200)
  ## is already selected -- those record their own dense history.
  if (isTRUE(rxode2::rxModelVars(nlmixr2global$nlmEnv$model)$flags[["hasDelay"]] == 1L)) {
    .env$rxControl$dense <- TRUE
    if (is.null(.env$rxControl$method) || .env$rxControl$method < 200L) {
      .env$rxControl$method <- 0L
      .env$rxControl$stiff2 <- 0L
    }
  }
  # matExp-native sensitivities (#860): same reasoning as the hasDelay guard
  # above -- nlm.cpp calls rxSolve_() directly, bypassing rxSolve.default()'s
  # S3-dispatched matExp()/indLin() auto-detection that would otherwise force
  # method="indLin".  Left at the ordinary-ODE default, a matExp model's
  # dydt() is a no-op stub and every solve silently free-runs a zero
  # derivative (states frozen at their post-dose values).
  .mv <- rxode2::rxModelVars(nlmixr2global$nlmEnv$model)
  if (is.list(.mv$indLin) && length(.mv$indLin) == 4L && !isTRUE(unname(.env$rxControl$method) == 3L)) {
    .env$rxControl$method <- 3L
  }
  .env$param <- setNames(par, sprintf("THETA[%d]", seq_along(par)))
  .nlmFitDataSetup(data)
  .env$needFD <- .f$eventTheta
  # Ship a single xform sub-list so C wires log/logit/probit back-transforms
  # via scaleAttachXform (src/scale.h), like focei/saem; nlm-family prints thetas in `par` order.
  .ctl$xform <- .iterPrintXParFromUi(ui, names(par))
  .env$control <- .ctl
  .env$data <- nlmixr2global$nlmEnv$data
  .Call(`_nlmixr2est_nlmSetup`, .env)
  ## Activate event-jump sensitivity (if eventSens="jump") now, before the
  ## scaleC solve -- doing it after would mis-scale dosing params from a jump-less gradient.
  ## Deactivated in .nlmFreeEnv; no-op for "fd".
  .env$esActive <- FALSE
  if (!is.null(.env$thetaGrad)) {
    .env$esActive <- isTRUE(tryCatch(
      rxode2::rxEventSensLoadModel(.env$thetaGrad),
      error = function(e) FALSE
    ))
  }
  if (is.null(.ctl$scaleC)) {
    .ctl$scaleC <- ui$scaleCtheta
    if (.ctl$scaleType == 2L && .ctl$gradTo > 0) {
      # the derivative-based |gradTo/gradient_i(par)| blows up for a near-zero
      # starting gradient (issue #994), so it gets FOCEi's band guard, falling
      # back to the transform-aware ui$scaleCtheta
      .gradScaleC <- .Call(`_nlmixr2est_nlmGetScaleC`, par, .ctl$gradTo)
      if (length(.gradScaleC) > 0L) {
        .ctl$scaleC <- mapply(.guardScaleC, .gradScaleC, .ctl$scaleC, USE.NAMES = FALSE)
      }
    }
  }
  # also replaces the unguarded values nlmGetScaleC() leaves in the C++ buffer
  .Call(`_nlmixr2est_nlmSetScaleC`, .ctl$scaleC)
  .env$scaleC <- .ctl$scaleC
  .p <- .Call(`_nlmixr2est_nlmScalePar`, par)
  .env$par.ini <- .p
  .env$.ctl <- .ctl
  .env$upper <- .env$lower <- NULL
  if (!is.null(upper)) {
    .env$upper <- .Call(`_nlmixr2est_nlmScalePar`, upper)
  }
  if (is.null(.env$upper)) {
    .env$upper <- rep(Inf, length(.p))
  }
  if (!is.null(lower)) {
    .env$lower <- .Call(`_nlmixr2est_nlmScalePar`, lower)
  }
  if (is.null(.env$lower)) {
    .env$lower <- rep(-Inf, length(.p))
  }
  .Call(`_nlmixr2est_nlmPrintHeader`)
  .env
}

#' Frees nlm environment
#'
#' @return Nothing, called for side effects
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.nlmFreeEnv <- function() {
  .Call(`_nlmixr2est_nlmFree`)
  rxode2::rxSolveFree()
  ## Deactivate any event-jump sensitivity injection from .nlmSetupEnv; no-op if never activated.
  tryCatch(rxode2::rxEventSensDeactivate(), error = function(e) NULL)
}
#' Finalizes output list
#'
#' @param env nlm environment
#' @param lst output list
#' @param par parameter name of final estimate in output
#' @param printLine Print the final line when print is nonzero
#' @param hessianCov boolean indicating a hessian should be
#'   used/calculated for covariance
#' @return modified list with `$cov`
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.nlmFinalizeList <- function(env, lst, par = "par", printLine = TRUE, hessianCov = TRUE) {
  .ret <- lst
  .ctl <- env$.ctl
  .ret$scaleC <- env$scaleC
  .ret$parHistData <- .Call(`_nlmixr2est_nlmGetParHist`, printLine)
  .name <- env$thetaNames
  if (inherits(lst, "nls")) {
    .ret[[par]] <- coef(lst)
  }
  .parScaled <- setNames(.ret[[par]], .name)
  .ret[[paste0(par, ".scaled")]] <- .parScaled
  .par <- .Call(`_nlmixr2est_nlmUnscalePar`, .parScaled)
  .ret[[par]] <- setNames(.par, .name)
  # if using hessian to caluclate covariance
  if (!any(names(.ctl) == "covMethod")) {
    .ctl$covMethod <- "r"
  }
  # the residual degrees of freedom of a least-squares (nls) fit, NA otherwise
  .rdf <- if (inherits(lst, "nls")) {
    length(stats::residuals(lst)) - length(.parScaled)
  } else if (inherits(lst, "nls.lm")) {
    length(lst$fvec) - length(lst$par)
  } else {
    NA_integer_
  }
  if (!is.na(.rdf) && .rdf <= 0 && (inherits(lst, "nls") || (hessianCov && .ctl$covMethod != ""))) {
    # sigma^2 = RSS / (n - p) does not exist
    warning(
      sprintf(
        "nls has %d residual degrees of freedom, no residual variance; covariance step failed",
        as.integer(.rdf)
      ),
      call. = FALSE
    )
    .ret$covMethod <- "failed"
  } else if (inherits(lst, "nls")) {
    # sigma^2 (J'J)^-1, the residual variance times summary()$cov.unscaled
    .g <- .covGuard(stats::vcov(lst))
    if (.g$ok) {
      .ret$cov.scaled <- .g$cov
      .ret$cov <- .Call(`_nlmixr2est_nlmAdjustCov`, .ret$cov.scaled, .parScaled)
      .ret$covMethod <- "nls"
    } else {
      warning(sprintf("nls covariance %s; covariance step failed", .g$reason), call. = FALSE)
      .ret$covMethod <- "failed"
    }
  } else if (hessianCov && .ctl$covMethod != "") {
    .malert("calculating covariance")
    if (!any(names(.ret) == "hessian")) {
      .p <- setNames(.parScaled, NULL)
      # the stencil differences objectives that agree in their last digits, so its
      # solves run at the covariance probe tolerances, as FOCEi's do
      .tol <- .Call(`_nlmixr2est_covProbeSolveTolSet_`)
      .hess <- tryCatch(
        nlmixr2Hess(.p, nlmixr2est::.nlmixrNlmFunC),
        finally = .Call(`_nlmixr2est_covProbeSolveTolRestore_`, .tol)
      )
      .ret$hessian <- .hess
    }
    dimnames(.ret$hessian) <- list(.name, .name)
    # r matrix: nlm-family objectives (`.nlmixrNlmFunC`/optim's `fn`) are
    # built as a plain -1*LL, not the -2*LL scale FOCEI/SAEM/etc use -- so
    # the Hessian is already the Fisher information (unlike the FOCEI R matrix,
    # which halves a -2*LL Hessian to get there).  Do not rescale here.
    .r <- .ret$hessian
    if (inherits(lst, "nls.lm")) {
      # minpack.lm's hessian is J'J of the residuals; the information is
      # J'J / sigma^2, sigma^2 = RSS / (n - p) as in its vcov.nls.lm()
      .r <- .r / (lst$deviance / (length(lst$fvec) - length(lst$par)))
    }
    .rc <- .nlmCovFromHessian(.r)
    if (!is.null(.rc$warning)) {
      warning(.rc$warning, call. = FALSE)
    }
    if (is.null(.rc$r)) {
      .ret$covMethod <- "failed"
    } else {
      .rinv <- rxode2::rxInv(.rc$u)
      .cov <- .rinv %*% t(.rinv)
      dimnames(.cov) <- list(.name, .name)
      .ret$covMethod <- if (.ctl$covMethod != "r") paste0(.rc$type, " (", .ctl$covMethod, ")") else .rc$type
      .ret$cov.scaled <- .cov
      .ret$cov <- .Call(`_nlmixr2est_nlmAdjustCov`, .ret$cov.scaled, .parScaled)
    }
    .ret$r <- .r
    .msuccess("done")
  }
  .ret$censInformation <- .Call(`_nlmixr2est_nlmCensInfo`)
  .Call(`_nlmixr2est_nlmWarnings`)
  .nlmFreeEnv()
  .ret
}

#' The nlm-family Hessian from central differences of the analytic gradient
#'
#' Column `k` is the central difference of `.nlmixrNlminbGradC()` over the
#' central step `nlmixr2Gill83()` finds for parameter `k` (the step source
#' `nlmixr2Hess()` uses), at the covariance probe tolerances.  It needs a loaded
#' problem with a gradient solve type (`solveType` `"grad"` or `"hessian"`).
#' @param par scaled parameters at the estimates
#' @return the symmetric Hessian of the objective
#' @noRd
.nlmGradHessian <- function(par) {
  .tol <- .Call(`_nlmixr2est_covProbeSolveTolSet_`)
  on.exit(.Call(`_nlmixr2est_covProbeSolveTolRestore_`, .tol))
  .h <- nlmixr2Gill83(nlmixr2est::.nlmixrNlminbFunC, par)$aEpsC
  .n <- length(par)
  .hess <- matrix(0, .n, .n)
  for (.k in seq_len(.n)) {
    .e <- replace(numeric(.n), .k, .h[.k])
    .hess[, .k] <- (nlmixr2est::.nlmixrNlminbGradC(par + .e) -
      nlmixr2est::.nlmixrNlminbGradC(par - .e)) / (2 * .h[.k])
  }
  (.hess + t(.hess)) / 2
}
#' The information matrix an nlm-family covariance is inverted from
#'
#' The Hessian is factored and, when needed, repaired by the covariance step's one
#' acceptance rule (`covAccept_()`, `covAcceptRule()` in src/cholse.cpp), under the
#' FOCEi labels: "r" as it is; "r+" when Schnabel-Eskow's modified Cholesky
#' factorization adds diagonals within `foceiControl()`'s default `cholAccept`;
#' "|r|" for `sqrtm(R %*% R)`.  A numerically rank-deficient R is not repaired.
#' @param hess Hessian of the -LL objective (the R matrix)
#' @return list(r = the (repaired) R matrix and u = its upper Cholesky factor,
#'   both `NULL` when none is usable; type = "r", "r+", "|r|" or "failed";
#'   warning = what was done, `NULL` for "r")
#' @noRd
.nlmCovFromHessian <- function(hess) {
  .g <- .covGuard(hess)
  if (is.null(.g$cov)) {
    return(list(type = "failed", warning = sprintf("R matrix %s; covariance step failed", .g$reason)))
  }
  .r <- .g$cov
  .a <- covAccept_(.r, (.Machine$double.eps)^(1 / 3), formals(foceiControl)$cholAccept, TRUE)
  .reason <- if (.g$ok) "is nearly singular" else .g$reason
  switch(
    .a$type,
    "+" = list(r = .r, u = .a$U, type = "r+", warning = sprintf("R matrix %s; corrected as \"r+\"", .reason)),
    "|" = list(r = .a$M, u = .a$U, type = "|r|", warning = sprintf("R matrix %s; corrected as \"|r|\"", .reason)),
    singular = list(type = "failed", warning = "R matrix is singular; covariance step failed"),
    failed = list(type = "failed", warning = sprintf("R matrix %s; covariance step failed", .reason)),
    list(r = hess, u = .a$U, type = "r")
  )
}

#' Adjust nlm and family output environment
#'
#' Will take information like `$censInformation`, `$parHistData`,
#' `$cov` and `$covMethod` from the ret[[str]] and put it directly in
#' the environment `ret`
#'
#' @param ret environment for fit output that needs to be adjusted
#' @param str string for the fit output
#' @return updated environment
#' @keywords internal
#' @export
#' @author Matthew L. Fidler
.nlmFamilyAdjustOutput <- function(ret, str) {
  .nlm <- ret[[str]]
  .censInformation <- ret$censInformation
  if (
    is.null(.censInformation) &&
      !is.null(.nlm$censInformation)
  ) {
    .censInformation <- .nlm$censInformation
    ret[[str]][["censInformation"]] <- NULL
  }
  ret$censInformation <- .censInformation

  .parHistData <- ret$parHistData
  if (
    is.null(.parHistData) &&
      !is.null(.nlm$parHistData)
  ) {
    .parHistData <- .nlm$parHistData
    ret[[str]][["parHistData"]] <- NULL
  }
  ret$parHistData <- .parHistData

  .cov <- ret$cov
  if (
    is.null(.cov) &&
      !is.null(.nlm$cov)
  ) {
    .cov <- .nlm$cov
    ret[[str]][["cov"]] <- NULL
  }
  ret$cov <- .cov

  .covMethod <- ret$covMethod
  if (
    is.null(.covMethod) &&
      !is.null(.nlm$covMethod)
  ) {
    .covMethod <- .nlm$covMethod
    ret[[str]][["covMethod"]] <- NULL
  }
  ret$covMethod <- .covMethod

  ret
}

#' Adjust covariance matrix based on scaling parameters
#'
#' @param cov Covariance of scaled parameters
#' @param parScaled The final scaled parameter value
#' @return The adjusted covariance matrix based on the scaling
#' @export
#' @keywords internal
#' @author Matthew L. Fidler
.nlmAdjustCov <- function(cov, parScaled) {
  .Call(`_nlmixr2est_nlmAdjustCov`, cov, parScaled)
}

#' Uppercase data column names except the model covariates
#'
#' @param nms character vector of column names
#' @param covNames character vector of model covariate names (kept as-is)
#' @return character vector of names, upper-cased except those in `covNames`
#' @author Matthew L. Fidler
#' @noRd
.nmUpcaseNonCov <- function(nms, covNames) {
  if (is.null(covNames)) {
    covNames <- character(0)
  }
  vapply(
    nms,
    function(.x) {
      if (.x %in% covNames) .x else toupper(.x)
    },
    character(1),
    USE.NAMES = FALSE
  )
}

#' Detect the time-varying covariate columns for mu-referenced estimators (SAEM/NLME)
#'
#' @param dataSav preprocessed event-table data (from `.foceiPreProcessData()`)
#' @param ui rxode2 ui model (uses `ui$mv0`)
#' @param rxControl rxode2 control (for `addlKeepsCov`/`addlDropSs`/`ssAtDoseTime`)
#' @return character vector of time-varying covariate column names (possibly empty)
#' @author Matthew L. Fidler
#' @noRd
.nlmixrTimeVaryingCovariates <- function(dataSav, ui, rxControl) {
  .et <- rxode2::etTrans(
    dataSav,
    ui$mv0,
    addCmt = TRUE,
    addlKeepsCov = rxControl$addlKeepsCov,
    addlDropSs = rxControl$addlDropSs,
    ssAtDoseTime = rxControl$ssAtDoseTime
  )
  .nTv <- attr(class(.et), ".rxode2.lst")$nTv
  # nTv == 0 means no time-varying covariates; otherwise they follow the first 6 columns
  if (!is.null(.nTv) && .nTv == 0L) {
    return(character(0))
  }
  names(.et)[-seq_len(6)]
}

#' Stage the mu-referenced covariate split (time-varying vs not) into the ui env
#'
#' Splits the mu-referenced covariates: non-time-varying ones stay in
#' `muRefFinal` so they are absorbed into the phi term by the mu-ref drop
#' (`.saemDropMuRefFromModel` -> `$saemModel0` collapses the model to `phi +
#' timeVaryingCovariate*beta_cov`); time-varying ones are removed from
#' `muRefFinal` so they remain in the model as `beta_cov` regressors.  Both
#' `muRefFinal` and `timeVaryingCovariates` are assigned into the ui env and MUST
#' be removed on exit with `.nlmixrRmMuRefTimeVarying()`.  Shared by saem, the
#' mu-referenced focei family and vae so they all detect time-varying covariates
#' the same way.  Only the time-varying *split* is shared here; the model
#' expansion differs by method -- saem collapses lone etas into phi (theta forced
#' to 0), while vae and mu-referenced focei keep the etas as the inner problem
#' needs them.
#'
#' @param ui rxode2 ui (an environment) to stage the split into
#' @param timeVaryingCovariates character vector from
#'   `.nlmixrTimeVaryingCovariates()`
#' @return `ui`, invisibly (called for the side-effect assignments)
#' @noRd
.nlmixrSetMuRefTimeVarying <- function(ui, timeVaryingCovariates) {
  .muRefCovariateDataFrame <- ui$muRefCovariateDataFrame
  if (length(timeVaryingCovariates) > 0) {
    # A log-scale (exp-transformed) time-varying mu covariate can fit better
    # untransformed; the historical warning is left disabled but the detection
    # is kept so the behavior is easy to restore.
    .w <- which(.muRefCovariateDataFrame$covariate %in% timeVaryingCovariates)
    .covPar <- .muRefCovariateDataFrame[.w, "theta"]
    .w2 <- which(ui$muRefCurEval$parameter %in% .covPar)
    if (length(.w2) > 0) {
      .w3 <- which("exp" == ui$muRefCurEval$curEval[.w2])
      if (length(.w3) > 0) {
        .w2 <- .w2[.w3]
        .texp <- ui$muRefCurEval$parameter[.w2]
        .pars <- .muRefCovariateDataFrame$covariateParameter[.muRefCovariateDataFrame$theta %in% .texp]
        ## warning(paste0("log-scale mu referenced time varying covariates (",
        ##                paste(.pars, collapse=", "), ") may have better results ...
      }
    }
    # keep only non-time-varying covariates in the absorbed (mu-ref) set
    .muRefCovariateDataFrame <-
      .muRefCovariateDataFrame[!(.muRefCovariateDataFrame$covariate %in% timeVaryingCovariates), ]
  }
  assign("muRefFinal", .muRefCovariateDataFrame, ui)
  assign("timeVaryingCovariates", timeVaryingCovariates, ui)
  invisible(ui)
}

#' Remove the staged mu-ref time-varying covariate info from the ui env
#'
#' Undoes `.nlmixrSetMuRefTimeVarying()`; call from the estimation method's
#' `on.exit()` so the shared ui object is left unmodified after the fit.
#' @noRd
.nlmixrRmMuRefTimeVarying <- function(ui) {
  if (is.environment(ui) && exists("muRefFinal", envir = ui, inherits = FALSE)) {
    rm(list = "muRefFinal", envir = ui)
  }
  if (is.environment(ui) && exists("timeVaryingCovariates", envir = ui, inherits = FALSE)) {
    rm(list = "timeVaryingCovariates", envir = ui)
  }
  invisible(ui)
}

#' Integer code of an nlm-family control option given as a name or a code
#'
#' @param value the option as given: one of the names of `idx` (matched as
#'   `match.arg()` does, the first being the default) or one of its codes
#' @param idx the name -> code map, its names in the order of the option's
#'   choices
#' @param name the option's name, for the error
#' @return the integer code
#' @noRd
.nlmCtlCode <- function(value, idx, name) {
  if (!is.numeric(value)) {
    return(setNames(idx[match.arg(value, names(idx))], NULL))
  }
  if (length(value) != 1L || is.na(value) || !(value %in% idx)) {
    stop(
      "'",
      name,
      "' must be one of ",
      paste0(sprintf("\"%s\" (%d)", names(idx), idx), collapse = ", "),
      call. = FALSE
    )
  }
  as.integer(value)
}

#' The covMethod of an nlm-family control
#'
#' `""` (no covariance step) is one of the choices, which `match.arg()` cannot
#' match.
#' @param covMethod the argument as given
#' @param choice `match.arg(covMethod)` in the calling control; it is a promise,
#'   forced only when `covMethod` is not `""`
#' @return the name, or `""`
#' @noRd
.nlmCtlCovMethod <- function(covMethod, choice) {
  if (identical(covMethod, "")) "" else choice
}

#' Shared control setup for the nlm-family estimation methods
#'
#' @param env dispatch environment (provides `ui` and `control`)
#' @param controlFn the method's `*Control()` constructor (e.g. `nlmControl`)
#' @param controlClass the control object's S3 class (e.g. `"nlmControl"`)
#' @return Nothing; assigns the resolved control onto `env$ui`
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
.nlmFamilyControlGeneric <- function(env, controlFn, controlClass) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- controlFn()
  }
  if (!inherits(.control, controlClass)) {
    .control <- do.call(controlFn, .control)
  }
  assign("control", .control, envir = .ui)
}

#' The foceiControl that finalizes an nlm-family fit
#'
#' The optimizer has already run, so the FOCEi pass only builds the tables: no
#' outer or inner iterations, no covariance step, no interaction and no
#' scaling.  The settings that shape the model and the tables come from the
#' method's control; one it does not have (`sensMethod`, say) keeps the
#' `foceiControl()` default.
#' @param env fit environment holding the method's control
#' @param ctl name of the control in `env` (e.g. `"nlmControl"`)
#' @param assign when `TRUE`, also store the result as `env$control`
#' @param literalFixRes `literalFixRes` of the finalization (the control's own
#'   by default)
#' @return the `foceiControl()` object
#' @noRd
.nlmFamilyControlToFoceiControl <- function(env, ctl, assign = TRUE, literalFixRes = env[[ctl]]$literalFixRes) {
  .ctl <- env[[ctl]]
  .ret <- foceiControl(
    rxControl = .ctl$rxControl,
    maxOuterIterations = 0L,
    maxInnerIterations = 0L,
    covMethod = 0L,
    sumProd = .ctl$sumProd,
    optExpression = .ctl$optExpression,
    literalFix = .ctl$literalFix,
    literalFixRes = literalFixRes,
    scaleTo = 0,
    calcTables = .ctl$calcTables,
    addProp = .ctl$addProp,
    interaction = 0L,
    compress = .ctl$compress,
    ci = .ctl$ci,
    sigdigTable = .ctl$sigdigTable,
    indTolRelax = .ctl$indTolRelax,
    eventSens = .ctl$eventSens,
    sensMethod = .ctl$sensMethod
  )
  if (assign) {
    env$control <- .ret
  }
  .ret
}

#' The full theta vector of an nlm-family fit
#'
#' A fixed theta keeps its `ini()` value; an estimated one comes from the
#' optimizer's estimates, which are named by theta.
#' @param fit the optimizer result, as `.nlmFinalizeList()` returns it
#' @param ui rxode2 ui
#' @param par name of the estimates in `fit`
#' @return the thetas named and ordered as `ui$iniDf`
#' @noRd
.nlmFamilyGetTheta <- function(fit, ui, par) {
  .iniDf <- ui$iniDf
  .est <- fit[[par]]
  setNames(
    vapply(
      seq_along(.iniDf$name),
      function(i) {
        if (.iniDf$fix[i]) .iniDf$est[i] else .est[.iniDf$name[i]]
      },
      double(1),
      USE.NAMES = FALSE
    ),
    .iniDf$name
  )
}

#' Shared fit driver for the nlm-family estimation methods
#'
#' @param env dispatch environment (provides `ui`, `control`, `data`, `table`)
#' @param method estimation-method string; also the slot the raw fit is stored
#'   under (e.g. `"nlm"` -> `.ret[["nlm"]]`)
#' @param fitModel `function(ui, dataSav)` running the optimizer
#' @param getTheta `function(fit, ui)` returning the full theta vector, or the
#'   name of the optimizer's estimates in the fit (e.g. `"par"`), which
#'   `.nlmFamilyGetTheta()` completes with the fixed thetas
#' @param controlToFocei `function(env)` translating the control to a
#'   focei-style control for output assembly
#' @param returnFlag rxode2 control flag name that short-circuits and returns the
#'   raw optimizer result (e.g. `"returnNlm"`)
#' @param message `function(fit)` returning the `$message` (default `fit$message`)
#' @param emitFitWarnings when TRUE (the default), re-emit the warnings
#'   collected from `fitModel` (the optimizer, the covariance step and
#'   `nlmWarnings()`) via `warning()`, so they reach the fit's `$runInfo`;
#'   `FALSE` drops them
#' @param extra `$extra` print string, or a `function(control)` returning it
#' @param adjustOutput when TRUE, run `.nlmFamilyAdjustOutput()`
#' @param objective optional `function(fit)` returning the raw objective, or
#'   the name of the fit's minimized -log-likelihood, which is doubled; when
#'   `NULL` the driver does not set `$objective` (a `postSetup` closure did)
#' @param postSetup optional `function(ret, ui, fitList)` returning a modified
#'   `ret`, run right after the raw fit is stored and before
#'   `.nlmFamilyAdjustOutput()` (for methods that set cov/covMethod/objective
#'   with custom values)
#' @return the assembled nlmixr2 fit (or the raw optimizer result if `returnFlag`)
#' @author Matthew L. Fidler
#' @export
.nlmFamilyFitGeneric <- function(
  env,
  method,
  fitModel,
  getTheta,
  controlToFocei,
  returnFlag,
  objective = NULL,
  message = function(fit) fit$message,
  emitFitWarnings = TRUE,
  extra = "",
  adjustOutput = TRUE,
  postSetup = NULL
) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table
  nlmixrWithTiming("setup", {
    .foceiPreProcessData(.data, .ret, .ui, .control$rxControl)
  })
  # fitModel builds the symengine/sensitivity model and runs the optimizer;
  # time it as "optimize" so the work is not left in the "other" bucket (the
  # nlm-family model build and iterative solve are intertwined -- the
  # sensitivity model is the optimization model).
  .fit <- nlmixrWithTiming("optimize", {
    .collectWarn(fitModel(.ui, .ret$dataSav), lst = TRUE)
  })
  .ret[[method]] <- .fit[[1]]
  if (!is.null(postSetup)) {
    .ret <- postSetup(.ret, .ui, .fit)
  }
  if (adjustOutput) {
    .ret <- .nlmFamilyAdjustOutput(.ret, method)
  }
  .ret$message <- NULL
  if (emitFitWarnings) {
    lapply(.fit[[2]], function(.w) warning(.w, call. = FALSE))
  }
  if (rxode2::rxGetControl(.ui, returnFlag, FALSE)) {
    return(.ret[[method]])
  }
  .ret$message <- message(.ret[[method]])
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- if (is.character(getTheta)) {
    .nlmFamilyGetTheta(.ret[[method]], .ui, getTheta)
  } else {
    getTheta(.ret[[method]], .ui)
  }
  .ret$control <- .control
  .ret$extra <- if (is.function(extra)) extra(.control) else extra
  .nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject(.ret$control, .ret)
  .ret$est <- method
  if (is.character(objective)) {
    .ret$objective <- 2 * as.numeric(.ret[[method]][[objective]])
  } else if (!is.null(objective)) {
    .ret$objective <- objective(.ret[[method]])
  }
  # building the EBE model is another symengine model build; time it as "setup"
  .ret$model <- nlmixrWithTiming("setup", {
    .ui$ebe
  })
  # The control must stay on the ui until the EBE model is built; the build reads
  # optExpression/sumProd/eventSens off of it (#864)
  if (exists("control", .ui)) {
    rm(list = "control", envir = .ui)
  }
  .ret$ofvType <- method
  controlToFocei(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(
    .ret$ui,
    data = .ret$origData,
    control = .ret$control,
    table = .ret$table,
    env = .ret,
    est = method
  )
  .env <- .ret$env
  .env$method <- method
  .ret
}
