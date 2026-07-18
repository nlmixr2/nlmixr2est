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
.nlmSetupEnv <- function(par, ui, data, modelInfo, control,
                         lower=NULL, upper=NULL) {
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
  if (!any(names(.ctl) == "shiErr")) {
    .ctl$shiErr <-   (.Machine$double.eps)^(1/3)
  }

  if (!any(names(.ctl) == "hessErr")) {
    .ctl$hessErr <-   (.Machine$double.eps)^(1/3)
  }

  .env <- new.env(parent=emptyenv())
  .env$rxControl <- .ctl$rxControl
  .env$thetaNames <- names(par)
  .f <- modelInfo
  if (any(names(.f) == "thetaGrad")) {
    .env$predOnly <- .f$predOnly
    nlmixr2global$nlmEnv$model <- .env$thetaGrad <- .f$thetaGrad
    ## sensMethod="adjoint": the thetaGrad model carries the rx__adjFX_* sweep
    ## lhs, so it must be solved with the matching in-engine discrete-adjoint
    ## (s) method.  nlm.cpp calls rxSolve_ directly (bypassing rxSolve's method
    ## auto-upgrade), so set the s-method on rxControl here for every nlm-family
    ## method that consumes thetaGrad.
    .adj <- .nlmAdjointResolve(ui)
    if (isTRUE(.adj$useAdjoint)) {
      .env$rxControl$method <- .adj$sMethodInt
    }
  } else {
    nlmixr2global$nlmEnv$model <- .env$predOnly <- .f$predOnly
  }
  ## Delay differential equation models need a dense-output solver so delay()
  ## history is recorded and interpolated (also by the forward-sensitivity /
  ## jump-sensitivity states).  nlm.cpp calls rxSolve_ directly, bypassing
  ## rxSolve()'s hasDelay enforcement, so replicate it here: dense dop853 (which
  ## needs no analytic Jacobian) unless a discrete-adjoint method (>=200) is
  ## already selected -- those record their own dense history.
  if (isTRUE(rxode2::rxModelVars(nlmixr2global$nlmEnv$model)$flags[["hasDelay"]] == 1L)) {
    .env$rxControl$dense <- TRUE
    if (is.null(.env$rxControl$method) || .env$rxControl$method < 200L) {
      .env$rxControl$method <- 0L
      .env$rxControl$stiff2 <- 0L
    }
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
  ## scaleC solve -- doing it after would mis-scale dosing params from a jump-less gradient. Deactivated in .nlmFreeEnv; no-op for "fd".
  .env$esActive <- FALSE
  if (!is.null(.env$thetaGrad)) {
    .env$esActive <- isTRUE(tryCatch(
      rxode2::rxEventSensLoadModel(.env$thetaGrad),
      error = function(e) FALSE))
  }
  if (is.null(.ctl$scaleC) && .ctl$scaleType == 2L && .ctl$gradTo > 0) {
    .tmp <- .Call(`_nlmixr2est_nlmGetScaleC`, par, .ctl$gradTo)
    if (length(.tmp) == 0L) {
      .ctl$scaleC <- ui$scaleCtheta
      .Call(`_nlmixr2est_nlmSetScaleC`, .ctl$scaleC)
    } else {
      .ctl$scaleC <- .tmp
    }
  } else if (is.null(.ctl$scaleC)) {
    .ctl$scaleC <- ui$scaleCtheta
    .Call(`_nlmixr2est_nlmSetScaleC`, .ctl$scaleC)
  } else if (!is.null(.ctl$scaleC)) {
    .Call(`_nlmixr2est_nlmSetScaleC`, .ctl$scaleC)
  }
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
#' Split a foreign ("sa"/"imp") covariance request off an nlm-family covMethod.
#'
#' The nlm-family kernels compute only their own Hessian-based ("r") or
#' optimizer-native covariance; "sa"/"imp" are recomputed post-fit at the
#' converged estimates by the decoupled engine (see `.covRecompute`).  The
#' covMethod is left as "sa"/"imp" (a valid match.arg choice so the control
#' round-trips); `.nlmFinalizeList` skips the native cov step for it and the
#' central post-fit hook reads the deferred request off the control.
#' @param covMethod resolved nlm-family covMethod
#' @return list(covMethod, deferred)
#' @noRd
.nlmCovMethodDefer <- function(covMethod) {
  if (is.character(covMethod) && length(covMethod) == 1L &&
        covMethod %in% c("sa", "imp")) {
    return(list(covMethod = covMethod, deferred = covMethod))
  }
  list(covMethod = covMethod, deferred = NA_character_)
}

.nlmFinalizeList <- function(env, lst, par="par", printLine=TRUE,
                             hessianCov=TRUE) {
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
  if (inherits(lst, "nls")) {
    .cov <- summary(lst)$cov.unscaled
    .ret$cov <- .Call(`_nlmixr2est_nlmAdjustCov`, .cov, .parScaled)
  } else if (hessianCov && !(.ctl$covMethod %in% c("", "sa", "imp"))) {
    .malert("calculating covariance")
    if (!any(names(.ret) == "hessian")) {
      .p <- setNames(.parScaled, NULL)
      .hess <- nlmixr2Hess(.p, nlmixr2est::.nlmixrNlmFunC)
      .ret$hessian <- .hess
    }
    dimnames(.ret$hessian) <- list(.name, .name)
    .hess <- .ret$hessian
    if (any(is.na(.hess))) {
      .ret$covMethod <- "failed"
    } else {
      # r matrix
      .r <- 0.5 * .hess
      .ch <- try(cholSE(.r), silent = TRUE)
      .covType <- "r"
      if (inherits(.ch, "try-error")) {
        .r2 <- .r %*% .r
        .r2 <- try(sqrtm(.r2), silent=TRUE)
        .covType <- "|r|"
        if (!inherits(.r2, "try-error")) {
          .ch <- try(cholSE(.r), silent=TRUE)
          if (inherits(.ch, "try-error")) {
            .r2 <- .ch # switch to nearPD
          }
        }
        if (inherits(.r2, "try-error")) {
          .covType <- "r+"
          .r2 <- try(nmNearPD(.r), silent=TRUE)
          if (!inherits(.r2, "try-error")) {
            .ch <- try(cholSE(.r), silent=TRUE)
          }
        } else {
          .ch <- try(cholSE(.r), silent=TRUE)
        }
      }
      if (!inherits(.ch, "try-error")) {
        .rinv <- rxode2::rxInv(.ch)
        .rinv <- .rinv %*% t(.rinv)
        .cov <- 2*.rinv
        dimnames(.cov) <- list(.name, .name)
        dimnames(.rinv) <- list(.name, .name)
        .ret$covMethod <- .covType
        if (.ctl$covMethod != "r") {
          .ret$covMethod <- paste0(.covType, " (", .ctl$covMethod, ")")
        } else {
          .ret$covMethod <- .covType
        }
        .ret$cov.scaled <- .cov
        .ret$cov <- .Call(`_nlmixr2est_nlmAdjustCov`, .ret$cov.scaled, .parScaled)
      } else {
        .ret$covMethod <- "failed"
      }
      dimnames(.r) <- list(.name, .name)
      .ret$r <- .r
    }
    .msuccess("done")
  }
  .ret$censInformation <- .Call(`_nlmixr2est_nlmCensInfo`)
  .Call(`_nlmixr2est_nlmWarnings`)
  .nlmFreeEnv()
  .ret
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
  if (is.null(.censInformation) &&
        !is.null(.nlm$censInformation)) {
    .censInformation <- .nlm$censInformation
    ret[[str]][["censInformation"]] <- NULL
  }
  ret$censInformation <- .censInformation

  .parHistData <- ret$parHistData
  if (is.null(.parHistData) &&
        !is.null(.nlm$parHistData)) {
    .parHistData <- .nlm$parHistData
    ret[[str]][["parHistData"]] <- NULL
  }
  ret$parHistData <- .parHistData

  .cov <- ret$cov
  if (is.null(.cov) &&
        !is.null(.nlm$cov)) {
    .cov <- .nlm$cov
    ret[[str]][["cov"]] <- NULL
  }
  ret$cov <- .cov

  .covMethod <- ret$covMethod
  if (is.null(.covMethod) &&
        !is.null(.nlm$covMethod)) {
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
  if (is.null(covNames)) covNames <- character(0)
  vapply(nms, function(.x) {
    if (.x %in% covNames) .x else toupper(.x)
  }, character(1), USE.NAMES = FALSE)
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
  .et <- rxode2::etTrans(dataSav, ui$mv0, addCmt = TRUE,
                         addlKeepsCov = rxControl$addlKeepsCov,
                         addlDropSs = rxControl$addlDropSs,
                         ssAtDoseTime = rxControl$ssAtDoseTime)
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

#' Shared fit driver for the nlm-family estimation methods
#'
#' @param env dispatch environment (provides `ui`, `control`, `data`, `table`)
#' @param method estimation-method string; also the slot the raw fit is stored
#'   under (e.g. `"nlm"` -> `.ret[["nlm"]]`)
#' @param fitModel `function(ui, dataSav)` running the optimizer
#' @param getTheta `function(fit, ui)` returning the full theta vector
#' @param controlToFocei `function(env)` translating the control to a
#'   focei-style control for output assembly
#' @param returnFlag rxode2 control flag name that short-circuits and returns the
#'   raw optimizer result (e.g. `"returnNlm"`)
#' @param message `function(fit)` returning the `$message` (default `fit$message`)
#' @param emitFitWarnings when TRUE, re-emit the warnings collected from
#'   `fitModel` via `warning()` (nlm does this; the others do not)
#' @param extra `$extra` print string, or a `function(control)` returning it
#' @param adjustOutput when TRUE, run `.nlmFamilyAdjustOutput()`
#' @param objective optional `function(fit)` returning the raw objective; when
#'   `NULL` the driver does not set `$objective` (a `postSetup` closure did)
#' @param postSetup optional `function(ret, ui, fitList)` returning a modified
#'   `ret`, run right after the raw fit is stored and before
#'   `.nlmFamilyAdjustOutput()` (for methods that set cov/covMethod/objective
#'   with custom values)
#' @return the assembled nlmixr2 fit (or the raw optimizer result if `returnFlag`)
#' @author Matthew L. Fidler
#' @export
.nlmFamilyFitGeneric <- function(env, method, fitModel, getTheta,
                                 controlToFocei, returnFlag,
                                 objective = NULL,
                                 message = function(fit) fit$message,
                                 emitFitWarnings = FALSE,
                                 extra = "",
                                 adjustOutput = TRUE,
                                 postSetup = NULL) {
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
  .ret$fullTheta <- getTheta(.ret[[method]], .ui)
  .ret$control <- .control
  .ret$extra <- if (is.function(extra)) extra(.control) else extra
  .nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui)) {
    rm(list = "control", envir = .ui)
  }
  .ret$est <- method
  if (!is.null(objective)) {
    .ret$objective <- objective(.ret[[method]])
  }
  # building the EBE model is another symengine model build; time it as "setup"
  .ret$model <- nlmixrWithTiming("setup", {
    .ui$ebe
  })
  .ret$ofvType <- method
  controlToFocei(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData,
                                    control = .ret$control, table = .ret$table,
                                    env = .ret, est = method)
  .env <- .ret$env
  .env$method <- method
  .ret
}
