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
#' @param modelInfo A list containing the following elements:
#'
#' - `predOnly` -- A model with only predictions calculated.  These
#' predictions should be in terms of `THETA[#]` and `DV`.  The
#'
#' - `eventTheta` is an indicator if the `THETA[#]` is related to an
#' event (like `dur(x)` `f(x)`).  These variables will use Shi2021
#' finite differences and need to be indicated when setting up the
#' solving environment.  When finite differences are required, this is
#' `1L` when they are not it should be `0L`.  This should match the
#' length of `par`
#'
#' - `thetaGrad` -- needed when solveType != 1; a model that gives the
#' value and gradient of each `THETA[#]`
#'
#' An example can be found with `ui$nlmSensModel` or `ui$nlmRxModel`
#'
#' @param control is a control structure with a few required elements:
#'
#' - `rxControl` represents the rxode2 solving options
#' - `solveType` integer indicating the solveType (optional)
#' - `stickyRecalcN`
#' - `maxOdeRecalc`
#' - `odeRecalcFactor`
#' - `eventType` (optional)
#' - `shi21maxFD` (optional)
#' - `shiErr` (optional)
#' - `optimHessType` (optional)
#' - `shi21maxHess` (optional)
#' - `hessErr` (optional)
#' - `useColor`
#' - `printNcol`
#' - `print`
#' - `normType`
#' - `scaleType`
#' - `scaleCmin`
#' - `scaleCmax`
#' - `scaleTo`
#' - `scaleC`
#' - `gradTo` (optional); if missing assumed gradTo=0
#' @param lower lower bounds, will be scaled if present
#' @param upper upper bounds, will be scaled if present
#' @return nlm solve environment; of interest
#'
#' `$par.ini` -- scaled parameter initial value
#'
#' `$lower` -- scaled parameter lower value
#'
#' `$upper` -- scaled parameter upper value
#'
#' `$.ctl`  -- control structure
#'
#' @details
#'
#' In between using this, rxode2 solving should not be called.
#'
#' This will also print the header for solving (if print != 0)
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
  } else {
    nlmixr2global$nlmEnv$model <- .env$predOnly <- .f$predOnly
  }
  .env$param <- setNames(par, sprintf("THETA[%d]", seq_along(par)))
  .nlmFitDataSetup(data)
  .env$needFD <- .f$eventTheta
  # Iteration-print transforms: ship a single xform sub-list so the C
  # side wires log/logit/probit back-transforms through one helper
  # (scaleAttachXform in src/scale.h), identical to every other method
  # (focei: env$xform, saem: .cfg$xform, nlm-family: .ctl$xform).
  # nlm-family fits are population-only, so the printed parameters are
  # just thetas in `par` order.
  .ctl$xform <- .iterPrintXParFromUi(ui, names(par))
  .env$control <- .ctl
  .env$data <- nlmixr2global$nlmEnv$data
  .Call(`_nlmixr2est_nlmSetup`, .env)
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
  } else if (hessianCov && .ctl$covMethod != "") {
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
#' nlmixr2 normalizes input/event-table column names to upper case, but leaves
#' model covariate columns at their declared case.  This is the single shared
#' implementation used by `.foceiPreProcessData()`, the covariate-present
#' pre-process hook, and the output-table re-insertion helpers.
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

#' Detect the time-varying covariate columns for mu-referenced estimators
#'
#' SAEM and NLME both support mu-referenced covariates and must know which
#' covariate columns vary within a subject.  This is the single shared
#' implementation; it translates the (already preprocessed) data with the base
#' model vars and returns the time-varying covariate column names.
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
  # nTv == 0 -> no time-varying covariates; otherwise (or, defensively, when the
  # attribute is absent) the time-varying covariate columns follow the first 6
  # event-table columns.  This preserves the prior saem behavior exactly.
  if (!is.null(.nTv) && .nTv == 0L) {
    return(character(0))
  }
  names(.et)[-seq_len(6)]
}

#' Generic family-control setup for the nlm-family estimation methods
#'
#' Every nlm-family method (`nlm`, `nlminb`, `bobyqa`, `newuoa`, `uobyqa`,
#' `n1qn1`, `lbfgsb3c`, `optim`, `nls`) sets up its control identically: take the
#' control from the dispatch env, default it, coerce a plain list to the proper
#' control object, and assign it onto the ui.  This is that shared body; the
#' per-method `.<m>FamilyControl` wrappers pass their `*Control()` function and
#' its class name.
#'
#' @param env dispatch environment (provides `ui` and `control`)
#' @param controlFn the method's `*Control()` constructor (e.g. `nlmControl`)
#' @param controlClass the control object's S3 class (e.g. `"nlmControl"`)
#' @return Nothing; assigns the resolved control onto `env$ui`
#' @author Matthew L. Fidler
#' @noRd
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

#' Generic family-fit driver for the nlm-family estimation methods
#'
#' Every nlm-family method shares the same fit spine: preprocess the data, run
#' the method's optimizer (`fitModel`) collecting warnings, assemble the fit
#' environment (objective, theta, control, message, ...), and hand it to
#' `nlmixr2CreateOutputFromUi()`.  The per-method `.<m>FamilyFit` wrappers supply
#' only what genuinely differs.  The fit environment built here needs (and this
#' driver populates): table, origData/dataSav/idLvl/covLvl (via
#' `.foceiPreProcessData`), ui, adjObf, fullTheta, control, extra, est,
#' objective, model, ofvType, message, theta.
#'
#' @param env dispatch environment (provides `ui`, `control`, `data`, `table`)
#' @param method estimation-method string; also the slot the raw fit is stored
#'   under (e.g. `"nlm"` -> `.ret[["nlm"]]`)
#' @param fitModel `function(ui, dataSav)` running the optimizer
#' @param getTheta `function(fit, ui)` returning the full theta vector
#' @param objective `function(fit)` returning the raw objective (driver does not
#'   multiply; pass e.g. `function(f) 2 * as.numeric(f$minimum)`)
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
#'   `.nlmFamilyAdjustOutput()` -- for methods that must set cov/covMethod/
#'   objective with custom values that adjustOutput's `is.null` guards then keep
#' @return the assembled nlmixr2 fit (or the raw optimizer result if `returnFlag`)
#' @author Matthew L. Fidler
#' @noRd
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
  .foceiPreProcessData(.data, .ret, .ui, .control$rxControl)
  .fit <- .collectWarn(fitModel(.ui, .ret$dataSav), lst = TRUE)
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
  .ret$model <- .ui$ebe
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
