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
