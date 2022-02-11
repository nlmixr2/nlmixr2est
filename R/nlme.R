#' Control Values for nlme Fit with extra options for nlmixr
#'
#' The values supplied in the function call replace the defaults and
#' a list with all possible arguments is returned.  The returned list
#' is used as the ‘control’ argument to the ‘nlme’ function.
#'
#' @inheritParams nlme::nlmeControl
#' @inheritParams nlme::nlme
#' @param sens Calculate gradients using forward sensitivity
#' @inheritParams foceiControl
#' @export
#' @examples
#' nlmixr2::nlmeControl()
nlmeControl <- function(maxIter = 50, pnlsMaxIter = 7, msMaxIter = 50, minScale = 0.001,
    tolerance = 1e-05, niterEM = 25, pnlsTol = 0.001, msTol = 1e-06,
    returnObject = FALSE, msVerbose = FALSE, msWarnNoConv = TRUE,
    gradHess = TRUE, apVar = TRUE, .relStep = .Machine$double.eps^(1/3),
    minAbsParApVar = 0.05, opt = c("nlminb", "nlm"), natural = TRUE,
    sigma = NULL, optExpression=TRUE, sumProd=FALSE,
    rxControl=rxode2::rxControl(atol=1e-4, rtol=1e-4),
    random=NULL, fixed=NULL, sens=TRUE, ...) {

  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnObject, len=1, any.missing=FALSE)
  checkmate::assertLogical(msVerbose, len=1, any.missing=FALSE)
  checkmate::assertLogical(msWarnNoConv, len=1, any.missing=FALSE)
  checkmate::assertLogical(gradHess, len=1, any.missing=FALSE)
  checkmate::assertLogical(apVar, len=1, any.missing=FALSE)
  checkmate::assertLogical(natural, len=1, any.missing=FALSE)
  checkmate::assertLogical(sens, len=1, any.missing=FALSE)

  checkmate::assertIntegerish(pnlsMaxIter, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(msMaxIter, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(niterEM, len=1, any.missing=FALSE, lower=1)

  checkmate::assertNumeric(minScale, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(pnlsTol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(msTol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(tolerance, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(.relStep, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(minAbsParApVar, len=1, any.missing=FALSE, lower=0)

  if (!inherits(rxControl, "rxControl")) rxControl <- do.call(rxode2::rxControl, rxControl)

  if (is.null(sigma))
    sigma <- 0
  else if (!is.finite(sigma) || length(sigma) != 1 || sigma < 0)
    stop("Within-group std. dev. must be a positive numeric value")
  .ret <- list(maxIter = maxIter, pnlsMaxIter = pnlsMaxIter, msMaxIter = msMaxIter,
               minScale = minScale, tolerance = tolerance, niterEM = niterEM,
               pnlsTol = pnlsTol, msTol = msTol, returnObject = returnObject,
               msVerbose = msVerbose, msWarnNoConv = msWarnNoConv, gradHess = gradHess,
               apVar = apVar, .relStep = .relStep, minAbsParApVar = minAbsParApVar,
               opt = match.arg(opt), natural = natural, sigma = sigma,
               optExpression=optExpression, sumProd=sumProd,
               rxControl=rxControl, sens=sens,
               ...)
  class(.ret) <- "nlmeControl"
  .ret
}
#' Get the nlme family control
#'
#' @param env nlme optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmeFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2::nlmeControl()
  }
  if (!inherits(.control, "nlmeControl")){
    .control <- do.call(nlmixr2::nlmeControl, .control)
  }
  assign("control", .control, envir=.ui)
}

.nlmeFitDataObservations <- NULL
.nlmeFitDataAll   <- NULL
.nlmeFitRxModel   <- NULL
.nlmeFitRxControl <- NULL
.nlmeFitFunction <- NULL
.nlmeGradDimnames <- NULL
.nlmeFitSens <- FALSE


#' A surrogate function for nlme to call for ode solving
#'
#' @param pars Parameters that will be estimated
#' @param id The patient identifiers for the estimated data.
#' @return Predictions
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.nlmixrNlmeFun <- function(pars, id) {
  .ids <- as.character(unique(id))
  .datF <- do.call(rbind, lapply(seq_along(.ids), function(i){
    .datF <- .nlmeFitDataAll[.nlmeFitDataAll$ID == .ids[i], ]
    .datF$ID <- i
    .datF
  }))
  .pars <- as.data.frame(c(pars, list(ID=id)))
  .pars <- .pars[!duplicated(.pars$ID),]
  .pars$ID <- seq_along(.pars$ID)
  row.names(.pars) <- NULL
  .ctl <- .nlmeFitRxControl
  .ctl$returnType <- "data.frame"
  .retF <- do.call(rxode2::rxSolve, c(list(object=.nlmeFitRxModel, params=.pars, events=.datF),
                                      .nlmeFitRxControl))
  .ret <- .retF$rx_pred_
  if (.nlmeFitSens) {
    .grad <- .retF[, paste0("rxD_", names(pars))]
    .grad <- as.matrix(.grad)
    dimnames(.grad) <- .nlmeGradDimnames
    attr(.ret, "gradient") <- .grad
  }
  .ret
}

#' A surrogate function for nlme to call for ode solving
#'
#' @return User function for the saved model
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.nlmixrNlmeUserFun <- function() {
  .nlmeFitFunction
}

#' Setup the data for nlme estimation
#'
#' @param dataSav Formatted Data
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmeFitDataSetup <- function(dataSav) {
  .dsAll <- dataSav[dataSav$EVID != 2, ] # Drop EVID=2 for estimation
  assignInMyNamespace(".nlmeFitDataObservations", nlme::groupedData(DV ~ TIME | ID, .dsAll[.dsAll$EVID == 0, ]))
  assignInMyNamespace(".nlmeFitDataAll", .dsAll)
}

.nlmeFitModel <- function(ui, dataSav, timeVaryingCovariates=.tv) {
  .nlmeFitDataSetup(dataSav)
  assignInMyNamespace(".nlmeFitRxModel", rxode2::rxode2(ui$nlmeRxModel))
  assignInMyNamespace(".nlmeFitFunction", ui$nlmeFunction)
  assignInMyNamespace(".nlmeFitSens", rxode2::rxGetControl(ui, "sens", TRUE))
  assignInMyNamespace(".nlmeGradDimnames", ui$nlmeGradDimnames)
  assignInMyNamespace(".nlmeFitRxControl",  rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl()))

  .ctl <- ui$control
  class(.ctl) <- NULL
  eval(bquote(nlme::nlme(model=.(ui$nlmeModel), data=.nlmeFitDataObservations,
             fixed=.(ui$nlmeFixedFormula), random=.(ui$nlmePdOmega),
             start=.(ui$nlmeStart), weights=.(ui$nlmeWeights),
             control=.(.ctl), na.action=function(object, ...) {
               return(object)
             })))
}


.nlmeFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  .ret$table <- env$table
  .foceiPreProcessData(.data, .ret, .ui)
  .et <- rxode2::etTrans(.ret$dataSav, .ui$mv0, addCmt=TRUE)
  # Just like saem, nlme can use mu-referenced covariates
  .nTv <- attr(class(.et), ".rxode2.lst")$nTv
  if (is.null(.nTv)) .nTv <- 0
  .tv <- character(0)
  if (.nTv != 0) {
    .tv <- names(.et)[-seq(1, 6)]
  }
  .ret$nlme <- .nlmeFitModel(.ui, .ret$dataSav, timeVaryingCovariates=.tv)
  assign("nlmeFit", .ret$nlme, envir=globalenv())
  return(.ret)
  .ret$control <- .control
  nmObjHandleControlObject(.ret$control, .ret)
  .ret$ui <- .ui
  .saemCalcCov(.ret)
  .getSaemTheta(.ret)
  .getSaemOmega(.ret)
  .nlmixr2FitUpdateParams(.ret)
  .saemAddParHist(.ret)
  .saemCalcLikelihood(.ret)
   if (exists("control", .ui)) {
    rm(list="control", envir=.ui)
   }
  .ret$theta <- .ui$saemThetaDataFrame
  .ret$model <- .ui$saemModelPred
  .ret$message <- "" # no message for now
  .ret$est <- "saem"
  .saemControlToFoceiControl(.ret)
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="nlme")
  .env <- .ret$env
  .env$method <- "nlme "
  .ret
}

nlmixr2Est.nlme <- function(env, ...) {
  .ui <- env$ui
  .nlmeFamilyControl(env, ...)
  on.exit({rm("control", envir=.ui)})
  .nlmeFamilyFit(env,  ...)
}


