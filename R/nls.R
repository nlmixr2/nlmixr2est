#' nlmixr2 defaults controls for nls
#'
#' @inheritParams stats::nls
#' @inheritParams stats::nls.control
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @return nls control object
#' @export
#' @author Matthew L. Fidler
#' @examples
nlsControl <- function(maxiter=10000,
                       tol = 1e-05,
                       minFactor = 1/1024,
                       printEval = FALSE,
                       warnOnly = FALSE,
                       scaleOffset = 0,
                       nDcentral = FALSE,
                       algorithm = c("default", "plinear", "port"),
                       trace = TRUE,
                       rxControl=NULL,
                       optExpression=TRUE, sumProd=FALSE,
                       returnNls=FALSE,
                       addProp = c("combined2", "combined1"),
                       calcTables=TRUE, compress=TRUE,
                       adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {
  algorithm <- match.arg(algorithm)
  checkmate::assertLogical(trace, len=1, any.missing=FALSE)
  checkmate::assertLogical(nDcentral, len=1, any.missing=FALSE)
  checkmate::assertNumeric(scaleOffset, any.missing=FALSE, finite=TRUE)
  checkmate::assertLogical(warnOnly, len=1, any.missing = FALSE)
  checkmate::assertLogical(printEval, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(maxiter, len=1, any.missing=FALSE, lower=1)
  checkmate::assertNumeric(tol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(minFactor, len=1, any.missing=FALSE, lower=0)
  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnNls, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(ndigit)) {
    ndigit <- sigdig
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- rxode2::rxControl(sigdig=sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol=1e-4, rtol=1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'", call=FALSE)
  }
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower=1, finite=TRUE, any.missing=TRUE, len=1)
    if (is.null(sigdigTable)) {
      sigdigTable <- round(sigdig)
    }
  }
  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower=1, len=1, any.missing=FALSE)

  .ret <- list(algorithm=algorithm, maxiter=maxiter,
               tol=tol,
               minFactor=minFactor,
               printEval=printEval,
               warnOnly=warnOnly,
               scaleOffset=scaleOffset,
               nDcentral=nDcentral,
               optExpression=optExpression,
               sumProd=sumProd,
               rxControl=rxControl,
               returnNls=returnNls, addProp=addProp, calcTables=calcTables,
               compress=compress,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "nlsControl"
  .ret
}

#' Get the nls family control
#'
#' @param env nlme optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlsFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2est::nlsControl()
  }
  if (!inherits(.control, "nlsControl")){
    .control <- do.call(nlmixr2est::nlsControl, .control)
  }
  assign("control", .control, envir=.ui)
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.nlsControl <- function(control, env) {
  assign("nlsControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.nls <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nlsControl", .env)) {
    .control <- get("nlsControl", .env)
    if (inherits(.control, "nlsControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nlsControl")) return(.control)
  }
  stop("cannot find nls related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.nls <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- nlsControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("nlsControl", .ctl)
  if (!inherits(.ctl, "nlsControl")) {
    .minfo("invalid control for `est=\"nls\"`, using default")
    .ctl <- nlsControl()
  } else {
    .ctl <- do.call(nlsControl, .ctl)
  }
  .ctl
}

.nlsEnv <- new.env(parent=emptyenv())

#' A surrogate function for nls to call for ode solving
#'
#' @param pars Parameters that will be estimated
#' @return Predictions
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.nlmixrNlsFun <- function(...) {
  .retF <- do.call(rxode2::rxSolve,
                   c(list(object=.nlsEnv$model,
                          params=.nlsEnv$parFun(...),
                          events=.nlsEnv$data),
                     .nlsEnv$rxControl))
  # rx_dv_ applies TBS to DV
  # rx_pred_ is the predicted value
  # rx_r_ is the variance
  #
  # Give the weighted residual for nls
  (.retF$rx_dv_ - .retF$rx_pred_) / sqrt(.retF$rx_r_)
}

#'@export
rxUiGet.nlsModel0 <- function(x, ...) {
  .ui <- rxode2::rxUiDecompress(x[[1]])
  .ui$foceiModel0ll
}


#' This is a S3 method for getting the distribution lines for a base rxode2 nls problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the focei. This is based
#'   on the idea that the focei parameters are defined
#' @author Matthew Fidler
#' @keywords internal
#' @export
rxGetDistributionNlsLines <- function(line) {
  UseMethod("rxGetDistributionNlsLines")
}

#' @rdname rxGetDistributionNlsLines
#' @export
rxGetDistributionNlsLines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  .errNum <- line[[3]]
  .line <- rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1, .errNum,
                                                        rxPredLlik=.getRxPredLlikOption())
  .yj <- as.double(pred1$transform) - 1
  if (.yj == 2) {
    .lineExtra <- quote(rx_dv_ ~ DV)
  } else if (.yj == 3) {
    .lineExtra <- quote(rx_dv_ ~ log(DV))
  } else {
    .lineExtra <- quote(rx_dv_ ~ rxTBS(DV, rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  }
  c(.line, list(.lineExtra))
}

#' @rdname rxGetDistributionNlsLines
#' @export
rxGetDistributionNlsLines.default <- function(line) {
  stop("only normally related endoints can be used with 'nls', try 'nlm'", call.=FALSE)
}

#' @export
rxGetDistributionNlsLines.rxUi <- function(line) {
  .predDf <- rxUiGet.predDfFocei(list(line, TRUE))
  lapply(seq_along(.predDf$cond), function(c) {
    .mod <- .createFoceiLineObject(line, c)
    rxGetDistributionNlsLines(.mod)
  })
}

#' @export
rxUiGet.nlsModel0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=rxGetDistributionNlsLines(.f),
                              prefixLines=.uiGetThetaEta(.f),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE)
}


#summary(fm1DNase1)$cov.unscaled
