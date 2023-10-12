#' nlmixr2 defaults controls for nlm
#'
#' @inheritParams stats::nlm
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @return
#' @export
#' @author Matthew L. Fidler
#' @examples
nlmControl <- function(hessian = TRUE, typsize = NULL,
                       fscale = 1, print.level = 2, ndigit = NULL, gradtol = 1e-6,
                       stepmax = NULL,
                       steptol = 1e-6, iterlim = 10000, check.analyticals = TRUE,
                       rxControl=NULL,
                       optExpression=TRUE, sumProd=FALSE,
                       returnNlm=FALSE,
                       addProp = c("combined2", "combined1"),
                       calcTables=TRUE, compress=TRUE,
                       adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {
  checkmate::assertLogical(hessian, len=1, any.missing=FALSE)
  checkmate::assertNumeric(stepmax, lower=0, len=1, null.ok=TRUE, any.missing=FALSE)
  checkmate::assertIntegerish(print.level, lower=1, upper=2, any.missing=FALSE)
  checkmate::assertNumeric(ndigit, lower=0, len=1, any.missing=FALSE, null.ok=TRUE)
  checkmate::assertNumeric(gradtol, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(steptol, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(iterlim, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(check.analyticals, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnNlm, len=1, any.missing=FALSE)
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

  .ret <- list(hessian = hessian, typsize = typsize,
               fscale = fscale, print.level = print.level, ndigit=ndigit, gradtol = gradtol,
               stepmax = stepmax,
               steptol = steptol, iterlim = iterlim,
               check.analyticals = check.analyticals,
               optExpression=optExpression,
               sumProd=sumProd,
               rxControl=rxControl,
               returnNlm=returnNlm, addProp=addProp, calcTables=calcTables,
               compress=compress,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "nlmControl"
  .ret
}

#' Get the nlm family control
#'
#' @param env nlme optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2est::nlmControl()
  }
  if (!inherits(.control, "nlmControl")){
    .control <- do.call(nlmixr2est::nlmControl, .control)
  }
  assign("control", .control, envir=.ui)
}


.nlmEnv <- new.env(parent=emptyenv())

#' A surrogate function for nlm to call for ode solving
#'
#' @param dv The observations for the `nlm` function
#' @param pars Parameters that will be estimated
#' @return Predictions
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.nlmixrNlmFun <- function(pars) {
  row.names(.pars) <- NULL
  .retF <- do.call(rxode2::rxSolve, c(list(object=.nlmEnv$model, params=.pars,
                                           events=.nlmEnv$data),
                                      .nlmEnv$rxControl))
  sum(.retF$rx_pred_)
}

#' Setup the data for nlm estimation
#'
#' @param dataSav Formatted Data
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmFitDataSetup <- function(dataSav, model) {
  .dsAll <- dataSav[dataSav$EVID != 2, ] # Drop EVID=2 for estimation
  .nlmEnv$data <- .dsAll
}

#'@export
rxUiGet.nlmModel0 <- function(x, ...) {
  .ui <- rxode2::rxUiDecompress(x[[1]])
  .predDf <- .ui$predDf
  .save <- .predDf
  .predDf[.predDf$distribution == "norm", "distribution"] <- "dnorm"
  assign("predDf", .predDf, envir=.ui)
  on.exit(assign("predDf", .save, envir=.ui))
  .ui$foceiModel0ll
}

#' Load the saem model into symengine
#'
#' @param x rxode2 UI object
#' @return String for loading into symengine
#' @author Matthew L. Fidler
#' @noRd
.nlmPrune <- function(x) {
  .x <- x[[1]]
  .x <- .x$nlmModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .malert("pruning branches ({.code if}/{.code else}) of nlm model...")
  .ret <- rxode2::.rxPrune(.x, envir = .env)
  .mv <- rxode2::rxModelVars(.ret)
  ## Need to convert to a function
  if (rxode2::.rxIsLinCmt() == 1L) {
    .vars <- c(.mv$params, .mv$lhs, .mv$slhs)
    .mv <- rxode2::.rxLinCmtGen(length(.mv$state), .vars)
  }
  .msuccess("done")
  rxode2::rxNorm(.mv)
}

#' @export
rxUiGet.loadPruneNlm <- function(x, ...) {
  .loadSymengine(.nlmPrune(x), promoteLinSens = FALSE)
}

#' @export
rxUiGet.nlmRxModel <- function(x, ...) {
  .s <- rxUiGet.loadPruneNlm(x, ...)
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  ## .lhs0 <- .s$..lhs0
  ## if (is.null(.lhs0)) .lhs0 <- ""
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- ""
  .ret <- paste(c(
    #.s$..stateInfo["state"],
    #.lhs0,
    .ddt,
    .prd,
    #.s$..stateInfo["statef"],
    #.s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  if (.sumProd) {
    .malert("stabilizing round off errors in nlm model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "nlm model")
    .msuccess("done")
  }
  paste(c(rxUiGet.foceiParams(x, ...), rxUiGet.foceiCmtPreModel(x, ...),
          .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")
}


#' @export
rxUiGet.nlmParNameFun <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .env <- new.env(parent=emptyenv())
  .env$i <- 1
  eval(str2lang(
    paste0("function(p) {c(",
           paste(vapply(seq_along(.iniDf$ntheta), function(t) {
             if (.iniDf$fix[t]) {
               paste0("'THETA[", t, "]'=", .iniDf$est[t])
             } else {
               .ret <- paste0("'THETA[", t, "]'=p[", .env$i, "]")
               .env$i <- .env$i + 1
               .ret
             }
           }, character(1), USE.NAMES=FALSE), collapse=","), ")}")))
}
