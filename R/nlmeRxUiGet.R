#' This is a S3 method for getting the distribution lines for a base rxode2 nlme problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the focei. This is based
#'   on the idea that the focei parameters are defined
#' @author Matthew Fidler
#' @keywords internal
#' @export
rxGetDistributionNlmeLines <- function(line) {
  UseMethod("rxGetDistributionNlmeLines")
}
#' Get the DV transformation
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @param yj The transformation number for the current error
#' @return The transformation expression
#' @author Matthew Fidler
#' @keywords internal
#' @export
.rxGetDVFTransform <- function(env, pred1, yj) {
  if (yj == 2) {
    return(quote(DVIN))
  } else if (yj == 3) {
    return(quote(log(DVIN)))
  } else {
    return(quote(rxTBS(DVIN, rx_lambda_, rx_yj_, rx_low_, rx_hi_)))
  }
}
#' @export
rxGetDistributionNlmeLines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  .ret <- vector("list", 7)
  .yj <- as.double(pred1$transform) - 1
  .ret[[1]] <- bquote(rx_yj_ ~ .(.yj))
  .ret[[2]] <- bquote(rx_lambda_~.(rxode2::.rxGetLambdaFromPred1AndIni(env, pred1)))
  .ret[[3]] <- bquote(rx_low_ ~ .(rxode2::.rxGetLowBoundaryPred1AndIni(env, pred1)))
  .ret[[4]] <- bquote(rx_hi_ ~ .(rxode2::.rxGetHiBoundaryPred1AndIni(env, pred1)))
  .ret[[5]] <- bquote(rx_pred_ ~ .(rxode2::.rxGetPredictionFTransform(env, pred1, .yj)))
  .ret[[6]] <- bquote(rx_dv_ ~ .(.rxGetDVFTransform(env, pred1, .yj)))
  .ret[[7]] <- quote(rx_diff_ <- rx_pred_ - rx_dv_)
  .ret
}

#' @export
rxGetDistributionNlmeLines.t <- function(line) {
  stop("t isn't supported yet")
}

#' @export
rxGetDistributionNlmeLines.default  <- function(line) {
  print(line[[1]])
  print(line[[2]])
  stop("Distribution not supported")
}

#' @export
rxGetDistributionNlmeLines.rxUi <- function(line) {
  .predDf <- get("predDf", line)
  if (length(.predDf$cond) != 1) stop("nlme does not support multiple endpoint models")
  lapply(seq_along(.predDf$cond), function(c){
    .mod <- .createFoceiLineObject(line, c)
    rxGetDistributionNlmeLines(.mod)
  })
}


#' Load the nlme model into symengine
#'
#' @param x rxode2 UI object
#' @return String for loading into symengine
#' @author Matthew L. Fidler
#' @noRd
.nlmePrune <- function(x) {
  .x <- x[[1]]
  .x <- .x$nlmeModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .malert("pruning branches ({.code if}/{.code else}) of nlme model...")
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
rxUiGet.nlmeModel0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=rxGetDistributionNlmeLines(.f),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE,
                              lstExpr=.saemDropMuRefFromModel(.f))
}
#attr(rxUiGet.nlmeModel, "desc") <- "nlmixr nlme model, equivalent to saem rxode2 model"



#' @export
rxUiGet.loadPruneNlme <- function(x, ...) {
  .loadSymengine(.nlmePrune(x), promoteLinSens = FALSE)
}

#' @export
rxUiGet.nlmeFunction <- function(x, ...) {
  .ui <- x[[1]]
  .estPar <- rxUiGet.saemParamsToEstimate(x, ...)
  #.par <- c(.estPar, .ui$covariates)
  .par <- .estPar
  eval(parse(text=paste0("function(", paste(.par, collapse=","), ", TIME, ID) {\n",
                         "nlmixr2::.nlmixrNlmeFun(c(", paste(paste0(.estPar, "=", .estPar), collapse=","), "))\n",
                         "}")))
}

#' @export
rxUiGet.nlmeModel <- function(x, ...) {
  .s <- rxUiGet.loadPruneNlme(x, ...)
  .prd <- get("rx_diff_", envir = .s)
  .prd <- paste0("rx_diff_=", rxode2::rxFromSE(.prd))
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
    .malert("stabilizing round off errors in nlme model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "nlme model")
     .msuccess("done")
  }
  paste(c(rxUiGet.saemParams(x, ...), rxUiGet.foceiCmtPreModel(x, ...),
          .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")
}

