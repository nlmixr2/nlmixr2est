
#'  Determine if the parameter is a mu-referenced covariate
#'
#' @param expr Expression to check
#' @param muRefCovariateDataFrame Mu Ref data frame
#' @return A boolaen that tells if the expression is a mu-ref covariate
#' @author Matthew L. Fidler
#' @noRd
.saemDropParametersIsMuRefCovariate <- function(expr, muRefCovariateDataFrame) {
  if (length(expr) == 3) {
    if (identical(expr[[1]], quote(`*`))) {
      if (length(expr[[2]]) == 1 &&
            length(expr[[3]]) == 1) {
        .cov1 <- as.character(expr[[2]])
        .cov2 <- as.character(expr[[3]])
        .w <- which(muRefCovariateDataFrame$covariate == .cov1 &
                      muRefCovariateDataFrame$covariateParameter == .cov2)
        if (length(.w) == 1) return(TRUE)
        .w <- which(muRefCovariateDataFrame$covariate == .cov2 &
                      muRefCovariateDataFrame$covariateParameter == .cov1)
        if (length(.w) == 1) return(TRUE)
      }
    }
  }
  FALSE
}
#' Drop mu-referenced parameters
#'
#' @param line Line to change
#' @param muRefDataFrame Mu-referenced data frame
#' @param muRefCovariateDataFrame Mu Referenced Covariates
#' @return Remove mu-referenced etas and covariates
#' @author Matthew L. Fidler
#' @noRd
.saemDropParameters <- function(line, muRefDataFrame, muRefCovariateDataFrame) {
  f <- function(x) {
    if (is.name(x) || is.atomic(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`+`))) {
        if (.saemDropParametersIsMuRefCovariate(x[[2]], muRefCovariateDataFrame)) {
          return(f(x[[3]]))
        }
        if (.saemDropParametersIsMuRefCovariate(x[[3]], muRefCovariateDataFrame)) {
          return(f(x[[2]]))
        }
        if (length(x[[2]]) == 1) {
          .char <- as.character(x[[2]])
          if (.char %in% muRefDataFrame$eta) {
            return(f(x[[3]]))
          }
        }
        if (length(x[[3]]) == 1) {
          .char <- as.character(x[[3]])
          if (.char %in% muRefDataFrame$eta) {
            return(f(x[[2]]))
          }
        }
      }
      as.call(lapply(x, f))
    } else {
      return(x)
    }
  }
  f(line)
}
#' Drop mu referenced etas and covariates
#'
#' @param ui rxode2 ui
#' @return model line expression with mu referenced information dropped.
#' @author Matthew L. Fidler
#' @noRd
.saemDropMuRefFromModel <- function(ui) {
  if (exists("muRefFinal", ui)) {
    .muRefFinal <- ui$muRefFinal
  } else {
    .muRefFinal <- ui$muRefCovariateDataFrame
  }
  .muRefDataFrame <- ui$muRefDataFrame
  lapply(ui$lstExpr, function(line){
    .saemDropParameters(line, .muRefDataFrame, .muRefFinal)
  })
}

#' This is a S3 method for getting the distribution lines for a base rxode2 saem problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the estimation of saem
#' @author Matthew Fidler
#' @keywords internal
#' @export
nmGetDistributionSaemLines <- function(line) {
  UseMethod("nmGetDistributionSaemLines")
}
#' Creates a saem line object from a predDf line
#'
#' @param x rxode2 ui object
#' @param line Line number for saem error line object
#' @param len Number of prediction statements
#' @return nmGetDistributionSaemLines object
#' @author Matthew L. Fidler
#' @noRd
.createSaemLineObject <- function(x, line) {
  .predDf <- get("predDf", x)
  if (line > nrow(.predDf)) {
    return(NULL)
  }
  .predLine <- .predDf[line, ]
  .ret <- list(x, .predLine)
  class(.ret) <- c(paste(.predLine$distribution), "nmGetDistributionSaemLines")
  .ret
}

#' @rdname nmGetDistributionSaemLines
#' @export
nmGetDistributionSaemLines.rxUi <- function(line) {
  .predDf <- get("predDf", line)
  lapply(seq_along(.predDf$cond), function(c){
    .mod <- .createSaemLineObject(line, c)
    nmGetDistributionSaemLines(.mod)
  })
}

#' @rdname nmGetDistributionSaemLines
#' @export
nmGetDistributionSaemLines.norm <- function(line) {
  .rx <- line[[1]]
  .pred1 <- line[[2]]
  if (.pred1[["linCmt"]]) {
    .var <- quote(linCmt())
  } else {
    .var <- .enQuote(.pred1[["var"]])
  }
  return(list(bquote(rx_pred_ <- .(.var))))
}

#' @export
nmGetDistributionSaemLines.t <- function(line) {
  stop("t isn't supported yet")
}

#' @export
nmGetDistributionSaemLines.default  <- function(line) {
  stop("Distribution not supported")
}

#' @export
rxUiGet.saemModel0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=nmGetDistributionSaemLines(.f),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE,
                              lstExpr=.saemDropMuRefFromModel(.f))
}
#attr(rxUiGet.saemModel0, "desc") <- "saem initial model"

#'@export
rxUiGet.saemModelPred0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=rxGetDistributionFoceiLines(.f),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE,
                              lstExpr=.saemDropMuRefFromModel(.f))
}
# attr(rxUiGet.saemModel0, "desc") <- "saem predOnly for use in calculating residuals with focei engine"



#' Load the saem model into symengine
#'
#' @param x rxode2 UI object
#' @return String for loading into symengine
#' @author Matthew L. Fidler
#' @noRd
.saemPrune <- function(x) {
  .x <- x[[1]]
  .x <- .x$saemModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .malert("pruning branches ({.code if}/{.code else}) of saem model...")
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

#' Load the saem predOnly model into symengine
#'
#' @param x rxode2 UI object
#' @return String for loading into symengine
#' @author Matthew L. Fidler
#' @noRd
.saemPrunePred <- function(x) {
  .x <- x[[1]]
  .x <- .x$saemModelPred0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .malert("pruning branches ({.code if}/{.code else}) of saem model...")
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
rxUiGet.loadPruneSaem <- function(x, ...) {
  .loadSymengine(.saemPrune(x), promoteLinSens = FALSE)
}
#attr(rxUiGet.loadPruneSaem, "desc") <- "load the saem model into symengine"

#' @export
rxUiGet.loadPruneSaemPred <- function(x, ...) {
  .loadSymengine(.saemPrunePred(x), promoteLinSens = FALSE)
}
#attr(rxUiGet.loadPruneSaem, "desc") <- "load the saem model into symengine"


#' @export
rxUiGet.saemParamsToEstimate <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  c(.iniDf$name[!is.na(.iniDf$ntheta) & is.na(.iniDf$err)], .ui$nonMuEtas)
}
#attr(rxUiGet.saemParamsToEstimate, "desc") <- "Get the parameters to estimate"

#' @export
rxUiGet.saemThetaName <- rxUiGet.saemParamsToEstimate
#attr(rxUiGet.saemParamsToEstimate, "desc") <- "Get the parameters to estimate"

#' @export
rxUiGet.saemParams <- function(x, ...) {
  .ui <- x[[1]]
  .par <- c(rxUiGet.saemParamsToEstimate(x, ...), .ui$covariates)
  paste0("params(", paste(.par, collapse=","), ")")
}
attr(rxUiGet.saemParams, "desc") <- "Get the params() for a saem model"

#' @export
rxUiGet.saemModel <- function(x, ...) {
  .s <- rxUiGet.loadPruneSaem(x, ...)
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
    .malert("stabilizing round off errors in saem model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "saem model")
     .msuccess("done")
  }
  paste(c(rxUiGet.saemParams(x, ...), rxUiGet.foceiCmtPreModel(x, ...),
          .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")
}

#'@export
rxUiGet.saemModelPredReplaceLst <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .thetaNames <- .iniDf[!is.na(.iniDf$ntheta) & is.na(.iniDf$err), ]
  .thetaValue <- setNames(paste0("THETA[", .thetaNames$ntheta, "]"), .thetaNames$name)
  if (length(.ui$nonMuEtas) > 0) {
    .nonMuThetas <- setNames(rep("", length(.ui$nonMuEtas)), .ui$nonMuEtas)
    .thetaValue <- c(.thetaValue, .nonMuThetas)
  }
  .thetaErrNames <- .iniDf[!is.na(.iniDf$ntheta) & !is.na(.iniDf$err), ]

  .thetaValueErr <- setNames(paste0("THETA[", .thetaErrNames$ntheta, "]"), .thetaErrNames$name)
  .thetaValue <- c(.thetaValue, .thetaValueErr)

  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  for (.e in seq_along(.etaTrans)) {
    .eta <- paste0("ETA[", .e, "]")
    .tn <- .etaTrans[.e]
    if (.thetaValue[.tn] == "") {
      .thetaValue[.tn] <- .eta
    } else {
      .thetaValue[.tn] <- paste0(.thetaValue[.tn], " + ", .eta)
    }
  }
  if (exists("muRefFinal", .ui)) {
    .muRefFinal <- .ui$muRefFinal
  } else {
    .muRefFinal <- .ui$muRefCovariateDataFrame
  }
  for (.c in seq_along(.muRefFinal$theta)) {
    .tv <- .muRefFinal$theta[.c]
    .w <- which(.thetaNames$name == .muRefFinal$covariateParameter[.c])
    if (length(.w) == 1L) {
      .tcov <- paste0("THETA[", .thetaNames$ntheta[.w], "]")
      .tcov <- paste0(.muRefFinal$covariate[.c], " * ", .tcov)
      .cur <- c(.thetaValue[.tv], .tcov)
      .cur <- .cur[.cur != ""]
      .thetaValue[.tv] <- paste(.cur, collapse=" + ")
    }
  }
  .thetaValue
}
#attr(rxUiGet.saemModelPredReplaceLst, "desc") <- "Replace the mu referenced thetas with these values"

.saemModelPredSymengineEnvironment <- NULL

#' @export
rxUiGet.saemModelPred <- function(x, ...) {
  .s <- rxUiGet.loadPruneSaemPred(x, ...)
  .replaceLst <- rxUiGet.saemModelPredReplaceLst(x, ...)
  assignInMyNamespace(".saemModelPredSymengineEnvironment", .s)
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .r <- get("rx_r_", envir = .s)
  .r <- paste("rx_r_=", rxode2::rxFromSE(.r))
  .yj <- paste(get("rx_yj_", envir = .s))
  .yj <- paste0("rx_yj_~", rxode2::rxFromSE(.yj))
  .lambda <- paste(get("rx_lambda_", envir = .s))
  .lambda <- paste0("rx_lambda_~", rxode2::rxFromSE(.lambda))
  .hi <- paste(get("rx_hi_", envir = .s))
  .hi <- paste0("rx_hi_~", rxode2::rxFromSE(.hi))
  .low <- paste(get("rx_low_", envir = .s))
  .low <- paste0("rx_low_~", rxode2::rxFromSE(.low))
  ## if (is.null(.lhs0)) .lhs0 <- ""
  .ui <- x[[1]]
  .lhsIn <- .ui$mv0$lhs
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- ""

  .ret <- paste(c(
    .ddt,
    #.yj,
    #.lambda,
    #.hi,
    #.low,
    .prd#,
    #.r,
    #.s$..lhs,
    #"tad=tad()",
    #"dosenum=dosenum()"
  ), collapse = "\n")
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .ret2 <- c(.r,
             .s$..lhs,
             "tad=tad()",
             "dosenum=dosenum()")

  if (.sumProd) {
    .malert("stabilizing round off errors in saem predOnly model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .ret2 <- rxode2::rxSumProdModel(.ret2)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "saem predOnly model")
    .ret2 <- rxode2::rxOptExpr(.ret2, "saem predOnly model")
    .msuccess("done")
  }
  .ret <- paste(c(
    .yj,
    .lambda,
    .hi,
    .low,
    .ret,
    .ret2
  ), collapse = "\n")
  .ret <- c(rxUiGet.foceiParams(x, ...),
            rxUiGet.foceiCmtPreModel(x, ...),
            "rx_pred_=NA\nrx_r_=NA\n",
            paste(names(.replaceLst), "<-", .replaceLst),
            .ret,
            vapply(.uiGetThetaEta(x[[1]]), deparse1, character(1), USE.NAMES=FALSE),
            .foceiToCmtLinesAndDvid(x[[1]]))
  .ret <- .ret[.ret != ""]
  .ret <- list(predOnly=rxode2::rxode2(paste(.ret, collapse="\n")))
  class(.ret) <- "saemModelList"
  .ret
}