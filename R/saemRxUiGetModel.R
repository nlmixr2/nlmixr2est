#'  Determine if the parameter is a mu-referenced covariate
#'
#' @param expr Expression to check
#' @param muRefCovariateDataFrame Mu Ref data frame
#' @param noCovs Do not look for covariates
#' @return A boolean that tells if the expression is a mu-ref
#'   covariate
#' @author Matthew L. Fidler
#' @noRd
.saemDropParametersIsMuRefCovariate <- function(expr, muRefCovariateDataFrame, noCovs=FALSE) {
  if (noCovs) return(FALSE)
  if (length(expr) == 3) {
    if (identical(expr[[1]], quote(`*`))) {
      if (length(expr[[2]]) == 1 &&
            length(expr[[3]]) == 1) {
        .cov1 <- as.character(expr[[2]])
        .cov2 <- as.character(expr[[3]])
        .w <- which(muRefCovariateDataFrame$covariate == .cov1)
        if (length(.w) >= 1) {
          return(any(muRefCovariateDataFrame$covariateParameter[.w] == .cov2))
        }
        .w <- which(muRefCovariateDataFrame$covariate == .cov2)
        if (length(.w) >= 1) {
          return(any(muRefCovariateDataFrame$covariateParameter[.w] == .cov1))
        }
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
#' @param noCovs Do not look for covariates
#' @return Remove mu-referenced etas and covariates
#' @author Matthew L. Fidler
#' @noRd
.saemDropParameters <- function(line, muRefDataFrame, muRefCovariateDataFrame, noCovs=FALSE) {
  f <- function(x) {
    if (is.name(x) || is.atomic(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`+`)) &&
            length(x) == 2) {
        return(f(x[[2]]))
      }
      if (identical(x[[1]], quote(`+`))) {
        if (.saemDropParametersIsMuRefCovariate(x[[2]], muRefCovariateDataFrame, noCovs=noCovs)) {
          return(f(x[[3]]))
        }
        if (.saemDropParametersIsMuRefCovariate(x[[3]], muRefCovariateDataFrame, noCovs=noCovs)) {
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
#' @param noCovs Do not look for covariates
#' @return model line expression with mu referenced information dropped.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.saemDropMuRefFromModel <- function(ui, noCovs=FALSE) {
  .muRefFinal <- ui$saemMuRefCovariateDataFrame
  .muRefDataFrame <- ui$muRefDataFrame
  lapply(ui$lstExpr, function(line){
    .saemDropParameters(line, .muRefDataFrame, .muRefFinal, noCovs=noCovs)
  })
}

#' This is a S3 method for getting the distribution lines for a base
#' rxode2 saem problem
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
rxUiGet.saemParamsLine <- function(x, ...) {
  .x <- x[[1]]
  .names <- .x$iniDf[!is.na(.x$iniDf$ntheta) & is.na(.x$iniDf$err), "name"]
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  .names <- .names[!(.names %in% .cov$covariateParameter)]
  str2lang(paste0("param(", paste(.names, collapse=", "), ")"))
}
attr(rxUiGet.saemParamsLine, "rstudio") <- quote(param(tcl))

#' @export
rxUiGet.saemModel0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=nmGetDistributionSaemLines(.f),
                              paramsLine=NA,
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE,
                              lstExpr=.saemDropMuRefFromModel(.f))
}
#attr(rxUiGet.saemModel0, "desc") <- "saem initial model"
attr(rxUiGet.saemModel0, "rstudio") <- quote(rxModelVars({}))

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
attr(rxUiGet.saemModelPred0, "rstudio") <- quote(rxModelVars({}))


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
  .ret <- rxode2::.rxPrune(.x, envir = .env,
                           strAssign=rxode2::rxModelVars(x[[1]])$strAssign)
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
  .ui0 <- .x
  .x <- .x$saemModelPred0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .malert("pruning branches ({.code if}/{.code else}) of saem model...")
  .ret <- rxode2::.rxPrune(.x, envir = .env,
                           strAssign=rxode2::rxModelVars(.ui0)$strAssign)
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
attr(rxUiGet.loadPruneSaem, "rstudio") <- emptyenv()

#' @export
rxUiGet.loadPruneSaemPred <- function(x, ...) {
  .loadSymengine(.saemPrunePred(x), promoteLinSens = FALSE)
}
#attr(rxUiGet.loadPruneSaem, "desc") <- "load the saem model into symengine"
attr(rxUiGet.loadPruneSaemPred, "rstudio") <- emptyenv()


#' @export
rxUiGet.saemParamsToEstimate <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .ret <- c(.iniDf$name[!is.na(.iniDf$ntheta) & is.na(.iniDf$err)])
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  if (length(.cov$theta) > 0) {
    .theta <- .ret
    .theta <- .theta[!(.theta %in% .cov$covariateParameter)]
    .allCovs <- rxUiGet.saemCovars(x, ...)
    .lc <- length(.allCovs)
    .m <- matrix(rep(NA_character_, .lc * length(.theta)), ncol = .lc)
    dimnames(.m) <- list(.theta, .allCovs)
    for (.c in seq_along(.cov$covariateParameter)) {
      .curTheta <- .cov[.c, "theta"]
      .curCov <- .cov[.c, "covariate"]
      .curPar <- .cov[.c, "covariateParameter"]
      .m[.curTheta, .curCov] <- .curPar
    }
    .m <- cbind(matrix(.theta, ncol=1), .m)
    .m <- as.vector(t(.m))
    .ret <- .m[!is.na(.m)]
  }
  c(.ret, .ui$nonMuEtas)
}
#attr(rxUiGet.saemParamsToEstimate, "desc") <- "Get the parameters to estimate"
attr(rxUiGet.saemParamsToEstimate, "rstudio") <- "tcl"

#' @export
rxUiGet.saemParamsToEstimateCov <- function(x, ...) {
  .pars <- rxUiGet.saemParamsToEstimate(x, ...)
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
  .pars[!(.pars %in% .cov$covariateParameter)]
}
attr(rxUiGet.saemParamsToEstimateCov, "rstudio") <- "tcl"

#' @export
rxUiGet.saemThetaName <- rxUiGet.saemParamsToEstimate
#attr(rxUiGet.saemParamsToEstimate, "desc") <- "Get the parameters to estimate"

#' @export
rxUiGet.saemParams <- function(x, ...) {
  .ui <- x[[1]]
  .par <- c(rxUiGet.saemParamsToEstimateCov(x, ...), .ui$covariates)
  paste0("params(", paste(.par, collapse=","), ")")
}
attr(rxUiGet.saemParams, "desc") <- "Get the params() for a saem model"
attr(rxUiGet.saemParams, "rstudio") <- "params(tka)"

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
  .cmt <-  rxUiGet.foceiCmtPreModel(x, ...)
  .interp <- rxUiGet.interpLinesStr(x, ...)
  if (.interp != "") {
    .cmt <-paste0(.cmt, "\n", .interp)
  }
  paste(c(rxUiGet.saemParams(x, ...), .cmt,
          .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")
}
attr(rxUiGet.saemModel, "rstudio") <- "params(tcl)"

#'@export
rxUiGet.saemModelPredReplaceLst <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .thetaNames <- .iniDf[!is.na(.iniDf$ntheta) & is.na(.iniDf$err), ]
  .etas <- .iniDf[which(.iniDf$neta1 == .iniDf$neta2),"name"]
  if (length(.thetaNames$name) == 0L) {
    .thetaValue <- character(0L)
  } else {
    .thetaValue <- setNames(paste0("THETA[", .thetaNames$ntheta, "]"), .thetaNames$name)
  }
  if (length(.ui$nonMuEtas) > 0) {
    .nonMuThetas <- setNames(rep("", length(.ui$nonMuEtas)), .ui$nonMuEtas)
    .thetaValue <- c(.thetaValue, .nonMuThetas)
  }
  .thetaErrNames <- .iniDf[!is.na(.iniDf$ntheta) & !is.na(.iniDf$err), ]

  .thetaValueErr <- setNames(paste0("THETA[", .thetaErrNames$ntheta, "]"), .thetaErrNames$name)
  .thetaValue <- c(.thetaValue, .thetaValueErr)
  .etaTrans <- rxUiGet.saemEtaTransPred(x, ...)
  for (.e in seq_along(.etaTrans)) {
    .eta <- paste0("ETA[", .e, "]")
    .tn <- .etaTrans[.e]
    if (.tn < 0) {
      .thetaValue[.etas[-.tn]] <- .eta
    } else {
      if (.thetaValue[.tn] == "") {
        .thetaValue[.tn] <- .eta
      } else {
        .thetaValue[.tn] <- paste0(.thetaValue[.tn], " + ", .eta)
      }
    }
  }
  .muRefFinal <- rxUiGet.saemMuRefCovariateDataFrame(x, ...)
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
attr(rxUiGet.saemModelPredReplaceLst, "rstudio") <- c(tka="THETA[1] + ETA[1]")

.saemModelEnv <- new.env(parent = emptyenv())
.saemModelEnv$symengine <- NULL
.saemModelEnv$predSymengine <- NULL

#' @export
rxUiGet.interpLinesStr <- function(x, ...) {
  .ui <- x[[1]]
  .interp <- x[[1]]$interpLines
  if (is.null(.interp)) {
    .interp <- ""
  } else {
    .interp <- vapply(.interp, deparse1, character(1), USE.NAMES = FALSE)
  }
  .interp
}
attr(rxUiGet.interpLinesStr, "rstudio") <- ""

#' @export
rxUiGet.saemModelPred <- function(x, ...) {
  .ui0 <- x[[1]]
  .levels  <- .ui0$levels
  if (!is.null(.levels)) {
    .levels <- vapply(seq_along(.levels),
                      function(i){
                        deparse1(.levels[[i]])
                      },
                      character(1), USE.NAMES=FALSE)
    .levels <- paste(.levels, collapse="\n")
  }
  .s <- rxUiGet.loadPruneSaemPred(x, ...)
  .saemModelEnv$symengine <- .s
  .replaceLst <- rxUiGet.saemModelPredReplaceLst(x, ...)

  .saemModelEnv$predSymengine <- .s
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
  .ret0 <- c(.yj,
             .lambda,
             .hi,
             .low)
  .ret2 <- c(.r,
             .s$..lhs,
             "tad=tad()",
             "dosenum=dosenum()")

  if (.sumProd) {
    .malert("stabilizing round off errors in saem predOnly model...")
    .ret0 <- rxode2::rxSumProdModel(.ret0)
    .ret <- rxode2::rxSumProdModel(.ret)
    .ret2 <- rxode2::rxSumProdModel(.ret2)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret0 <- gsub("rx_expr_", "rx_expr", rxode2::rxOptExpr(.ret0, "saem predOnly model 0"))
    .ret <- rxode2::rxOptExpr(.ret, "saem predOnly model 1")
    .ret2 <- gsub("rx_expr_", "rx_expr__", rxode2::rxOptExpr(.ret2, "saem predOnly model 2"))
    .msuccess("done")
  }
  .ret <- paste(c(
    .ret0,
    .ret,
    .ret2
  ), collapse = "\n")
  .interp <- rxUiGet.interpLinesStr(x, ...)
  .ret <- c(rxUiGet.foceiParams(x, ...),
            rxUiGet.foceiCmtPreModel(x, ...),
            .interp,
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
