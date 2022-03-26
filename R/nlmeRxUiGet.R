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
    return(quote(rxDV))
  } else if (yj == 3) {
    return(quote(log(rxDV)))
  } else {
    return(quote(rxTBS(rxDV, rx_lambda_, rx_yj_, rx_low_, rx_hi_)))
  }
}
#' @export
rxGetDistributionNlmeLines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  .ret <- vector("list", 6)
  .yj <- as.double(pred1$transform) - 1
  .ret[[1]] <- bquote(rx_yj_ ~ .(.yj))
  .ret[[2]] <- bquote(rx_lambda_~.(rxode2::.rxGetLambdaFromPred1AndIni(env, pred1)))
  .ret[[3]] <- bquote(rx_low_ ~ .(rxode2::.rxGetLowBoundaryPred1AndIni(env, pred1)))
  .ret[[4]] <- bquote(rx_hi_ ~ .(rxode2::.rxGetHiBoundaryPred1AndIni(env, pred1)))
  .ret[[5]] <- bquote(rx_pred_f_ ~ .(rxode2::.rxGetPredictionF(env, pred1)))
  .ret[[6]] <- bquote(rx_pred_ ~ .(rxode2::.rxGetPredictionFTransform(env, pred1, .yj)))
  .ret
}

#' @export
rxGetDistributionNlmeLines.t <- function(line) {
  stop("t isn't supported yet", call.=FALSE)
}

#' @export
rxGetDistributionNlmeLines.default  <- function(line) {
  stop("distribution not supported", call.=FALSE)
}

#' @export
rxGetDistributionNlmeLines.rxUi <- function(line) {
  .predDf <- get("predDf", line)
  if (length(.predDf$cond) != 1) stop("nlme does not support multiple endpoint models", call.=FALSE)
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
rxUiGet.nlmeFD <- function(x, ...) {
  .loadSymengine(.nlmePrune(x), promoteLinSens = FALSE)
}
#attr(rxUiGet.nlmeFD, "desc") <- "Symengine environment pred-only nlme environment"

#' @export
rxUiGet.nlmeFunction <- function(x, ...) {
  .ui <- x[[1]]
  .estPar <- rxUiGet.saemParamsToEstimate(x, ...)
  #.par <- c(.estPar, .ui$covariates)
  .par <- .estPar
  eval(parse(text=paste0("function(", paste(.par, collapse=","), ", ID) {\n",
                         "nlmixr2est::.nlmixrNlmeFun(list(", paste(paste0(.estPar, "=", .estPar), collapse=","), "), ID)\n",
                         "}")))



}

#' @export
rxUiGet.nlmeRxModelFD <- function(x, ...) {
  .s <- rxUiGet.nlmeFD(x, ...)
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  ## .lhs0 <- .s$..lhs0
  ## if (is.null(.lhs0)) .lhs0 <- ""
  .ddt <- .s$..ddt
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

#' @export
rxUiGet.nlmeRxModel <- function(x, ...) {
  rxUiGet.nlmeRxModelFD(x, ...)
}


#attr(rxUiGet.nlmeRxModel, "desc") <- "nlme rxode2 text"

#' @export
rxUiGet.nlmeModel <- function(x, ...) {
  .ui <- x[[1]]
  # This includes both fixed effects and non-mu-referenced random effects
  saem.par <- rxUiGet.saemParamsToEstimate(x, ...)
  muRef.par <-
    setNames(
      paste(x[[1]]$muRefDataFrame$theta, x[[1]]$muRefDataFrame$eta, sep="+"),
      x[[1]]$muRefDataFrame$theta
    )
  nonMuRef.par <-
    setNames(nm=setdiff(saem.par, x[[1]]$muRefDataFrame$theta))
  all.par <- c(muRef.par, nonMuRef.par)
  as.formula(
    sprintf(
      "DV~.nlmixrNlmeFun(pars=list(%s), id=ID)",
      paste(names(all.par), all.par, sep="=", collapse=", ")
    )
  )
}
#attr(rxUiGet.nlmeModel, "desc") <- "nlme formula for nlmixr model"

#' @export
rxUiGet.nlmeGradDimnames <- function(x, ...) {
  .estPar <- rxUiGet.saemParamsToEstimate(x, ...)
  eval(parse(text=paste0("list(NULL, list(", paste(paste0(.estPar, "=quote(", .estPar, ")"), collapse=","), "))")))
}


#' @export
rxUiGet.nlmePdOmega <- function(x, ...) {
  .ui <- x[[1]]
  .omega <- .ui$omega
  .omega2 <- .omega
  diag(.omega2) <- 0
  .name <- dimnames(.omega)[[1]]
  .name <- .nlmeGetNonMuRefNames(.name, .ui)
  dimnames(.omega) <- list(.name, .name)
  if (all(.omega2 == 0)) {
    nlme::pdDiag(value=.omega, form=as.formula(paste(paste(.name, collapse="+"), "~1")))
  } else {
    .omega <- as.matrix(Matrix::nearPD(.omega)$mat)
    dimnames(.omega) <- list(.name, .name)
    warning("nlme will estimate a full omega matrix if any covariances are estimated", call.=FALSE)
    nlme::pdSymm(value=.omega, form=as.formula(paste(paste(.name, collapse="+"), "~1")))
  }
}
#attr(rxUiGet.nlmePdOmega, "desc") <- "nlme omega matrix form"


#' @export
rxUiGet.nlmeStart <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .w <- which(!is.na(.iniDf$ntheta) & is.na(.iniDf$err))
  setNames(.iniDf$est[.w], .iniDf$name[.w])
}
#attr(rxUiGet.nlmeStart, "desc") <- "nlme starting estimates for fixed effects"


#' @export
rxUiGet.nlmeFixedFormula <- function(x, ...) {
  .start <- rxUiGet.nlmeStart(x, ...)
  as.formula(paste(paste(names(.start), collapse="+"), "~1"))
}
#attr(rxUiGet.nlmeStart, "desc") <- "nlme starting estimates for fixed effects"
#' @export
rxUiGet.nlmeWeights <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  if (length(.predDf$cond) != 1) stop("cannot apply to multiple endpoint models", call.=FALSE)
  if (.predDf$distribution != "norm") stop("nlme not supported for this unexplained error structure", call.=FALSE)
  .errType <- .predDf$errType
  if (.errType == "prop") {
    return(nlme::varPower(fixed=list(power=1)))
  } else if (.errType == "pow") {
    return(nlme::varPower())
  } else if (.errType == "add") {
    return(NULL)
  }
  .addProp <- .predDf$addProp
  if (.addProp == "default") {
    .addProp <- rxode2::rxGetControl(.ui, "addProp", "combined2")
  }
  if (.addProp == "combined1") {
    if (.errType == "add + prop") {
      return(nlme::varConstPower(fixed=list(power=1)))
    } else {
      return(nlme::varConstPower())
    }
  } else {
    if (.errType == "add + prop") {
      return(nlme::varConstProp())
    } else {
      stop("add+prop combined2 does not support nlme power currently",
           call.=FALSE)
    }
  }
}

