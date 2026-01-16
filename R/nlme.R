#' Control Values for nlme Fit with extra options for nlmixr
#'
#' The values supplied in the function call replace the defaults and
#' a list with all possible arguments is returned.  The returned list
#' is used as the ‘control’ argument to the ‘nlme’ function.
#'
#' @inheritParams nlme::nlmeControl
#' @inheritParams nlme::nlme
#' @param returnNlme Returns the nlme object instead of the nlmixr
#'   object (by default FALSE).  If any of the nlme specific options
#'   of `random`, `fixed`, `sens`, the nlme object is returned
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @return a nlmixr-nlme list
#' @examples
#' nlmeControl()
#' nlmixr2NlmeControl()
#' @family Estimation control
#' @export
nlmixr2NlmeControl <- function(maxIter = 100, pnlsMaxIter = 100, msMaxIter = 100, minScale = 0.001,
    tolerance = 1e-05, niterEM = 25, pnlsTol = 0.001, msTol = 1e-06,
    returnObject = FALSE, msVerbose = FALSE, msWarnNoConv = TRUE,
    gradHess = TRUE, apVar = TRUE, .relStep = .Machine$double.eps^(1/3),
    minAbsParApVar = 0.05, opt = c("nlminb", "nlm"), natural = TRUE,
    sigma = NULL, optExpression=TRUE, literalFix=TRUE, sumProd=FALSE,
    rxControl=NULL,
    method=c("ML", "REML"),
    random=NULL, fixed=NULL, weights=NULL, verbose=TRUE, returnNlme=FALSE,
    addProp = c("combined2", "combined1"), calcTables=TRUE, compress=TRUE,
    adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, muRefCovAlg=TRUE, ...) {

  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnObject, len=1, any.missing=FALSE)
  checkmate::assertLogical(msVerbose, len=1, any.missing=FALSE)
  checkmate::assertLogical(msWarnNoConv, len=1, any.missing=FALSE)
  checkmate::assertLogical(gradHess, len=1, any.missing=FALSE)
  checkmate::assertLogical(apVar, len=1, any.missing=FALSE)
  checkmate::assertLogical(natural, len=1, any.missing=FALSE)
  checkmate::assertLogical(verbose, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnNlme, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)

  checkmate::assertIntegerish(pnlsMaxIter, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(msMaxIter, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(niterEM, len=1, any.missing=FALSE, lower=1)
  checkmate::assertNumeric(minScale, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(pnlsTol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(msTol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(tolerance, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(.relStep, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(minAbsParApVar, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(ci, lower=0, upper=1, any.missing=FALSE, len=1)
  checkmate::assertLogical(muRefCovAlg, any.missing=FALSE, len=1)

  method <- match.arg(method)
  addProp <- match.arg(addProp)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl", "covMethod"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (is.null(sigdig)) {
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

  if (is.null(sigma))
    sigma <- 0
  else if (!is.finite(sigma) || length(sigma) != 1 || sigma < 0)
    stop("Within-group std. dev. must be a positive numeric value")

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

  .ret <- list(maxIter = maxIter, pnlsMaxIter = pnlsMaxIter, msMaxIter = msMaxIter,
               minScale = minScale, tolerance = tolerance, niterEM = niterEM,
               pnlsTol = pnlsTol, msTol = msTol, returnObject = returnObject,
               msVerbose = msVerbose, msWarnNoConv = msWarnNoConv, gradHess = gradHess,
               apVar = apVar, .relStep = .relStep, minAbsParApVar = minAbsParApVar,
               opt = match.arg(opt), natural = natural, sigma = sigma,
               optExpression=optExpression, literalFix=literalFix, sumProd=sumProd,
               rxControl=rxControl, method=method,verbose=verbose,
               returnNlme=returnNlme, addProp=addProp, calcTables=calcTables,
               compress=compress, random=random, fixed=fixed, weights=weights,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable, muRefCovAlg=muRefCovAlg,
               genRxControl=.genRxControl)
  class(.ret) <- "nlmeControl"
  .ret
}

#' @export
rxUiDeparse.nlmeControl <- function(object, var) {
  .default <- nlmeControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}


#' @rdname nlmixr2NlmeControl
#' @export
nlmeControl <- nlmixr2NlmeControl
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
    .control <- nlmixr2est::nlmeControl()
  }
  if (!inherits(.control, "nlmeControl")){
    .control <- do.call(nlmixr2est::nlmeControl, .control)
  }
  assign("control", .control, envir=.ui)
}

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
    .datF <- nlmixr2global$nlmeFitDataAll[nlmixr2global$nlmeFitDataAll$ID == .ids[i], ]
    .datF$ID <- i
    .datF
  }))
  .pars <- as.data.frame(c(pars, list(ID=id)))
  .pars <- .pars[!duplicated(.pars$ID),]
  .pars$ID <- seq_along(.pars$ID)
  row.names(.pars) <- NULL
  .retF <- do.call(rxode2::rxSolve, c(list(object=nlmixr2global$nlmeFitRxModel, params=.pars, events=.datF),
                                      nlmixr2global$nlmeFitRxControl))
  .ret <- .retF$rx_pred_
  .ret
}

#' Setup the data for nlme estimation
#'
#' @param dataSav Formatted Data
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmeFitDataSetup <- function(dataSav) {
  .dsAll <- dataSav[dataSav$EVID != 2, ] # Drop EVID=2 for estimation
  nlmixr2global$nlmeFitDataAll <- .dsAll
}

.nlmeFitModel <- function(ui, dataSav, timeVaryingCovariates) {
  .nlmeFitDataSetup(dataSav)
  nlmixr2global$nlmeFitRxModel <- rxode2::rxode2(ui$nlmeRxModel)
  nlmixr2global$nlmeFitRxControl <- rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl())

  .ctl <- ui$control
  class(.ctl) <- NULL
  .fixed <- rxode2::rxGetControl(ui, "fixed", NULL)
  if (is.null(.fixed)) {
    .fixed <- ui$nlmeFixedFormula
  } else {
    rxode2::rxAssignControlValue(ui, "returnNlme", TRUE)
  }
  .random <- rxode2::rxGetControl(ui, "random", NULL)
  if (is.null(.random)) {
    .random <- ui$nlmePdOmega
  } else {
    rxode2::rxAssignControlValue(ui, "returnNlme", TRUE)
  }
  .verbose <- rxode2::rxGetControl(ui, "verbose", TRUE)
  .method <- rxode2::rxGetControl(ui, "method", "ML")
  .weights <- rxode2::rxGetControl(ui, "weights", NULL)
  if (is.null(.weights)) {
    .weights <- ui$nlmeWeights
  } else {
    rxode2::rxAssignControlValue(ui, "returnNlme", TRUE)
  }
  ret <-
    eval(bquote(nlme::nlme(
      model=.(ui$nlmeModel),
      data=nlme::groupedData(DV ~ TIME | ID, dataSav[dataSav$EVID == 0, ]),
      method=.(.method),
      fixed=.(.fixed),
      random=.(.random),
      start=.(ui$nlmeStart),
      weights=.(.weights),
      control=.(.ctl),
      verbose=.(.verbose),
      na.action=function(object, ...) {
        object
      }
    )))
  ret
}

#' Get the theta estimates from nlme using roxde2 ui
#'
#' @param nlme nlme object
#' @param ui rxode2 object
#' @return named theta vector
#' @author Matthew L. Fidler
#' @noRd
.nlmeGetTheta <- function(nlme, ui) {
  .f <- nlme::fixef(nlme)
  .predDf <- ui$predDf
  .iniDf <- ui$iniDf
  .errType <- .predDf$errType
  if (.errType == "prop") {
    .w <- which(ui$iniDf$err == "prop")
    .prop <- setNames(nlme$sigma, ui$iniDf$name[.w])
    return(c(.f, .prop))
  } else if (.errType == "pow") {
    .nlmePars <- coef(nlme$modelStruct$varStruct)
    .w <- which(ui$iniDf$err == "pow2")
    .pw <- setNames(.nlmePars["power"], ui$iniDf$name[.w])
    .w <- which(ui$iniDf$err == "pow")
    .coef <- setNames(nlme$sigma, ui$iniDf$name[.w])
    return(c(.f, .coef, .pw))
  } else if (.errType == "add") {
    .w <- which(ui$iniDf$err == "add")
    .add <- setNames(nlme$sigma, ui$iniDf$name[.w])
    return(c(.f, .add))
  }
  .addProp <- .predDf$addProp
  if (.addProp == "default") {
    .addProp <- rxode2::rxGetControl(ui, "addProp", "combined2")
  }
  if (.addProp == "combined1") {
    if (.errType == "add + prop") {
       .nlmePars <- coef(nlme$modelStruct$varStruct)
      .w <- which(ui$iniDf$err == "add")
      .add <- setNames(exp(.nlmePars["const"]), ui$iniDf$name[.w])
      .w <- which(ui$iniDf$err == "prop")
      .prop <- setNames(nlme$sigma, ui$iniDf$name[.w])
      c(.f, .add, .prop)
    } else {
      .nlmePars <- coef(nlme$modelStruct$varStruct)
      .w <- which(ui$iniDf$err == "add")
      .add <- setNames(exp(.nlmePars["const"]), ui$iniDf$name[.w])
      .w <- which(ui$iniDf$err == "pow")
      .prop <- setNames(nlme$sigma, ui$iniDf$name[.w])
      .w <- which(ui$iniDf$err == "pow2")
      .pow <- setNames(.nlmePars["power"], ui$iniDf$name[.w])
      c(.f, .add, .prop, .pow)
    }
  } else {
    if (.errType == "add + prop") {
      .nlmePars <- coef(nlme$modelStruct$varStruct)
      .w <- which(ui$iniDf$err == "add")
      .add <- setNames(exp(.nlmePars["const"]), ui$iniDf$name[.w])
      .w <- which(ui$iniDf$err == "prop")
      .prop <- setNames(.nlmePars["prop"], ui$iniDf$name[.w])
      c(.f, .add, .prop)
    } else {
      stop("add+prop combined2 does not support nlme power currently",
           call.=FALSE)
    }
  }
}

#' Get non mu referenced names from mu referenced theta
#'
#' @param names Names to translate
#' @param ui rxode2 ui
#' @return non mu referenced names
#' @author Matthew L. Fidler
#' @noRd
.nlmeGetNonMuRefNames <- function(names, ui) {
  .muRef <- ui$muRefDataFrame
  vapply(names, function(n) {
    .w <- which(.muRef$theta == n)
    if (length(.w) == 1) return(.muRef$eta[.w])
    n
  }, character(1), USE.NAMES=FALSE)
}

#' Get the nlme eta matrix as expected for focei
#'
#' @param nlme nlme object
#' @param ui rxode2 ui object
#' @return matrix of eta estimates, ordered by ID and named by the eta names in focei.
#' @author Matthew L. Fidler
#' @noRd
.nlmeGetEtaMat <- function(nlme, ui) {
  .etaMat <- nlme::ranef(nlme)
  .etaMat <- .etaMat[order(as.numeric(row.names(.etaMat))),, drop = FALSE]
  names(.etaMat) <- .nlmeGetNonMuRefNames(names(.etaMat), ui)
  row.names(.etaMat) <- NULL
  as.matrix(.etaMat)
}

#' Get the covariance from nlme
#'
#' @param nlme nlme object
#' @author Matthew L. Fidler
#' @noRd
.nlmeGetCov <- function(nlme) {
  .snt <- summary(nlme)$tTable
  .se <- .snt[,"Std.Error"]
  if (length(.se) == 1) {
    matrix(.se * .se, 1, 1, dimnames=list(rownames(.snt), rownames(.snt)))
  } else {
    .cov <- diag(.se * .se)
    dimnames(.cov) <- list(rownames(.snt), rownames(.snt))
    .cov
  }
}

#' Get the omega matrix from nlme
#'
#' @param nlme nlme object
#' @param ui rxode2 object
#' @return Named omega matrix
#' @author Matthew L. Fidler
#' @noRd
.nlmeGetOmega <- function(nlme, ui) {
  .omega <- ui$omega
  diag(.omega) <- 0
  .vc <- nlme::VarCorr(nlme)
  .var <- as.matrix(.vc[,"Variance", drop = FALSE])
  .rn <- rownames(.var)
  .name <- .nlmeGetNonMuRefNames(.rn, ui)
  .var <- setNames(suppressWarnings(as.numeric(.var)), .name)
  .var <- .var[names(.var) != "Residual"]
  if (length(.var) == 1) {
    .ome <- matrix(.var, 1, 1)
  } else {
   .ome <- diag(.var)
  }
  .name <- names(.var)
  dimnames(.ome) <- list(.name, .name)
  if (all(.omega == 0)) {
    return(.ome)
  }
  .cor2 <- as.data.frame(.vc[-length(.rn), -(1:2), drop = FALSE])
  .cor2$extra <- ""
  names(.cor2) <- rownames(.cor2)
  .cor2 <- as.matrix(.cor2)
  diag(.cor2) <- "1"
  .cor2[upper.tri(.cor2)] <- .cor2[lower.tri(.cor2)]
  .cor2 <- matrix(suppressMessages(as.numeric(.cor2)), nrow(.cor2), ncol(.cor2), dimnames=dimnames(.ome))
  diag(.ome) <- sqrt(diag(.ome))
  .ome <- .ome %*% .cor2 %*% .ome
  .ome <- as.matrix(Matrix::nearPD(ui$omega)$mat)
  dimnames(.ome) <- list(.name, .name)
  .ome
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.nlmeControl <- function(control, env) {
  eval(rxode2::rxUiDeparse(control, "control"))
  assign("nlmeControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.nlme <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nlmeControl", .env)) {
    .control <- get("nlmeControl", .env)
    if (inherits(.control, "nlmeControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nlmeControl")) return(.control)
  }
  stop("cannot find nlme related control object", call.=FALSE)
}

.nlmeControlToFoceiControl <- function(env, assign=TRUE) {
  .nlmeControl <- env$nlmeControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$nlmeControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                etaMat=env$etaMat,
                                sumProd=.nlmeControl$sumProd,
                                optExpression=.nlmeControl$optExpression,
                                literalFix=.nlmeControl$literalFix,
                                literalFixRes=FALSE,
                                scaleTo=0,
                                calcTables=.nlmeControl$calcTables,
                                addProp=.nlmeControl$addProp,
                                skipCov=.ui$foceiSkipCov,
                                interaction=1L,
                                compress=.nlmeControl$compress,
                                ci=.nlmeControl$ci,
                                sigdigTable=.nlmeControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @export
#' @rdname nmObjGetFoceiControl
nmObjGetFoceiControl.nlme <- function(x, ...) {
  .nlmeControlToFoceiControl(x[[1]])
}

.nlmeFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  # The environment needs:
  # - table for table options
  # - $origData -- Original Data
  # - $dataSav -- Processed data from .foceiPreProcessData
  # - $idLvl -- Level information for ID factor added
  # - $covLvl -- Level information for items to convert to factor
  # - $ui for ui fullTheta Full theta information
  # - $etaObf data frame with ID, etas and OBJI
  # - $cov For covariance
  # - $covMethod for the method of calculating the covariance
  # - $adjObf Should the objective function value be adjusted
  # - $objective objective function value
  # - $extra Extra print information
  # - $method Estimation method (for printing)
  # - $omega Omega matrix
  # - $etaObf Eat objective function data frame
  # - $theta Is a theta data frame
  # - $model a list of model information for table generation.  Needs a `predOnly` model
  # - $message Message for display
  # - $est estimation method
  # - $ofvType (optional) tells the type of ofv is currently being used
  # When running the focei problem to create the nlmixr object, you also need a
  #  foceiControl object
  .ret$table <- env$table
  .foceiPreProcessData(.data, .ret, .ui, .control$rxControl)
  .et <- rxode2::etTrans(.ret$dataSav, .ui$mv0, addCmt=TRUE,
                         addlKeepsCov = .control$rxControl$addlKeepsCov,
                         addlDropSs = .control$rxControl$addlDropSs,
                         ssAtDoseTime = .control$rxControl$ssAtDoseTime)
  # Just like saem, nlme can use mu-referenced covariates
  .nTv <- attr(class(.et), ".rxode2.lst")$nTv
  if (is.null(.nTv)) .nTv <- 0
  .tv <- character(0)
  if (.nTv != 0) {
    .tv <- names(.et)[-seq(1, 6)]
  }
  .nlme <- .collectWarn(.nlmeFitModel(.ui, .ret$dataSav, timeVaryingCovariates=.tv), lst = TRUE)
  .ret$nlme <- .nlme[[1]]
  .ret$message <- NULL
  lapply(.nlme[[2]], function(x){
    warning(x, call.=FALSE)
    if (regexpr("PNLS", x) != -1) {
      .ret$message <- c(.ret$message, paste0(x, " (carefully review results)"))
    }
  })
  if (is.null(.ret$message)) {
    .ret$message <- ""
  } else {
    .ret$message <- paste(.ret$message, collapse="\n")
  }
  if (rxode2::rxGetControl(.ui, "returnNlme", FALSE)) {
    return(.ret$nlme)
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .nlmeGetTheta(.ret$nlme, .ui)
  .ret$cov <- .nlmeGetCov(.ret$nlme)
  .ret$covMethod <- "nlme"
  .ret$etaMat <- .nlmeGetEtaMat(.ret$nlme, .ui)
  .ret$etaObf <- data.frame(ID = seq_along(.ret$etaMat[, 1]),
                           as.data.frame(.ret$etaMat),
                           OBJI = NA)
  .ret$omega <- .nlmeGetOmega(.ret$nlme, .ui)
  .ret$control <- .control
  .ret$extra <- paste0(" by ", crayon::bold$yellow(ifelse(.control$method == "REML", "REML", "maximum likelihood")))
  .nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui)) {
    rm(list="control", envir=.ui)
  }
  .ret$est <- "nlme"
  # There is no parameter history for nlme
  .ret$objective <- -2 * as.numeric(logLik(.ret$nlme))
  .ret$model <- .ui$ebe
  .ret$est <- "nlme"
  .ret$ofvType <- "nlme"
  .nlmeControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="nlme")
  .env <- .ret$env
  .env$method <- "nlme"
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.nlme <- function(env, ...) {
  .model <- .uiApplyMu2(env)
  .ui <- env$ui
  rxode2::assertRxUiMixedOnly(.ui, " for the estimation routine 'nlme', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiNormal(.ui, " for the estimation routine 'nlme'", .var.name=.ui$modelName)
  rxode2::assertRxUiSingleEndpoint(.ui, " for the estimation routine 'nlme'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'nlme'", .var.name=.ui$modelName)
  rxode2::assertRxUiEstimatedResiduals(.ui, " for the estimation routine 'nlme'", .var.name=.ui$modelName)
  .nlmeFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .uiFinalizeMu2(.nlmeFamilyFit(env,  ...), .model)
}
attr(nlmixr2Est.nlme, "covPresent") <- TRUE
