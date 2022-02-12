#' Control Values for nlme Fit with extra options for nlmixr
#'
#' The values supplied in the function call replace the defaults and
#' a list with all possible arguments is returned.  The returned list
#' is used as the ‘control’ argument to the ‘nlme’ function.
#'
#' @inheritParams nlme::nlmeControl
#' @inheritParams nlme::nlme
#' @param sens Calculate gradients using forward sensitivity
#' @param returnNlme Returns the nlme object instead of the nlmixr
#'   object (by default FALSE).  If any of the nlme specific options
#'   of `random`, `fixed`, `sens`, the nlme object is returned
#' @inheritParams foceiControl
#' @export
#' @examples
#' nlmixr2::nlmeControl()
#' nlmixr2NlmeControl()
nlmixr2NlmeControl <- function(maxIter = 50, pnlsMaxIter = 7, msMaxIter = 50, minScale = 0.001,
    tolerance = 1e-05, niterEM = 25, pnlsTol = 0.001, msTol = 1e-06,
    returnObject = TRUE, msVerbose = FALSE, msWarnNoConv = TRUE,
    gradHess = TRUE, apVar = TRUE, .relStep = .Machine$double.eps^(1/3),
    minAbsParApVar = 0.05, opt = c("nlminb", "nlm"), natural = TRUE,
    sigma = NULL, optExpression=TRUE, sumProd=FALSE,
    rxControl=rxode2::rxControl(atol=1e-4, rtol=1e-4),
    method=c("ML", "REML"),
    random=NULL, fixed=NULL, weights=NULL, sens=FALSE, verbose=TRUE, returnNlme=FALSE,
    addProp = c("combined2", "combined1"), calcTables=TRUE, compress=TRUE,
    adjObf=TRUE, ...) {

  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnObject, len=1, any.missing=FALSE)
  checkmate::assertLogical(msVerbose, len=1, any.missing=FALSE)
  checkmate::assertLogical(msWarnNoConv, len=1, any.missing=FALSE)
  checkmate::assertLogical(gradHess, len=1, any.missing=FALSE)
  checkmate::assertLogical(apVar, len=1, any.missing=FALSE)
  checkmate::assertLogical(natural, len=1, any.missing=FALSE)
  checkmate::assertLogical(sens, len=1, any.missing=FALSE)
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
  method <- match.arg(method)
  addProp <- match.arg(addProp)

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
               rxControl=rxControl, sens=sens, method=method,verbose=verbose,
               returnNlme=returnNlme, addProp=addProp, calcTables=calcTables,
               compress=compress, random=random, fixed=fixed, weights=weights,
               ...)
  class(.ret) <- "nlmeControl"
  .ret
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
  eval(bquote(nlme::nlme(model=.(ui$nlmeModel), data=.nlmeFitDataObservations,
                         method=.(.method),fixed=.(.fixed), random=.(.random),
                         start=.(ui$nlmeStart), weights=.(.weights),
                         control=.(.ctl), verbose=.(.verbose), na.action=function(object, ...) {
                           return(object)
                         })))
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
      return(c(.f, .add, .prop))
    } else {
      .nlmePars <- coef(nlme$modelStruct$varStruct)
      .w <- which(ui$iniDf$err == "add")
      .add <- setNames(exp(.nlmePars["const"]), ui$iniDf$name[.w])
      .w <- which(ui$iniDf$err == "pow")
      .prop <- setNames(nlme$sigma, ui$iniDf$name[.w])
      .w <- which(ui$iniDf$err == "pow2")
      .pow <- setNames(.nlmePars["power"], ui$iniDf$name[.w])
      return(c(.f, .add, .prop, .pow))
    }
  } else {
    if (.errType == "add + prop") {
      .nlmePars <- coef(nlme$modelStruct$varStruct)
      .w <- which(ui$iniDf$err == "add")
      .add <- setNames(exp(.nlmePars["const"]), ui$iniDf$name[.w])
      .w <- which(ui$iniDf$err == "prop")
      .prop <- setNames(.nlmePars["prop"], ui$iniDf$name[.w])
      return(c(.f, .add, .prop))
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
  .ome <- as.matrix(Matrix::nearPD(f$omega)$mat)
  dimnames(.ome) <- list(.name, .name)
  .ome
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.nlmeControl <- function(control, env) {
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

.nlmeControlToFoceiControl <- function(env) {
  .nlmeControl <- env$nlmeControl
  .ui <- env$ui
  .ctl <- env$nlmeControl$rxControl
  names(.ctl) <- sub("maxsteps", "maxstepsOde", names(.ctl))
  .ctl <- .ctl[names(.ctl) != "scale"]
  .ctl$maxOuterIterations <- 0
  .ctl$maxInnerIterations <- 0
  .ctl$covMethod <- "" #.covMethod
  .ctl$etaMat <- env$etaMat
  .ctl$sumProd <- .nlmeControl$sumProd
  .ctl$optExpression <- .nlmeControl$optExpression
  .ctl$scaleTo <- 0
  .ctl$calcTables <- .nlmeControl$calcTables
  if (.nlmeControl$addProp == 1L) {
    .ctl$addProp <- "combined1"
  } else {
    .ctl$addProp <- "combined2"
  }
  .ctl$skipCov <- .ui$foceiSkipCov
  .ctl$interaction <- 1L
  .ctl$compress <- .nlmeControl$compress
  env$control <- do.call(foceiControl, .ctl)
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
  # - $ui for ui object
  # - $fullTheta Full theta information
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
  .foceiPreProcessData(.data, .ret, .ui)
  .et <- rxode2::etTrans(.ret$dataSav, .ui$mv0, addCmt=TRUE)
  # Just like saem, nlme can use mu-referenced covariates
  .nTv <- attr(class(.et), ".rxode2.lst")$nTv
  if (is.null(.nTv)) .nTv <- 0
  .tv <- character(0)
  if (.nTv != 0) {
    .tv <- names(.et)[-seq(1, 6)]
  }
  if (rxode2::rxGetControl(.ui, "sens", FALSE) && length(f$nonMuEtas) > 0) {
    stop("'sens=TRUE' requires mu-referenced etas", call.=FALSE)
  }
  .nlme <- .collectWarnings(.nlmeFitModel(.ui, .ret$dataSav, timeVaryingCovariates=.tv), lst = TRUE)
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
  .ui <- env$ui
  .nlmeFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)})
  .nlmeFamilyFit(env,  ...)
}


