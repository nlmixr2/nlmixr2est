#' nlmixr2 defaults controls for nlm
#'
#' @inheritParams stats::nlm
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @return nlme control object
#' @details
#'
#' Note the covariance is calculated by nlmixr instead of optimHess, so `hessian` is not a possible option
#'
#' @export
#' @author Matthew L. Fidler
#' @examples
nlmControl <- function(typsize = NULL,
                       fscale = 1, print.level = 2, ndigit = NULL, gradtol = 1e-6,
                       stepmax = NULL,
                       steptol = 1e-6, iterlim = 10000, check.analyticals = FALSE,
                       rxControl=NULL,
                       optExpression=TRUE, sumProd=FALSE,
                       returnNlm=FALSE,
                       addProp = c("combined2", "combined1"),
                       calcTables=TRUE, compress=TRUE,
                       covMethod=c("r", "nlm", ""),
                       adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {
  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
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
  covMethod <- match.arg(covMethod)
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

  .ret <- list(covMethod=covMethod,
               typsize = typsize,
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


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.nlmControl <- function(control, env) {
  assign("nlmControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.nlm <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nlmControl", .env)) {
    .control <- get("nlmControl", .env)
    if (inherits(.control, "nlmControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nlmControl")) return(.control)
  }
  stop("cannot find nlm related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.nlm <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- nlmControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("nlmControl", .ctl)
  if (!inherits(.ctl, "nlmControl")) {
    .minfo("invalid control for `est=\"nlm\"`, using default")
    .ctl <- nlmControl()
  } else {
    .ctl <- do.call(nlmControl, .ctl)
  }
  .ctl
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
  .retF <- do.call(rxode2::rxSolve,
                   c(list(object=.nlmEnv$model,
                          params=.nlmEnv$parTrans(pars),
                          events=.nlmEnv$data),
                     .nlmEnv$rxControl))
  -sum(.retF$rx_pred_)
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
               if (!is.na(.iniDf$err[t]) &&
                     !is.finite(.iniDf$upper[t]) &&
                      .iniDf$lower[t] == 0) {
                 .ret <- paste0("'THETA[", t, "]'=abs(p[", .env$i, "])")
               } else {
                 .ret <- paste0("'THETA[", t, "]'=p[", .env$i, "]")
                 .env$i <- .env$i + 1
               }
               .ret
             }
           }, character(1), USE.NAMES=FALSE), collapse=","), ")}")))
}

#' @export
rxUiGet.nlmParIni <- function(x, ...) {
  .ui <- x[[1]]
  .ui$iniDf$est[!.ui$iniDf$fix]
}

#' @export
rxUiGet.nlmParName <- function(x, ...) {
  .ui <- x[[1]]
  .ui$iniDf$name[!.ui$iniDf$fix]
}

#' Setup the data for nlm estimation
#'
#' @param dataSav Formatted Data
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmFitDataSetup <- function(dataSav) {
  .dsAll <- dataSav[dataSav$EVID != 2, ] # Drop EVID=2 for estimation
  .nlmEnv$data <- rxode2::etTrans(.dsAll, .nlmEnv$model)
}

.nlmFitModel <- function(ui, dataSav) {
  .nlmEnv$model <- rxode2::rxode2(ui$nlmRxModel)
  .nlmFitDataSetup(dataSav)
  .nlmEnv$rxControl <- rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl())
  .nlmEnv$rxControl$returnType <- 2L # use data.frame output
  .nlmEnv$parTrans <- ui$nlmParNameFun
  .ctl <- ui$control
  class(.ctl) <- NULL
  .p <- ui$nlmParIni
  .typsize <- .ctl$typsize
  if (is.null(.typsize)) {
    .typsize <- rep(1, length(.p))
  } else if (length(.typsize) == 1L) {
    .typsize <- rep(.typsize, length(.p))
  } else {
    stop("'typsize' needs to match the number of estimated parameters (or equal 1)", call.=FALSE)
  }
  .stepmax <- .ctl$stepmax
  if (is.null(.stepmax)) {
    .stepmax <- max(1000 * sqrt(sum((.p/.typsize)^2)), 1000)
  }
  .ret <- eval(bquote(stats::nlm(
    f=.(nlmixr2est::.nlmixrNlmFun),
    p=.(.p),
    hessian=.(.ctl$covMethod == "nlm"),
    typsize=.(.typsize),
    fscale=.(.ctl$fscale),
    print.level=.(.ctl$print.level),
    ndigit=.(.ctl$ndigit),
    gradtol=.(.ctl$gradtol),
    stepmax=.(.stepmax),
    steptol = .(.ctl$steptol),
    iterlim = .(.ctl$iterlim),
    check.analyticals = .(.ctl$check.analyticals)
  )))
  # be nice and name items
  .name <- ui$nlmParName
  names(.ret$estimate) <- .name
  names(.ret$gradient) <- .name
  if (any(.ctl$covMethod == c("r", "nlm"))) {
    .malert("calculating covariance")
    if (.ctl$covMethod != "nlm") {
      .p <- setNames(.ret$estimate, NULL)
      .hess <- nlmixr2Hess(.p, nlmixr2est::.nlmixrNlmFun)
      .ret$hessian <- .hess
    }
    dimnames(.ret$hessian) <- list(.name, .name)
    .hess <- .ret$hessian

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
      if (.ctl$covMethod == "nlm") {
        .ret$covMethod <- paste0(.covType, " (nlm)")
      } else {
        .ret$covMethod <- .covType
      }
      .ret$cov <- .cov

    } else {
      .ret$covMethod <- "failed"
    }
    dimnames(.r) <- list(.name, .name)
    .ret$r <- .r
    .msuccess("done")
  }
  .ret
}
#' Get the full theta for nlm methods
#'
#' @param nlm enhanced nlm return
#' @param ui ui object
#' @return named theta matrix
#' @author Matthew L. Fidler
#' @noRd
.nlmGetTheta <- function(nlm, ui) {
  .iniDf <- ui$iniDf
  setNames(vapply(seq_along(.iniDf$name),
         function(i) {
           if (.iniDf$fix[i]) {
             return(.iniDf$est[i])
           } else {
             return(nlm$estimate[.iniDf$name[i]])
           }
         }, double(1), USE.NAMES=FALSE),
         .iniDf$name)
}

.nlmControlToFoceiControl <- function(env, assign=TRUE) {
  .nlmControl <- env$nlmControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$nlmControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.nlmControl$sumProd,
                                optExpression=.nlmControl$optExpression,
                                scaleTo=0,
                                calcTables=.nlmControl$calcTables,
                                addProp=.nlmControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.nlmControl$compress,
                                ci=.nlmControl$ci,
                                sigdigTable=.nlmControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}


.nlmFamilyFit <- function(env, ...) {
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
  # - $theta Is a theta data frame
  # - $model a list of model information for table generation.  Needs a `predOnly` model
  # - $message Message for display
  # - $est estimation method
  # - $ofvType (optional) tells the type of ofv is currently being used
  # When running the focei problem to create the nlmixr object, you also need a
  #  foceiControl object
  .ret$table <- env$table
  .foceiPreProcessData(.data, .ret, .ui, .control$rxControl)
  .nlm <- .collectWarn(.nlmFitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$nlm <- .nlm[[1]]
  .ret$message <- NULL
  if (rxode2::rxGetControl(.ui, "returnNlm", FALSE)) {
    return(.ret$nlm)
  }
  if (.ret$nlm$code == 1) {
    .ret$message <- "relative gradient is close to zero, current iterate is probably solution"
  } else if (.ret$nlm$code == 2) {
    .ret$message <- "successive iterates within tolerance, current iterate is probably solution"
  } else if (.ret$nlm$code == 3) {
    .ret$message <- c("last global step failed to locate a point lower than 'estimate'",
                      "either 'estimate' is an approximate local minimum of the function or 'steptol' is too small")
  } else if (.ret$nlm$code == 4) {
    .ret$message <- "iteration limit exceeded"
  } else if (.ret$nlm$code == 5) {
    .ret$message <- c("maximum step size 'stepmax' exceeded five consecutive times",
                      "either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or 'stepmax' is too small")
  } else {
    .ret$message <- ""
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .nlmGetTheta(.ret$nlm, .ui)
  .ret$cov <- .ret$nlm$cov
  .ret$covMethod <- .ret$nlm$covMethod
  #.ret$etaMat <- NULL
  #.ret$etaObf <- NULL
  #.ret$omega <- NULL
  .ret$control <- .control
  .ret$extra <- ""
  .nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui)) {
    rm(list="control", envir=.ui)
  }
  .ret$est <- "nlm"
  # There is no parameter history for nlme
  .ret$objective <- -2 * as.numeric(.ret$nlm$minimum)
  .ret$model <- .ui$ebe
  .ret$ofvType <- "nlm"
  .nlmControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="nlm")
  .env <- .ret$env
  .env$method <- "nlm"
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.nlm <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'nlm', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'nlm'", .var.name=.ui$modelName)
  .nlmFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .nlmFamilyFit(env,  ...)
}
