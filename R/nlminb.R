#' nlmixr2 nlminb defaults
#'
#' @inheritParams stats::nlminb
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @param returnNlminb logical; when TRUE this will return the nlminb
#'   result instead of the nlmixr2 fit object
#' @param eval.max Maximum number of evaluations of the objective
#'   function allowed.  Defaults to 200.
#' @param iter.max Maximum number of iterations allowed.  Defaults to
#'   150.
#' @param trace The value of the objective function and the parameters
#'   is printed every trace'th iteration.  When 0 no trace information
#'   is to be printed
#'
#' @param abs.tol Absolute tolerance.  Defaults to 0 so the absolute
#'   convergence test is not used.  If the objective function is known
#'   to be non-negative, the previous default of `1e-20` would be more
#'   appropriate
#'
#' @param rel.tol Relative tolerance.  Defaults to `1e-10`.
#'
#' @param x.tol X tolerance.  Defaults to `1.5e-8`.
#'
#' @param xf.tol false convergence tolerance.  Defaults to `2.2e-14`.
#'
#' @param step.min  Minimum step size.  Default to ‘1.’.
#'
#' @param step.max Maximum step size.  Default to ‘1.’.
#'
#' @param sing.tol singular convergence tolerance; defaults to `rel.tol;.
#'
#' @param scale.init ... probably need to check PORT documentation
#'
#' @param diff.g an estimated bound on the relative error in the
#'   objective function value
#'
#' @inheritParams nlmControl
#'
#' @export
#' @author Matthew L. Fidler
#' @examples
#' \donttest{
#' # A logit regression example with emax model
#'
#' dsn <- data.frame(i=1:1000)
#' dsn$time <- exp(rnorm(1000))
#' dsn$DV=rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))
#'
#' mod <- function() {
#'  ini({
#'    E0 <- 0.5
#'    Em <- 0.5
#'    E50 <- 2
#'    g <- fix(2)
#'  })
#'  model({
#'    v <- E0+Em*time^g/(E50^g+time^g)
#'    ll(bin) ~ DV * v - log(1 + exp(v))
#'  })
#' }
#'
#' fit2 <- nlmixr(mod, dsn, est="nlminb")
#'
#' print(fit2)
#'
#' # you can also get the nlm output with fit2$nlminb
#'
#' fit2$nlminb
#' }
nlminbControl <- function(eval.max=200,
                          iter.max=150,
                          trace=0, # nolint
                          abs.tol=0,
                          rel.tol=1e-10,
                          x.tol=1.5e-8,
                          xf.tol=2.2e-14,
                          step.min=1,
                          step.max=1,
                          sing.tol=rel.tol,
                          scale = 1,
                          scale.init=NULL,
                          diff.g=NULL,
                          rxControl=NULL,
                          optExpression=TRUE, sumProd=FALSE,
                          literalFix=TRUE,
                          literalFixRes=TRUE,
                          returnNlminb=FALSE,
                          solveType=c("hessian", "grad", "fun"),

                          stickyRecalcN=4,
                          maxOdeRecalc=5,
                          odeRecalcFactor=10^(0.5),

                          eventType=c("central", "forward"),
                          shiErr=(.Machine$double.eps)^(1/3),
                          shi21maxFD=20L,

                          optimHessType=c("central", "forward"),
                          hessErr =(.Machine$double.eps)^(1/3),
                          shi21maxHess=20L,

                          useColor = crayon::has_color(),
                          printNcol = floor((getOption("width") - 23) / 12), #
                          print = 1L, #
                          normType = c("rescale2", "mean", "rescale", "std", "len", "constant"), #
                          scaleType = c("nlmixr2", "norm", "mult", "multAdd"), #
                          scaleCmax = 1e5, #
                          scaleCmin = 1e-5, #
                          scaleC=NULL,
                          scaleTo=1.0,
                          gradTo=1.0,

                          addProp = c("combined2", "combined1"),
                          calcTables=TRUE, compress=TRUE,
                          covMethod=c("r", "nlminb", ""),
                          adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {
  checkmate::assertIntegerish(eval.max, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(iter.max, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(trace, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(rel.tol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(x.tol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(xf.tol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(step.min, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(step.max, len=1, any.missing=FALSE, lower=0)
  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFixRes, len=1, any.missing=FALSE)

  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnNlminb, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)

  .solveTypeIdx <- c("hessian" = 3L, "grad" = 2L, "fun" = 1L)
  if (checkmate::testIntegerish(solveType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    solveType <- as.integer(solveType)
  } else {
    solveType <- setNames(.solveTypeIdx[match.arg(solveType)], NULL)
  }

  if (missing(covMethod) && any(solveType == 2:3)) {
    covMethod <- "nlminb"
  } else {
    covMethod <- match.arg(covMethod)
  }
  if (covMethod == "nlminb" && !any(solveType == 2:3)) {
    warning("using the Hessian function used during nlminb optimization requires a hessian or gradient solving type\n",
            "switching to covMethod='r'")
    covMethod <- "r"
  }

  .eventTypeIdx <- c("central" =2L, "forward"=1L)
  if (checkmate::testIntegerish(eventType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    eventType <- as.integer(eventType)
  } else {
    eventType <- setNames(.eventTypeIdx[match.arg(eventType)], NULL)
  }

  .optimHessTypeIdx <- c("central" =2L, "forward"=1L)
  if (checkmate::testIntegerish(optimHessType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    optimHessType <- as.integer(optimHessType)
  } else {
    optimHessType <- setNames(.optimHessTypeIdx[match.arg(optimHessType)], NULL)
  }

  checkmate::assertNumeric(shiErr, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(hessErr, lower=0, any.missing=FALSE, len=1)

  checkmate::assertIntegerish(shi21maxFD, lower=1, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(shi21maxHess, lower=1, any.missing=FALSE, len=1)

  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% "genRxControl")]
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


  checkmate::assertLogical(useColor, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(print, len=1, lower=0, any.missing=FALSE)
  checkmate::assertIntegerish(printNcol, len=1, lower=1, any.missing=FALSE)
  if (checkmate::testIntegerish(scaleType, len=1, lower=1, upper=4, any.missing=FALSE)) {
    scaleType <- as.integer(scaleType)
  } else {
    .scaleTypeIdx <- c("norm" = 1L, "nlmixr2" = 2L, "mult" = 3L, "multAdd" = 4L)
    scaleType <- setNames(.scaleTypeIdx[match.arg(scaleType)], NULL)
  }

  .normTypeIdx <- c("rescale2" = 1L, "rescale" = 2L, "mean" = 3L, "std" = 4L, "len" = 5L, "constant" = 6L)
  if (checkmate::testIntegerish(normType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    normType <- as.integer(normType)
  } else {
    normType <- setNames(.normTypeIdx[match.arg(normType)], NULL)
  }
  checkmate::assertNumeric(scaleCmax, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(scaleCmin, lower=0, any.missing=FALSE, len=1)
  if (!is.null(scaleC)) {
    checkmate::assertNumeric(scaleC, lower=0, any.missing=FALSE)
  }
  checkmate::assertNumeric(scaleTo, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(gradTo, len=1, lower=0, any.missing=FALSE)

  .ret <- list(eval.max=eval.max,
               iter.max=iter.max,
               trace=trace, # nolint
               abs.tol=abs.tol,
               rel.tol=rel.tol,
               x.tol=x.tol,
               xf.tol=xf.tol,
               step.min=step.min,
               step.max=step.max,
               sing.tol=sing.tol,
               scale.init=scale.init,
               diff.g=diff.g,
               scale=scale,
               solveType=solveType,
               stickyRecalcN=as.integer(stickyRecalcN),
               maxOdeRecalc=as.integer(maxOdeRecalc),
               odeRecalcFactor=odeRecalcFactor,

               eventType=eventType,
               shiErr=shiErr,
               shi21maxFD=as.integer(shi21maxFD),

               optimHessType=optimHessType,
               hessErr=hessErr,
               shi21maxHess=as.integer(shi21maxHess),

               useColor=useColor,
               print=print,
               printNcol=printNcol,
               scaleType=scaleType,
               normType=normType,

               scaleCmax=scaleCmax,
               scaleCmin=scaleCmin,
               scaleC=scaleC,
               scaleTo=scaleTo,
               gradTo=gradTo,

               covMethod=covMethod,
               optExpression=optExpression,
               literalFix=literalFix,
               literalFixRes=literalFixRes,
               sumProd=sumProd,
               rxControl=rxControl,
               returnNlminb=returnNlminb,
               addProp=match.arg(addProp),
               calcTables=calcTables,
               compress=compress,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "nlminbControl"
  .ret
}

#' @export
rxUiDeparse.nlminbControl <- function(object, var) {
  .default <- nlminbControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' A surrogate function for nlminb to call for ode solving
#'
#' @param pars Parameters that will be estimated
#' @return Predictions
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.nlmixrNlminbFunC <- function(pars) {
  .Call(`_nlmixr2est_nlminbFunC`, pars, 1L)
}
#' @rdname dot-nlmixrNlminbFunC
#' @export
.nlmixrNlminbGradC <- function(pars) {
  .Call(`_nlmixr2est_nlminbFunC`, pars, 2L)
}
#' @rdname dot-nlmixrNlminbFunC
#' @export
.nlmixrNlminbHessC <- function(pars) {
  .Call(`_nlmixr2est_nlminbFunC`, pars, 3L)
}
#' Get the nlminb family control
#'
#' @param env nlminb nlminbization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlminbFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2est::nlminbControl()
  }
  if (!inherits(.control, "nlminbControl")){
    .control <- do.call(nlmixr2est::nlminbControl, .control)
  }
  assign("control", .control, envir=.ui)
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.nlminbControl <- function(control, env) {
  eval(rxode2::rxUiDeparse(control, "control"))
  assign("nlminbControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.nlminb <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nlminbControl", .env)) {
    .control <- get("nlminbControl", .env)
    if (inherits(.control, "nlminbControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nlminbControl")) return(.control)
  }
  stop("cannot find nlminb related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.nlminb <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- nlminbControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("nlminbControl", .ctl)
  if (!inherits(.ctl, "nlminbControl")) {
    .minfo("invalid control for `est=\"nlminb\"`, using default")
    .ctl <- nlminbControl()
  } else {
    .ctl <- do.call(nlminbControl, .ctl)
  }
  .ctl
}

#' Setup the data for nlminb estimation
#'
#' @param dataSav Formatted Data
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlminbFitDataSetup <- function(dataSav) {
  .dsAll <- dataSav[dataSav$EVID != 2, ] # Drop EVID=2 for estimation
  if (any(names(.dsAll) == "CENS")) {
    if (!all(.dsAll$CENS == 0)) {
      stop("'nlminb' does not work with censored data", call. =FALSE)
    }
  }
  nlmixr2global$nlmEnv$data <- rxode2::etTrans(.dsAll, nlmixr2global$nlmEnv$model)
}

.nlminbFitModel <- function(ui, dataSav) {
  # Use nlmEnv and function for DRY principle
  .ctl <- ui$control
  .keep <- c("eval.max", "iter.max", "trace", "abs.tol", "rel.tol","x.tol", "xf.tol",
             "step.min", "step.max", "sing.tol", "diff.g", "scale.init")
  if (is.null(.ctl$diff.g)) {
    .keep <- .keep[.keep != "diff.g"]
  }
  if (is.null(.ctl$scale.init)) {
    .keep <- .keep[.keep != "scale.init"]
  }
  .oCtl <- setNames(lapply(.keep, function(x) {.ctl[[x]]}), .keep)
  class(.ctl) <- NULL
  if (.ctl$solveType == 1L) {
    .mi <-  ui$nlmRxModel
  } else {
    .mi <- ui$nlmSensModel
  }
  .p <- setNames(ui$nlmParIni, ui$nlmParName)
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl,
                       lower=ui$optimParLower, ui$optimParUpper)
  on.exit({.nlmFreeEnv()})
  if (.ctl$solveType == 1L) {
    # pred only
    .ret <- bquote(stats::nlminb(
      start=.(.env$par.ini),
      objective=.(nlmixr2est::.nlmixrNlminbFunC),
      scale = .(.ctl$scale),
      control = .(.oCtl),
      lower=.(.env$lower),
      upper=.(.env$upper)))
  } else {
    # grad/hessian added
    if (.ctl$solveType == 2L) {
      # Gradient
      .ret <- bquote(stats::nlminb(
        start=.(.env$par.ini),
        objective=.(nlmixr2est::.nlmixrNlminbFunC),
        gradient=.(nlmixr2est::.nlmixrNlminbGradC),
        scale = .(.ctl$scale),
        control = .(.oCtl),
        lower=.(.env$lower),
        upper=.(.env$upper)))
    } else {
      # Gradient / Hessian
      .ret <- bquote(stats::nlminb(
        start=.(.env$par.ini),
        objective=.(nlmixr2est::.nlmixrNlminbFunC),
        gradient=.(nlmixr2est::.nlmixrNlminbGradC),
        hessian=.(nlmixr2est::.nlmixrNlminbHessC),
        scale = .(.ctl$scale),
        control = .(.oCtl),
        lower=.(.env$lower),
        upper=.(.env$upper)))

    }
  }
  .ret <- eval(.ret)
  .nlmFinalizeList(.env, .ret, par="par", printLine=TRUE,
                   hessianCov=TRUE)
}

#' Get the full theta for nlminb methods
#'
#' @param nlm enhanced nlminb return
#' @param ui ui object
#' @return named theta matrix
#' @author Matthew L. Fidler
#' @noRd
.nlminbGetTheta <- function(nlm, ui) {
  .iniDf <- ui$iniDf
  setNames(vapply(seq_along(.iniDf$name),
                  function(i) {
                    if (.iniDf$fix[i]) {
                      .iniDf$est[i]
                    } else {
                      nlm$par[.iniDf$name[i]]
                    }
                  }, double(1), USE.NAMES=FALSE),
           .iniDf$name)
}
.nlminbControlToFoceiControl <- function(env, assign=TRUE) {
  .nlminbControl <- env$nlminbControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$nlminbControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.nlminbControl$sumProd,
                                optExpression=.nlminbControl$optExpression,
                                literalFix=.nlminbControl$literalFix,
                                literalFixRes=.nlminbControl$literalFixRes,
                                scaleTo=0,
                                calcTables=.nlminbControl$calcTables,
                                addProp=.nlminbControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.nlminbControl$compress,
                                ci=.nlminbControl$ci,
                                sigdigTable=.nlminbControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

.nlminbFamilyFit <- function(env, ...) {
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
  .nlminb <- .collectWarn(.nlminbFitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$nlminb <- .nlminb[[1]]
  .ret$parHistData <- .ret$nlminb$parHistData
  .ret$nlm$parHistData <- NULL
  .ret$message <- .ret$nlminb$message
  if (rxode2::rxGetControl(.ui, "returnNlminb", FALSE)) {
    return(.ret$nlminb)
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .nlminbGetTheta(.ret$nlminb, .ui)
  .ret$cov <- .ret$nlminb$cov
  .ret$covMethod <- .ret$nlminb$covMethod
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
  .ret$est <- "nlminb"
  # There is no parameter history for nlme
  .ret$objective <- 2 * as.numeric(.ret$nlminb$objective)
  .ret$model <- .ui$ebe
  .ret$ofvType <- "nlminb"
  .nlminbControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="nlminb")
  .env <- .ret$env
  .env$method <- "nlminb"
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.nlminb <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'nlminb', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'nlminb'", .var.name=.ui$modelName)
  .nlminbFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .nlminbFamilyFit(env,  ...)
}
attr(nlmixr2Est.nlminb, "covPresent") <- TRUE
