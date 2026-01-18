#' Control for bobyqa estimation method in nlmixr2
#'
#' @inheritParams foceiControl
#' @inheritParams saemControl
#'
#' @param returnBobyqa return the bobyqa output instead of the nlmixr2
#'   fit
#'
#' @param npt The number of points used to approximate the objective
#'   function via a quadratic approximation. The value of npt must be
#'   in the interval [n+2,(n+1)(n+2)/2] where n is the number of
#'   parameters in `par`. Choices that exceed 2*n+1 are not
#'   recommended.  If not defined, it will be set to min(n * 2, n+2).
#'
#' @param rhobeg `rhobeg` and `rhoend` must be set to the initial and
#'   final values of a trust region radius, so both must be positive
#'   with `0 < rhoend < rhobeg`. Typically `rhobeg` should be about
#'   one tenth of the greatest expected change to a variable.  If the
#'   user does not provide a value, this will be set to `min(0.95, 0.2
#'   * max(abs(par)))`.  Note also that smallest difference
#'   `abs(upper-lower)` should be greater than or equal to `rhobeg*2`.
#'   If this is not the case then `rhobeg` will be adjusted.
#' @param rhoend The smallest value of the trust region radius that is
#'   allowed. If not defined, then 1e-6 times the value set for
#'   `rhobeg` will be used.
#' @param iprint The value of `iprint` should be set to an integer
#'   value in `0, 1, 2, 3, ...`, which controls the amount of
#'   printing.  Specifically, there is no output if `iprint=0` and
#'   there is output only at the start and the return if `iprint=1`.
#'   Otherwise, each new value of `rho` is printed, with the best
#'   vector of variables so far and the corresponding value of the
#'   objective function. Further, each new value of the objective
#'   function with its variables are output if `iprint=3`.  If `iprint
#'   > 3`, the objective function value and corresponding variables
#'   are output every `iprint` evaluations.  Default value is `0`.
#' @param maxfun The maximum allowed number of function
#'   evaluations. If this is exceeded, the method will terminate.
#' @return bobqya control structure
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
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
#' fit2 <- nlmixr(mod, dsn, est="bobyqa")
#'
#' print(fit2)
#'
#' # you can also get the nlm output with
#'
#' fit2$bobyqa
#'
#' # The nlm control has been modified slightly to include
#' # extra components and name the parameters
#' }
bobyqaControl <- function(npt=NULL,
                          rhobeg=NULL,
                          rhoend=NULL,
                          iprint=0L,
                          maxfun=100000L,
                          returnBobyqa=FALSE,
                          stickyRecalcN=4,
                          maxOdeRecalc=5,
                          odeRecalcFactor=10^(0.5),

                          useColor = crayon::has_color(),
                          printNcol = floor((getOption("width") - 23) / 12), #
                          print = 1L, #

                          normType = c("rescale2", "mean", "rescale", "std", "len", "constant"), #
                          scaleType = c("nlmixr2", "norm", "mult", "multAdd"), #
                          scaleCmax = 1e5, #
                          scaleCmin = 1e-5, #
                          scaleC=NULL,
                          scaleTo=1.0,

                          rxControl=NULL,
                          optExpression=TRUE, sumProd=FALSE,
                          literalFix=TRUE,
                          literalFixRes=TRUE,
                          addProp = c("combined2", "combined1"),
                          calcTables=TRUE, compress=TRUE,
                          covMethod=c("r", ""),
                          adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {

  checkmate::assertIntegerish(npt, null.ok=TRUE, any.missing=FALSE, lower=2, len=1)
  checkmate::assertNumeric(rhobeg, null.ok=TRUE, any.missing=FALSE, lower=0, len=1)
  checkmate::assertNumeric(rhoend, null.ok=TRUE, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(iprint, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxfun, any.missing=FALSE, lower=10, len=1)

  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFixRes, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnBobyqa, len=1, any.missing=FALSE)
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

  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)

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
  checkmate::assertIntegerish(printNcol, len=1, lower=0, any.missing=FALSE)
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

  .ret <- list(npt=npt,
               rhobeg=rhobeg,
               rhoend=rhoend,
               iprint=iprint,
               maxfun=maxfun,
               covMethod=match.arg(covMethod),
               optExpression=optExpression,
               literalFix=literalFix,
               literalFixRes=literalFixRes,
               sumProd=sumProd,
               rxControl=rxControl,
               returnBobyqa=returnBobyqa,

               stickyRecalcN=as.integer(stickyRecalcN),
               maxOdeRecalc=as.integer(maxOdeRecalc),
               odeRecalcFactor=odeRecalcFactor,

               useColor=useColor,
               print=print,
               printNcol=printNcol,
               scaleType=scaleType,
               normType=normType,

               scaleCmax=scaleCmax,
               scaleCmin=scaleCmin,
               scaleC=scaleC,
               scaleTo=scaleTo,

               addProp=match.arg(addProp),
               calcTables=calcTables,
               compress=compress,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "bobyqaControl"
  .ret
}

#' @export
rxUiDeparse.bobyqaControl <- function(object, var) {
  .default <- bobyqaControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' Get the bobyqa family control
#'
#' @param env bobyqa optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.bobyqaFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2est::bobyqaControl()
  }
  if (!inherits(.control, "bobyqaControl")){
    .control <- do.call(nlmixr2est::bobyqaControl, .control)
  }
  assign("control", .control, envir=.ui)
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.bobyqaControl <- function(control, env) {
  eval(rxode2::rxUiDeparse(control, "control"))
  assign("bobyqaControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.bobyqa <- function(x, ...) {
  .env <- x[[1]]
  if (exists("bobyqaControl", .env)) {
    .control <- get("bobyqaControl", .env)
    if (inherits(.control, "bobyqaControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "bobyqaControl")) return(.control)
  }
  stop("cannot find bobyqa related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.bobyqa <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- bobyqaControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("bobyqaControl", .ctl)
  if (!inherits(.ctl, "bobyqaControl")) {
    .minfo("invalid control for `est=\"bobyqa\"`, using default")
    .ctl <- bobyqaControl()
  } else {
    .ctl <- do.call(bobyqaControl, .ctl)
  }
  .ctl
}

.bobyqaControlToFoceiControl <- function(env, assign=TRUE) {
  .bobyqaControl <- env$bobyqaControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$bobyqaControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.bobyqaControl$sumProd,
                                optExpression=.bobyqaControl$optExpression,
                                literalFix=.bobyqaControl$literalFix,
                                literalFixRes=.bobyqaControl$literalFixRes,
                                scaleTo=0,
                                calcTables=.bobyqaControl$calcTables,
                                addProp=.bobyqaControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.bobyqaControl$compress,
                                ci=.bobyqaControl$ci,
                                sigdigTable=.bobyqaControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

.bobyqaFitModel <- function(ui, dataSav) {
  # Use nlmEnv and function for DRY principle
  rxode2::rxReq("minqa")
  .ctl <- ui$control
  .keep <- c("npt", "rhobeg", "rhoend", "iprint", "maxfun")
  .keep <- .keep[vapply(.keep, function(opt) {
    !is.null(.ctl[[opt]])
  }, logical(1), USE.NAMES = FALSE)]

  .oCtl <- setNames(lapply(.keep, function(x) {.ctl[[x]]}), .keep)
  class(.ctl) <- NULL
  .p <- setNames(ui$nlmParIni, ui$nlmParName)
  .mi <-  ui$nlmRxModel
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl,
                       lower=ui$optimParLower, upper=ui$optimParUpper)
  on.exit({.nlmFreeEnv()})
  # support gradient
  .ret <- bquote(minqa::bobyqa(
    par=.(.env$par.ini),
    fn=.(nlmixr2est::.nlmixrOptimFunC),
    lower=.(.env$lower),
    upper=.(.env$upper),
    control=.(.oCtl)))
  .ret <- eval(.ret)
  .nlmFinalizeList(.env, .ret, par="par", printLine=TRUE,
                   hessianCov=TRUE)
}

#' Get the full theta for nlm methods
#'
#' @param optim enhanced nlm return
#' @param ui ui object
#' @return named theta matrix
#' @author Matthew L. Fidler
#' @noRd
.bobyqaGetTheta <- function(nlm, ui) {
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

.bobyqaFamilyFit <- function(env, ...) {
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
  .bobyqa <- .collectWarn(.bobyqaFitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$bobyqa <- .bobyqa[[1]]
  .ret$parHistData <- .ret$bobyqa$parHistData
  .ret$bobyqa$parHistData <- NULL
  .ret$message <- .ret$bobyqa$message
  if (rxode2::rxGetControl(.ui, "returnBobyqa", FALSE)) {
    return(.ret$bobyqa)
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .bobyqaGetTheta(.ret$bobyqa, .ui)
  .ret$cov <- .ret$bobyqa$cov
  .ret$covMethod <- .ret$bobyqa$covMethod
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
  .ret$est <- "bobyqa"
  # There is no parameter history for nlme
  .ret$objective <- 2 * as.numeric(.ret$bobyqa$fval)
  .ret$model <- .ui$ebe
  .ret$ofvType <- "bobyqa"
  .bobyqaControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="bobyqa")
  .env <- .ret$env
  .env$method <- "bobyqa"
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.bobyqa <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'bobyqa', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'bobyqa'", .var.name=.ui$modelName)
  .bobyqaFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .bobyqaFamilyFit(env,  ...)
}
attr(nlmixr2Est.bobyqa, "covPresent") <- TRUE

#minqa::bobyqa()
