#' Control for newuoa estimation method in nlmixr2
#'
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @inheritParams bobyqaControl
#'
#' @param returnNewuoa return the newuoa output instead of the nlmixr2
#'   fit
#' @return newuoa control structure
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
#' fit2 <- nlmixr(mod, dsn, est="newuoa")
#'
#' print(fit2)
#'
#' # you can also get the nlm output with
#'
#' fit2$newuoa
#'
#' # The nlm control has been modified slightly to include
#' # extra components and name the parameters
#' }
newuoaControl <- function(npt=NULL,
                          rhobeg=NULL,
                          rhoend=NULL,
                          iprint=0L,
                          maxfun=100000L,
                          returnNewuoa=FALSE,
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
                          calcTables=TRUE, compress=FALSE,
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
  checkmate::assertLogical(returnNewuoa, len=1, any.missing=FALSE)
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
               returnNewuoa=returnNewuoa,

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
  class(.ret) <- "newuoaControl"
  .ret
}

#' @export
rxUiDeparse.newuoaControl <- function(object, var) {
  .default <- newuoaControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' Get the newuoa family control
#'
#' @param env newuoa optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.newuoaFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2est::newuoaControl()
  }
  if (!inherits(.control, "newuoaControl")){
    .control <- do.call(nlmixr2est::newuoaControl, .control)
  }
  assign("control", .control, envir=.ui)
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.newuoaControl <- function(control, env) {
  ## eval(rxode2::rxUiDeparse(control, "control"))
  assign("newuoaControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.newuoa <- function(x, ...) {
  .env <- x[[1]]
  if (exists("newuoaControl", .env)) {
    .control <- get("newuoaControl", .env)
    if (inherits(.control, "newuoaControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "newuoaControl")) return(.control)
  }
  stop("cannot find newuoa related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.newuoa <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- newuoaControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("newuoaControl", .ctl)
  if (!inherits(.ctl, "newuoaControl")) {
    .minfo("invalid control for `est=\"newuoa\"`, using default")
    .ctl <- newuoaControl()
  } else {
    .ctl <- do.call(newuoaControl, .ctl)
  }
  .ctl
}

.newuoaControlToFoceiControl <- function(env, assign=TRUE) {
  .newuoaControl <- env$newuoaControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$newuoaControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.newuoaControl$sumProd,
                                optExpression=.newuoaControl$optExpression,
                                literalFix=.newuoaControl$literalFix,
                                literalFixRes=.newuoaControl$literalFixRes,
                                scaleTo=0,
                                calcTables=.newuoaControl$calcTables,
                                addProp=.newuoaControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.newuoaControl$compress,
                                ci=.newuoaControl$ci,
                                sigdigTable=.newuoaControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

.newuoaFitModel <- function(ui, dataSav) {
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
  .ret <- bquote(minqa::newuoa(
    par=.(.env$par.ini),
    fn=.(nlmixr2est::.nlmixrOptimFunC),
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
.newuoaGetTheta <- function(nlm, ui) {
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

.newuoaFamilyFit <- function(env, ...) {
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
  .newuoa <- .collectWarn(.newuoaFitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$newuoa <- .newuoa[[1]]
  .ret$parHistData <- .ret$newuoa$parHistData
  .ret$newuoa$parHistData <- NULL
  .ret$message <- .ret$newuoa$message
  if (rxode2::rxGetControl(.ui, "returnNewuoa", FALSE)) {
    return(.ret$newuoa)
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .newuoaGetTheta(.ret$newuoa, .ui)
  .ret$cov <- .ret$newuoa$cov
  .ret$covMethod <- .ret$newuoa$covMethod
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
  .ret$est <- "newuoa"
  # There is no parameter history for nlme
  .ret$objective <- 2 * as.numeric(.ret$newuoa$fval)
  .ret$model <- .ui$ebe
  .ret$ofvType <- "newuoa"
  .newuoaControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="newuoa")
  .env <- .ret$env
  .env$method <- "newuoa"
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.newuoa <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'newuoa', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'newuoa'", .var.name=.ui$modelName)
  .newuoaFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .newuoaFamilyFit(env,  ...)
}
attr(nlmixr2Est.newuoa, "covPresent") <- TRUE

#minqa::newuoa()
