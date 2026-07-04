#' Control for lbfgsb3c estimation method in nlmixr2
#'
#' @inheritParams iterPrintParams
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @inheritParams lbfgsb3c::lbfgsb3c
#' @inheritParams nlmControl
#'
#' @param returnLbfgsb3c return the lbfgsb3c output instead of the nlmixr2
#'   fit
#'
#' @param trace If positive, print tracing information; higher values give
#'   more detail (see source for "L-BFGS-B" trace levels).
#'
#' @param factr Convergence tolerance factor for "L-BFGS-B"; converges when
#'   the objective reduction is within this factor of machine tolerance
#'   (default 1e7, i.e. ~1e-8).
#'
#' @param pgtol Tolerance on the projected gradient for "L-BFGS-B"; 0
#'   (default) suppresses the check.
#'
#' @param abstol Absolute x-value tolerance for "L-BFGS-B"; 0 (default)
#'   suppresses the check.
#'
#' @param reltol Relative x-value tolerance for "L-BFGS-B"; 0 (default)
#'   suppresses the check.
#'
#' @param lmm Number of BFGS updates retained in "L-BFGS-B" (default 5).
#'
#' @param maxit maximum number of iterations.

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
#' fit2 <- nlmixr(mod, dsn, est="lbfgsb3c")
#'
#' print(fit2)
#'
#' # you can also get the nlm output with fit2$lbfgsb3c
#'
#' fit2$lbfgsb3c
#'
#' # The nlm control has been modified slightly to include
#' # extra components and name the parameters
#' }
lbfgsb3cControl <- function(trace=0,
                            factr=1e7,
                            pgtol=0,
                            abstol=0,
                            reltol=0,
                            lmm=5L,
                            maxit=10000L,
                            returnLbfgsb3c=FALSE,
                            stickyRecalcN=4,
                            maxOdeRecalc=5,
                            odeRecalcFactor=10^(0.5),
                            indTolRelax=TRUE,

                            useColor = NULL,
                            printNcol = NULL, #
                            print = 1L, #

                            normType = c("rescale2", "mean", "rescale", "std", "len", "constant"), #
                            scaleType = c("nlmixr2", "norm", "mult", "multAdd"), #
                            scaleCmax = 1e5, #
                            scaleCmin = 1e-5, #
                            scaleC=NULL,
                            scaleTo=1.0,
                            gradTo=1.0,

                            rxControl=NULL,
                            optExpression=TRUE, sumProd=FALSE,
                            literalFix=TRUE,
                            literalFixRes=TRUE,
                            addProp = c("combined2", "combined1"),
                            eventSens = c("jump", "fd"),
                            sensMethod = c("forward", "adjoint", "auto"),
                            calcTables=TRUE, compress=FALSE,
                            covMethod=c("r", ""),
                            adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {

  checkmate::assertIntegerish(trace, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(factr, len=1, any.missing=FALSE, lower=10)
  checkmate::assertNumeric(pgtol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(abstol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(reltol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertIntegerish(lmm, len=1, any.missing=FALSE, lower=1)

  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFixRes, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnLbfgsb3c, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl", "iterPrintControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)
  checkmate::assertLogical(indTolRelax, any.missing=FALSE, len=1)

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

  .iterPrintControl <- .absorbIterPrintControl(print = print,
                                               printNcol = printNcol,
                                               useColor = useColor,
                                               iterPrintControl = .xtra$iterPrintControl)
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

  .ret <- list(
    trace=trace,
    factr=factr,
    pgtol=pgtol,
    abstol=abstol,
    reltol=reltol,
    lmm=lmm,

    covMethod=match.arg(covMethod),
    optExpression=optExpression,
    literalFix=literalFix,
    literalFixRes=literalFixRes,
    sumProd=sumProd,
    rxControl=rxControl,
    returnLbfgsb3c=returnLbfgsb3c,
    gradTo=gradTo,

    stickyRecalcN=as.integer(stickyRecalcN),
    maxOdeRecalc=as.integer(maxOdeRecalc),
    odeRecalcFactor=odeRecalcFactor,
    indTolRelax=indTolRelax,

    iterPrintControl = .iterPrintControl,
    scaleType=scaleType,
    normType=normType,

    scaleCmax=scaleCmax,
    scaleCmin=scaleCmin,
    scaleC=scaleC,
    scaleTo=scaleTo,

    addProp=match.arg(addProp),
    eventSens=match.arg(eventSens),
    sensMethod=match.arg(sensMethod),
    calcTables=calcTables,
    compress=compress,
    ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
    genRxControl=.genRxControl)
  class(.ret) <- "lbfgsb3cControl"
  .ret
}

#' @export
rxUiDeparse.lbfgsb3cControl <- function(object, var) {
  .default <- lbfgsb3cControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' Get the lbfgsb3c family control
#'
#' @param env lbfgsb3c optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.lbfgsb3cFamilyControl <- function(env, ...) {
  .nlmFamilyControlGeneric(env, nlmixr2est::lbfgsb3cControl, "lbfgsb3cControl")
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.lbfgsb3cControl <- function(control, env) {
  assign("lbfgsb3cControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.lbfgsb3c <- function(x, ...) {
  .env <- x[[1]]
  if (exists("lbfgsb3cControl", .env)) {
    .control <- get("lbfgsb3cControl", .env)
    if (inherits(.control, "lbfgsb3cControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "lbfgsb3cControl")) return(.control)
  }
  stop("cannot find lbfgsb3c related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.lbfgsb3c <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- lbfgsb3cControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("lbfgsb3cControl", .ctl)
  if (!inherits(.ctl, "lbfgsb3cControl")) {
    .minfo("invalid control for `est=\"lbfgsb3c\"`, using default")
    .ctl <- lbfgsb3cControl()
  } else {
    .ctl <- do.call(lbfgsb3cControl, .ctl)
  }
  .ctl
}

.lbfgsb3cControlToFoceiControl <- function(env, assign=TRUE) {
  .lbfgsb3cControl <- env$lbfgsb3cControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$lbfgsb3cControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.lbfgsb3cControl$sumProd,
                                optExpression=.lbfgsb3cControl$optExpression,
                                literalFix=.lbfgsb3cControl$literalFix,
                                literalFixRes=.lbfgsb3cControl$literalFixRes,
                                scaleTo=0,
                                calcTables=.lbfgsb3cControl$calcTables,
                                addProp=.lbfgsb3cControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.lbfgsb3cControl$compress,
                                ci=.lbfgsb3cControl$ci,
                                sigdigTable=.lbfgsb3cControl$sigdigTable,
                                indTolRelax=.lbfgsb3cControl$indTolRelax,
                                eventSens=.lbfgsb3cControl$eventSens,
                                sensMethod=.lbfgsb3cControl$sensMethod)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

.lbfgsb3cFitModel <- function(ui, dataSav) {
  # Use nlmEnv and function for DRY principle
  rxode2::rxReq("lbfgsb3c")
  .ctl <- ui$control
  .keep <- c("trace","factr","pgtol", "abstol", "reltol", "lmm")
  .keep <- .keep[vapply(.keep, function(opt) {
    !is.null(.ctl[[opt]])
  }, logical(1), USE.NAMES = FALSE)]
  .oCtl <- setNames(lapply(.keep, function(x) {.ctl[[x]]}), .keep)
  class(.ctl) <- NULL

  .p <- setNames(ui$nlmParIni, ui$nlmParName)
  .mi <- ui$nlmSensModel
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl,
                       lower=ui$optimParLower, upper=ui$optimParUpper)
  on.exit({.nlmFreeEnv()})
  # support gradient
  .ret <- bquote(lbfgsb3c::lbfgsb3c(
    par=.(.env$par.ini),
    # fn is called every eval too, so use .nlmixrOptimFunC like gr does
    fn=.(nlmixr2est::.nlmixrOptimFunC),
    gr=.(nlmixr2est::.nlmixrOptimGradC),
    control=.(.oCtl),
    lower=.(.env$lower),
    upper=.(.env$upper)))
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
.lbfgsb3cGetTheta <- function(nlm, ui) {
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

.lbfgsb3cFamilyFit <- function(env, ...) {
  .nlmFamilyFitGeneric(
    env, "lbfgsb3c", .lbfgsb3cFitModel, .lbfgsb3cGetTheta,
    objective = function(.fit) 2 * as.numeric(.fit$value),
    controlToFocei = .lbfgsb3cControlToFoceiControl,
    returnFlag = "returnLbfgsb3c")
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.lbfgsb3c <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'lbfgsb3c', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'lbfgsb3c'", .var.name=.ui$modelName)
  .lbfgsb3cFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .lbfgsb3cFamilyFit(env,  ...)
}
attr(nlmixr2Est.lbfgsb3c, "covPresent") <- TRUE
attr(nlmixr2Est.lbfgsb3c, "unbounded") <- FALSE
