#' Control for n1qn1 estimation method in nlmixr2
#'
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @inheritParams n1qn1::n1qn1
#' @inheritParams nlmControl
#'
#' @param returnN1qn1 return the n1qn1 output instead of the nlmixr2
#'   fit
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
#' fit2 <- nlmixr(mod, dsn, est="n1qn1")
#'
#' print(fit2)
#'
#' # you can also get the nlm output with fit2$n1qn1
#'
#' fit2$n1qn1
#'
#' # The nlm control has been modified slightly to include
#' # extra components and name the parameters
#' }
n1qn1Control <- function(epsilon = (.Machine$double.eps) ^ 0.25,
                         max_iterations = 10000,
                         nsim = 10000,
                         imp = 0,
                         print.functions=FALSE,

                         returnN1qn1=FALSE,
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
                         gradTo=1.0,

                         rxControl=NULL,
                         optExpression=TRUE, sumProd=FALSE,
                         literalFix=TRUE,
                         literalFixRes=TRUE,
                         addProp = c("combined2", "combined1"),
                         calcTables=TRUE, compress=FALSE,
                         covMethod=c("r", "n1qn1", ""),
                         adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {

  checkmate::assertNumeric(epsilon, len=1, any.missing=FALSE, lower=0)
  checkmate::assertIntegerish(max_iterations, len=1, any.missing=FALSE, lower=10)
  checkmate::assertIntegerish(nsim, len=1, any.missing=FALSE, lower=10)
  checkmate::assertIntegerish(imp, len=1, any.missing=FALSE, lower=0)
  checkmate::assertLogical(print.functions, len=1, any.missing=FALSE)

  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFixRes, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnN1qn1, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% "genRxControl")]
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
  checkmate::assertNumeric(gradTo, len=1, lower=0, any.missing=FALSE)

  .ret <- list(
    epsilon = epsilon,
    max_iterations = max_iterations,
    nsim = nsim,
    imp = imp,
    print.functions=print.functions,
    covMethod=match.arg(covMethod),
    optExpression=optExpression,
    literalFix=literalFix,
    literalFixRes=literalFixRes,
    sumProd=sumProd,
    rxControl=rxControl,
    returnN1qn1=returnN1qn1,
    gradTo=gradTo,

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
  class(.ret) <- "n1qn1Control"
  .ret
}

#' @export
rxUiDeparse.n1qn1Control <- function(object, var) {
  .default <- n1qn1Control()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}


#' Get the n1qn1 family control
#'
#' @param env n1qn1 optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.n1qn1FamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2est::n1qn1Control()
  }
  if (!inherits(.control, "n1qn1Control")) {
    .control <- do.call(nlmixr2est::n1qn1Control, .control)
  }
  assign("control", .control, envir=.ui)
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.n1qn1Control <- function(control, env) {
  assign("n1qn1Control", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.n1qn1 <- function(x, ...) {
  .env <- x[[1]]
  if (exists("n1qn1Control", .env)) {
    .control <- get("n1qn1Control", .env)
    if (inherits(.control, "n1qn1Control")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "n1qn1Control")) return(.control)
  }
  stop("cannot find n1qn1 related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.n1qn1 <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- n1qn1Control()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("n1qn1Control", .ctl)
  if (!inherits(.ctl, "n1qn1Control")) {
    .minfo("invalid control for `est=\"n1qn1\"`, using default")
    .ctl <- n1qn1Control()
  } else {
    .ctl <- do.call(n1qn1Control, .ctl)
  }
  .ctl
}

.n1qn1ControlToFoceiControl <- function(env, assign=TRUE) {
  .n1qn1Control <- env$n1qn1Control
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$n1qn1Control$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.n1qn1Control$sumProd,
                                optExpression=.n1qn1Control$optExpression,
                                literalFix=.n1qn1Control$literalFix,
                                literalFixRes=.n1qn1Control$literalFixRes,
                                scaleTo=0,
                                calcTables=.n1qn1Control$calcTables,
                                addProp=.n1qn1Control$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.n1qn1Control$compress,
                                ci=.n1qn1Control$ci,
                                sigdigTable=.n1qn1Control$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

.n1qn1FitModel <- function(ui, dataSav) {
  # Use nlmEnv and function for DRY principle
  rxode2::rxReq("n1qn1")
  .ctl <- ui$control
  class(.ctl) <- NULL
  .p <- setNames(ui$nlmParIni, ui$nlmParName)
  .mi <- ui$nlmSensModel
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl,
                       lower=ui$optimParLower, upper=ui$optimParUpper)
  on.exit({.nlmFreeEnv()})
  # support gradient
  .ret <- bquote(n1qn1::n1qn1(
    # Calls grad with every function evaluation, use .nlmixrOptimFunC
    # which does as well
    call_eval=.(nlmixr2est::.nlmixrOptimFunC),
    #call_eval=.(nlmixr2est::.nlmixrNlminbFunC),
    call_grad=.(nlmixr2est::.nlmixrOptimGradC),
    #call_grad=.(nlmixr2est::.nlmixrNlminbGradC),
    vars=.(.env$par.ini),
    epsilon = .(.ctl$epsilon),
    max_iterations = .(.ctl$max_iterations),
    nsim = .(.ctl$nsim),
    imp = .(.ctl$imp)))
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
.n1qn1GetTheta <- function(nlm, ui) {
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

.n1qn1FamilyFit <- function(env, ...) {
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
  .n1qn1 <- .collectWarn(.n1qn1FitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$n1qn1 <- .n1qn1[[1]]
  .ret$parHistData <- .ret$n1qn1$parHistData
  .ret$n1qn1$parHistData <- NULL
  .ret$message <- .ret$n1qn1$message
  if (rxode2::rxGetControl(.ui, "returnN1qn1", FALSE)) {
    return(.ret$n1qn1)
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .n1qn1GetTheta(.ret$n1qn1, .ui)
  .ret$cov <- .ret$n1qn1$cov
  .ret$covMethod <- .ret$n1qn1$covMethod
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
  .ret$est <- "n1qn1"
  # There is no parameter history for nlme
  .ret$objective <- 2 * as.numeric(.ret$n1qn1$value)
  .ret$model <- .ui$ebe
  .ret$ofvType <- "n1qn1"
  .n1qn1ControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="n1qn1")
  .env <- .ret$env
  .env$method <- "n1qn1"
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.n1qn1 <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'n1qn1', try 'focei'",
                                   .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'n1qn1'",
                                   .var.name=.ui$modelName)
  rxode2::warnRxBounded(.ui, " which are ignored in 'n1qn1'",
                        .var.name=.ui$modelName)
  .n1qn1FamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .n1qn1FamilyFit(env,  ...)
}
attr(nlmixr2Est.n1qn1, "covPresent") <- TRUE
