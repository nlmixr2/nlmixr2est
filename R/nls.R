#' nlmixr2 defaults controls for nls
#'
#' @inheritParams stats::nls
#' @inheritParams stats::nls.control
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @inheritParams nlmControl
#' @inheritParams minpack.lm::nls.lm.control
#' @param returnNls logical; when TRUE, will return the nls object
#'   instead of the nlmixr object
#' @return nls control object
#' @export
#' @author Matthew L. Fidler
#' @examples
#' \donttest{
#'
#' one.cmt <- function() {
#'   ini({
#'    tka <- 0.45
#'    tcl <- log(c(0, 2.7, 100))
#'    tv <- 3.45
#'    add.sd <- 0.7
#'  })
#'  model({
#'    ka <- exp(tka)
#'    cl <- exp(tcl)
#'    v <- exp(tv)
#'    linCmt() ~ add(add.sd)
#'  })
#' }
#'
#' # Uses nlsLM from minpack.lm if available
#'
#' fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nls", nlsControl(algorithm="LM"))
#'
#' # Uses port and respect parameter boundaries
#' fit2 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nls", nlsControl(algorithm="port"))
#'
#' # You can access the underlying nls object with `$nls`
#' fit2$nls
#'
#' }
nlsControl <- function(maxiter=10000,
                       tol = 1e-05,
                       minFactor = 1/1024,
                       printEval = FALSE,
                       warnOnly = FALSE,
                       scaleOffset = 0,
                       nDcentral = FALSE,
                       algorithm = c("LM", "default", "plinear", "port"),
                       ############################################
                       ## minpack.lm
                       ftol = sqrt(.Machine$double.eps),
                       ptol = sqrt(.Machine$double.eps),
                       gtol = 0,
                       diag = list(),
                       epsfcn = 0,
                       factor = 100,
                       maxfev = integer(),
                       nprint = 0,

                       #### nlm C++ style style to give gradients
                       solveType=c("grad", "fun"),

                       stickyRecalcN=4,
                       maxOdeRecalc=5,
                       odeRecalcFactor=10^(0.5),

                       eventType=c("central", "forward"),
                       shiErr=(.Machine$double.eps)^(1/3),
                       shi21maxFD=20L,

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

                       ############################################
                       trace = FALSE, #nolint
                       rxControl=NULL,
                       optExpression=TRUE, sumProd=FALSE,
                       literalFix=TRUE,
                       returnNls=FALSE,
                       addProp = c("combined2", "combined1"),
                       calcTables=TRUE, compress=TRUE,
                       adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {
  algorithm <- match.arg(algorithm)
  if (algorithm == "LM" && !requireNamespace("minpack.lm", quietly = TRUE)) {
    .malert("to use the LM algorithm you must have minpack.lm installed")
    .malert("changing to default `nls` method")
    algorithm <- "default"
  }
  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)
  checkmate::assertNumeric(shiErr, lower=0, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(shi21maxFD, lower=1, any.missing=FALSE, len=1)

  .eventTypeIdx <- c("central" =2L, "forward"=1L)
  if (checkmate::testIntegerish(eventType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    eventType <- as.integer(eventType)
  } else {
    eventType <- setNames(.eventTypeIdx[match.arg(eventType)], NULL)
  }

  solveType <- match.arg(solveType)

  checkmate::assertNumeric(ftol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(ptol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(gtol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(epsfcn, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(factor, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(maxfev, min.len=0, max.len=1, any.missing=FALSE)

  checkmate::assertLogical(trace, len=1, any.missing=FALSE) # nolint
  checkmate::assertLogical(nDcentral, len=1, any.missing=FALSE)
  checkmate::assertNumeric(scaleOffset, any.missing=FALSE, finite=TRUE)
  checkmate::assertLogical(warnOnly, len=1, any.missing = FALSE)
  checkmate::assertLogical(printEval, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(maxiter, len=1, any.missing=FALSE, lower=1)
  checkmate::assertNumeric(tol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(minFactor, len=1, any.missing=FALSE, lower=0)
  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnNls, len=1, any.missing=FALSE)
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


  .ret <- list(algorithm=algorithm, maxiter=maxiter,
               tol=tol,
               trace=trace, #nolint
               minFactor=minFactor,
               printEval=printEval,
               warnOnly=warnOnly,
               scaleOffset=scaleOffset,
               nDcentral=nDcentral,
               solveType=solveType,
               stickyRecalcN=stickyRecalcN,
               maxOdeRecalc=maxOdeRecalc,
               odeRecalcFactor=odeRecalcFactor,
               eventType=eventType,
               shiErr=shiErr,
               shi21maxFD=shi21maxFD,
               ftol = ftol,
               ptol = ptol,
               gtol = gtol,
               diag = diag,
               epsfcn = epsfcn,
               factor = factor,
               maxfev = maxfev,
               nprint = nprint,
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

               optExpression=optExpression,
               literalFix=literalFix,
               sumProd=sumProd,
               rxControl=rxControl,
               returnNls=returnNls,
               addProp=match.arg(addProp),
               calcTables=calcTables,
               compress=compress,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "nlsControl"
  .ret
}

#' @export
rxUiDeparse.nlsControl <- function(object, var) {
  .default <- nlsControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' Get the nls family control
#'
#' @param env nlme optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlsFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2est::nlsControl()
  }
  if (!inherits(.control, "nlsControl")){
    .control <- do.call(nlmixr2est::nlsControl, .control)
  }
  assign("control", .control, envir=.ui)
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.nlsControl <- function(control, env) {
  eval(rxode2::rxUiDeparse(control, "control"))
  assign("nlsControl", control, envir=env)
}


#' @rdname nmObjGetControl
#' @export
nmObjGetControl.nls <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nlsControl", .env)) {
    .control <- get("nlsControl", .env)
    if (inherits(.control, "nlsControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nlsControl")) return(.control)
  }
  stop("cannot find nls related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.nls <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- nlsControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("nlsControl", .ctl)
  if (!inherits(.ctl, "nlsControl")) {
    .minfo("invalid control for `est=\"nls\"`, using default")
    .ctl <- nlsControl()
  } else {
    .ctl <- do.call(nlsControl, .ctl)
  }
  .ctl
}


#' A surrogate function for nls to call for ode solving
#'
#' @param DV dependent variable
#' @param ... Other parameters fed to prediction function
#' @return Predictions
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.nlmixrNlsFun <- function(DV, ...) {
  do.call(rxode2::rxSolve,
          c(list(object=nlmixr2global$nlsEnv$model,
                 params=nlmixr2global$nlsEnv$parFun(...),
                 events=nlmixr2global$nlsEnv$data),
            nlmixr2global$nlsEnv$rxControl))$rx_pred_
}
#' @rdname dot-nlmixrNlsFun
#' @export
.nlmixrNlsFunValGrad <- function(DV, ...) {
  .Call(`_nlmixr2est_solveGradNls`, c(...), 1L)
}
#' Internal nls functions for minpack.lm
#' @param x Parameter for estimate
#' @keywords internal
#' @export
.nlmixrNlsFunVal <- function(x) {
  .Call(`_nlmixr2est_solveGradNls`, x, 2L)
}

#' @rdname dot-nlmixrNlsFunVal
#' @export
.nlmixrNlsFunGrad <- function(x) {
  .Call(`_nlmixr2est_solveGradNls`, x, 3L)
}

#' Returns the data currently setup to run nls
#'
#' @return Returns the data currently setup to run nls
#' @export
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
.nlmixrNlsData <- function() {
  nlmixr2global$nlsEnv$dataNls
}


#' This is a S3 method for getting the distribution lines for a base rxode2 nls problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the focei. This is based
#'   on the idea that the focei parameters are defined
#' @author Matthew Fidler
#' @keywords internal
#' @export
rxGetDistributionNlsLines <- function(line) {
  UseMethod("rxGetDistributionNlsLines")
}

#' @rdname rxGetDistributionNlsLines
#' @export
rxGetDistributionNlsLines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  .errNum <- line[[3]]
  .line <- rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1, .errNum,
                                                        rxPredLlik=.getRxPredLlikOption())
  .yj <- as.double(pred1$transform) - 1
  if (.yj == 2) {
    .lineExtra <- quote(rx_dv_ ~ DV)
  } else if (.yj == 3) {
    .lineExtra <- quote(rx_dv_ ~ log(DV))
  } else {
    .lineExtra <- quote(rx_dv_ ~ rxTBS(DV, rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  }
  .lineExtra <- list(.lineExtra)
  if (pred1$dvid == 1) {
    # First estimated residual error is divided out, since it will be
    # estimated as the residual error by nls
    # add+prop and add+pow are not supported
    .errType <- as.character(pred1$errType)
    if (.errType == "add") {
      # In these cases you are simply dividing out the additive error
      # Simply force this to be one.
      .lineExtra <- c(.lineExtra, list(quote(rx_r_ ~ 1)))
    } else if (.errType == "prop") {
      #   rx_r_ ~ (rx_pred_f_ * prop.sd)^2
      .f <- pred1$f
      .type <- as.character(pred1$errTypeF)
      .lineExtra <- c(.lineExtra,
                      list(switch(.type, untransformed = quote(rx_r_ ~ (rx_pred_f_)^2),
                                  transformed = quote(rx_r_ ~ (rx_pred_)^2),
                                  f = bquote(rx_r_ ~ (.(str2lang(.f)))^2),
                                  none = quote(rx_r_ ~ (rx_pred_f_) ^2))))
    } else if (.errType == "pow") {
      .cnd <- pred1$cond
      if (!is.na(pred1$c)) {
        .p2 <- str2lang(pred1$c)
      } else {
        .w <- which(env$iniDf$err %in% c("pow2", "powF2", "powT2") & env$iniDf$condition == .cnd)
        if (length(.w) == 1L) {
          .p2 <- str2lang(env$iniDf$name[.w])
        } else {
          stop("cannot find exponent of power expression", call.=FALSE)
        }
      }
      .f <- pred1$f
      .type <- as.character(pred1$errTypeF)
      .lineExtra <- c(.lineExtra,
                      list(switch(.type, untransformed = bquote(rx_r_ ~ (rx_pred_f_)^(2 * .(.p2))),
                                  transformed = bquote(rx_r_ ~ (rx_pred_)^(2 * .(.p2))),
                                  f = bquote(rx_r_ ~ (.(str2lang(.f)))^(2 * .(.p2))),
                                  none = quote(rx_r_ ~ (rx_pred_f_) ^(2 * .(.p2))))))
    }
  }
  c(.line, .lineExtra)
}

#' @rdname rxGetDistributionNlsLines
#' @export
rxGetDistributionNlsLines.default <- function(line) {
  stop("only normally related endoints can be used with 'nls', try 'nlm'", call.=FALSE)
}

#' @export
rxGetDistributionNlsLines.rxUi <- function(line) {
  .predDf <- rxUiGet.predDfFocei(list(line, TRUE))
  lapply(seq_along(.predDf$cond), function(c) {
    .mod <- .createFoceiLineObject(line, c)
    rxGetDistributionNlsLines(.mod)
  })
}

#' Get the THETA lines from rxode2 UI for nlscontrol
#'
#' Will assign fixed values and remove error terms
#'
#' @param rxui This is the rxode2 ui object
#' @return The theta/eta lines
#' @author Matthew L. Fidler
#' @noRd
.uiGetNlsTheta <- function(rxui) {
  .iniDf <- rxui$iniDf
  #.w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  .env <- new.env(parent=emptyenv())
  .env$i <- 0
  .w <- which(!(rxui$iniDf$err %in% c("add", "prop", "pow")))
  lapply(.w, function(i) {
    if (rxui$iniDf$fix[i]) {
      return(eval(parse(text=paste0("quote(", .iniDf$name[i], " <- ", rxui$iniDf$est[i],")"))))
    }
    if (rxui$iniDf$err[i] %in% c("add", "prop", "pow")) {
      return(NULL)
    }
    .env$i <- .env$i + 1
    eval(parse(text=paste0("quote(", .iniDf$name[i], " <- THETA[", .env$i,"])")))
  })
}

#' @export
rxUiGet.nlsModel0 <- function(x, ...) {
  .f <- x[[1]]
  .ret <- rxode2::rxCombineErrorLines(.f, errLines=rxGetDistributionNlsLines(.f),
                                      prefixLines=.uiGetNlsTheta(.f),
                                      paramsLine=NA, #.uiGetThetaEtaParams(.f),
                                      modelVars=TRUE,
                                      cmtLines=FALSE,
                                      dvidLine=FALSE)
  ## pred <- (Vm * conc)/(K + conc)
  ## (resp - pred) / sqrt(pred)
  .ret <- .ret[[-1]]
  .w <- seq_along(.ret)
  .w <- .w[-1]
  as.call(c(list(quote(`rxModelVars`)),
            as.call(c(list(quote(`{`)),
                      lapply(.w, function(i){.ret[[i]]}),
                      list(quote(rx_pred_ <- (rx_dv_ - rx_pred_) / sqrt(rx_r_)))))))
}
attr(rxUiGet.nlsModel0, "rstudio") <- quote(rxModelVars({}))


#' Load the nls model into symengine
#'
#' @param x rxode2 UI object
#' @return String for loading into symengine
#' @author Matthew L. Fidler
#' @noRd
.nlsPrune <- function(x) {
  .x <- x[[1]]
  .x <- .x$nlsModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .malert("pruning branches ({.code if}/{.code else}) of nls model...")
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

#' @export
rxUiGet.loadPruneNls <- function(x, ...) {
  .loadSymengine(.nlsPrune(x), promoteLinSens = FALSE)
}
attr(rxUiGet.loadPruneNls, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsRxModel <- function(x, ...) {
  .s <- rxUiGet.loadPruneNls(x, ...)
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  ## .var <- get("rx_r_", envir = .s)
  ## .var <- paste0("rx_r_=", rxode2::rxFromSE(.var))

  ## .dv <- get("rx_dv_", envir = .s)
  ## .dv <- paste0("rx_dv_=", rxode2::rxFromSE(.dv))
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
  if (exists("..maxTheta", .s)) {
    .eventTheta <- rep(0L, .s$..maxTheta)
  } else {
    .eventTheta <- integer(0)
  }
  for (.v in .s$..eventVars) {
    .vars <- as.character(get(.v, envir = .s))
    .vars <- rxode2::rxGetModel(paste0("rx_lhs=", rxode2::rxFromSE(.vars)))$params
    for (.v2 in .vars) {
      .reg <- rex::rex(start, "THETA[", capture(any_numbers), "]", end)
      if (regexpr(.reg, .v2) != -1) {
        .num <- as.numeric(sub(.reg, "\\1", .v2))
        .eventTheta[.num] <- 1L
      }
    }
  }
  .s$.eventTheta <- .eventTheta

  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  if (.sumProd) {
    .malert("stabilizing round off errors in nls model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "nls model")
    .msuccess("done")
  }

  .cmt <-  rxUiGet.foceiCmtPreModel(x, ...)
  .interp <- rxUiGet.interpLinesStr(x, ...)
  if (.interp != "") {
    .cmt <-paste0(.cmt, "\n", .interp)
  }
  list(predOnly =rxode2::rxode2(paste(c(rxUiGet.nlsParams(x, ...), .cmt,
                                        .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")),
       eventTheta=.eventTheta)
}

#' @export
rxUiGet.loadPruneNlsSens <- function(x, ...) {
  .loadSymengine(.nlsPrune(x), promoteLinSens = TRUE)
}
attr(rxUiGet.loadPruneNlsSens, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsThetaS <- function(x, ...) {
  .s <- rxUiGet.loadPruneNlsSens(x, ...)
  .sensEtaOrTheta(.s, theta=TRUE)
}
attr(rxUiGet.nlsThetaS, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsHdTheta <- function(x, ...) {
  .s <- rxUiGet.nlsThetaS(x)
  .stateVars <- rxode2stateOde(.s)
  .predMinusDv <- rxode2::rxGetControl(x[[1]], "predMinusDv", TRUE)
  .grd <- rxode2::rxExpandFEta_(
    .stateVars, .s$..maxTheta,
    ifelse(.predMinusDv, 1L, 2L),
    isTheta=TRUE)
  if (rxode2::.useUtf()) {
    .malert("calculate \u2202(f)/\u2202(\u03B8)")
  } else {
    .malert("calculate d(f)/d(theta)")
  }
  rxode2::rxProgress(dim(.grd)[1])
  on.exit({
    rxode2::rxProgressAbort()
  })
  .any.zero <- FALSE
  .all.zero <- TRUE
  .ret <- apply(.grd, 1, function(x) {
    .l <- x["calc"]
    .l <- eval(parse(text = .l))
    .ret <- paste0(x["dfe"], "=", rxode2::rxFromSE(.l))
    .zErr <- suppressWarnings(try(as.numeric(get(x["dfe"], .s)), silent = TRUE))
    if (identical(.zErr, 0)) {
      .any.zero <<- TRUE
    } else if (.all.zero) {
      .all.zero <<- FALSE
    }
    rxode2::rxTick()
    .ret
  })
  if (.all.zero) {
    stop("none of the predictions depend on 'THETA'", call. = FALSE)
  }
  if (.any.zero) {
    warning("some of the predictions do not depend on 'THETA'", call. = FALSE)
  }
  .s$..HdTheta <- .ret
  .s$..pred.minus.dv <- .predMinusDv
  rxode2::rxProgressStop()
  .s
}
attr(rxUiGet.nlsHdTheta, "rstudio") <- emptyenv()

#' Finalize nls rxode2 based on symengine saved info
#'
#' @param .s Symengine/rxode2 object
#' @return Nothing
#' @author Matthew L Fidler
#' @noRd
.rxFinalizeNls <- function(.s, sum.prod = FALSE,
                           optExpression = TRUE) {
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .yj <- paste(get("rx_yj_", envir = .s))
  .yj <- paste0("rx_yj_~", rxode2::rxFromSE(.yj))
  .lambda <- paste(get("rx_lambda_", envir = .s))
  .lambda <- paste0("rx_lambda_~", rxode2::rxFromSE(.lambda))
  .hi <- paste(get("rx_hi_", envir = .s))
  .hi <- paste0("rx_hi_~", rxode2::rxFromSE(.hi))
  .low <- paste(get("rx_low_", envir = .s))
  .low <- paste0("rx_low_~", rxode2::rxFromSE(.low))
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- character(0)
  .sens <- .s$..sens
  if (is.null(.sens)) .sens <- character(0)
  .s$..nlsS <- paste(c(
    .s$params,
    .s$..stateInfo["state"],
    .ddt,
    .sens,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .s$..HdTheta,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")
  .lhs0 <- .s$..lhs0
  if (is.null(.lhs0)) .lhs0 <- ""
  .s$..pred.nolhs <- paste(c(
    .s$params,
    .s$..stateInfo["state"],
    .lhs0,
    .ddt,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")

  if (sum.prod) {
    .malert("stabilizing round off errors in nls gradient problem...")
    .s$..nlsS <- rxode2::rxSumProdModel(.s$..nlsS)
    .msuccess("done")
    .malert("stabilizing round off errors in nls pred-only problem...")
    .s$..pred.nolhs <- rxode2::rxSumProdModel(.s$..pred.nolhs)
    .msuccess("done")
  }
  if (optExpression) {
    .s$..nlsS <- rxode2::rxOptExpr(.s$..nlsS, "nls gradient")
    .s$..pred.nolhs <- rxode2::rxOptExpr(.s$..pred.nolhs, "nls pred-only")
  }
}

#' @export
rxUiGet.nlsEnv <- function(x, ...) {
  .s <- rxUiGet.nlsHdTheta(x, ...)
  .s$params <- rxUiGet.nlsParams(x, ...)
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizeNls(.s, .sumProd, .optExpression)
  .s$..outer <- NULL
  if (exists("..maxTheta", .s)) {
    .eventTheta <- rep(0L, .s$..maxTheta)
  } else {
    .eventTheta <- integer(0)
  }
  for (.v in .s$..eventVars) {
    .vars <- as.character(get(.v, envir = .s))
    .vars <- rxode2::rxGetModel(paste0("rx_lhs=", rxode2::rxFromSE(.vars)))$params
    for (.v2 in .vars) {
      .reg <- rex::rex(start, "THETA[", capture(any_numbers), "]", end)
      if (regexpr(.reg, .v2) != -1) {
        .num <- as.numeric(sub(.reg, "\\1", .v2))
        .eventTheta[.num] <- 1L
      }
    }
  }
  .s$.eventTheta <- .eventTheta
  .s
}
attr(rxUiGet.nlsEnv, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsSensModel <- function(x, ...) {
  .s <- rxUiGet.nlsEnv(x, ...)
  list(thetaGrad=rxode2::rxode2(.s$..nlsS),
       predOnly=rxode2::rxode2(.s$..pred.nolhs),
       eventTheta=.s$.eventTheta)
}


#' @export
rxUiGet.nlsParStart <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  setNames(lapply(.w, function(i){
    .ui$iniDf$est[i]
  }),
  .ui$iniDf$name[.w])
}

#' @export
rxUiGet.nlsParStartTheta <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  setNames(vapply(.w, function(i){
    .ui$iniDf$est[i]
  }, double(1), USE.NAMES = FALSE),
  paste0("THETA[", seq_along(.ui$iniDf$name[.w]), "]")
  )
}
attr(rxUiGet.nlsParStartTheta, "rstudio") <- c(`THETA[1]`=0.1)

#' @export
rxUiGet.nlsParams <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  paste0("params(", paste(c(paste0("THETA[", seq_along(.ui$iniDf$name[.w]), "]"), "DV"), collapse=", "), ")")
}
attr(rxUiGet.nlsParams, "rstudio") <- "params(THETA[1], DV)"

#' @export
rxUiGet.nlsParLower <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  setNames(vapply(.w, function(i){
    .ui$iniDf$lower[i]
  }, double(1), USE.NAMES=FALSE),
  .ui$iniDf$name[.w])
}
attr(rxUiGet.nlsParLower, "rstudio") <- c(`ka`=0.01)

#' @export
rxUiGet.nlsParUpper <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  setNames(vapply(.w, function(i){
    .ui$iniDf$upper[i]
  }, double(1), USE.NAMES=FALSE),
  .ui$iniDf$name[.w])
}
attr(rxUiGet.nlsParUpper, "rstudio") <- c(`ka`=1000)

#' @export
rxUiGet.nlsParNameFun <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .args <- vapply(seq_along(.iniDf$ntheta),
                  function(t) {
                    if (.iniDf$err[t] %in% c("add", "prop", "pow")) {
                      ""
                    } else if (.iniDf$fix[t]) {
                      ""
                    } else {
                      .iniDf$name[t]
                    }
                  }, character(1), USE.NAMES=FALSE)
  .args <- .args[.args != ""]
  eval(str2lang(
    paste0("function(",paste(.args, collapse=", "),
           ") {c(",
           paste(vapply(seq_along(.iniDf$ntheta), function(t) {
             if (.iniDf$err[t] %in% c("add", "prop", "pow")) {
               paste0("'THETA[", t, "]'=", .iniDf$est[t])
             } else if (.iniDf$fix[t]) {
               paste0("'THETA[", t, "]'=", .iniDf$est[t])
             } else {
               paste0("'THETA[", t, "]'=", .iniDf$name[t])
             }
           }, character(1), USE.NAMES=FALSE), collapse=","), ")}")))
}
attr(rxUiGet.nlsParNameFun, "rstudio") <- function() {}

.nlsFormulaArgs <- function(x) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .args <- vapply(seq_along(.iniDf$ntheta),
                  function(t) {
                    if (.iniDf$err[t] %in% c("add", "prop", "pow")) {
                      ""
                    } else if (.iniDf$fix[t]) {
                      ""
                    } else {
                      .iniDf$name[t]
                    }
                  }, character(1), USE.NAMES=FALSE)
  c("DV", .args[.args != ""])
}

#' @export
rxUiGet.nlsFormula <- function(x, ..., grad=FALSE) {
  .args <- .nlsFormulaArgs(x)
  str2lang(paste0("~nlmixr2est::.nlmixrNlsFunValGrad(",
                  paste(.args, collapse=", "),
                  ")"))
}
attr(rxUiGet.nlsFormula, "rstudio") <- quote(~nlmixr2est::.nlmixrNlsFunValGrad(DV, ka, V, CL))

#' Setup the data for nls estimation
#'
#' @param dataSav Formatted Data
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlsFitDataSetup <- function(dataSav) {
  .dsAll <- dataSav[dataSav$EVID != 2, ] # Drop EVID=2 for estimation
  if (any(names(.dsAll) == "CENS")) {
    if (!all(.dsAll$CENS == 0)) {
      stop("'nls' does not work with censored data", call. =FALSE)
    }
  }
  nlmixr2global$nlsEnv$dataNls <- .dsAll[.dsAll$EVID == 0, ] # only observations are passed to nls
  nlmixr2global$nlsEnv$data <- rxode2::etTrans(.dsAll, nlmixr2global$nlsEnv$model)
}

.nlsFitModel <- function(ui, dataSav) {
  .ctl <- ui$control
  if (.ctl$solveType != "fun") {
    .mi <- ui$nlsSensModel
    .ctl$solveType <- 10L
  } else {
    .mi <- ui$nlsRxModel
    .ctl$solveType <- 11L
    .ctl$gradTo <- 0.0
  }
  if (is.null(.ctl$scaleC)) {
    .ctl$scaleC <- ui$scaleCnls
  }
  .p <- unlist(ui$nlsParStart)
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl,
                       lower=ui$nlsParLower, upper=ui$nlsParUpper)
  .env$par.ini.list <- setNames(as.list(.env$par.ini), names(ui$nlsParStart))

  if (.ctl$algorithm == "LM") {
    .nls.control <- minpack.lm::nls.lm.control(ftol = .ctl$ftol, ptol = .ctl$ptol,
                                               gtol = .ctl$gtol,
                                               diag = .ctl$diag,
                                               epsfcn = .ctl$epsfcn,
                                               factor = .ctl$factor,
                                               maxfev = .ctl$maxfev,
                                               maxiter = .ctl$maxiter,
                                               nprint = .ctl$nprint)
    class(.ctl) <- NULL
    if (.ctl$solveType == 11L) {
      .ret <- bquote(minpack.lm::nls.lm(
        par=.(.env$par.ini),
        lower=.(.env$lower),
        upper=.(.env$upper),
        fn=nlmixr2est::.nlmixrNlsFunVal,
        control=.(.nls.control)))
    } else {
      .ret <- bquote(minpack.lm::nls.lm(
        par=.(.env$par.ini),
        lower=.(.env$lower),
        upper=.(.env$upper),
        fn=nlmixr2est::.nlmixrNlsFunVal,
        jac=nlmixr2est::.nlmixrNlsFunGrad,
        control=.(.nls.control)))
    }
    .ret <- eval(.ret)
    .ret <- .nlmFinalizeList(.env, .ret, par="par", printLine=TRUE,
                             hessianCov=TRUE)
    .ret$sd <- sd(.ret$fvec)
    .ret$logLik <- sum(stats::dnorm(.ret$fvec, log=TRUE))
  } else {
    nlmixr2global$nlsEnv$dataNls <- dataSav[dataSav$EVID == 0, ]
    .nls.control <- stats::nls.control(
      maxiter = .ctl$maxiter, tol = .ctl$tol, minFactor = .ctl$minFactor,
      printEval = .ctl$printEval, warnOnly = .ctl$warnOnly,
      scaleOffset = .ctl$scaleOffset,
      nDcentral = .ctl$nDcentral)
    class(.ctl) <- NULL
    .ret <- bquote(stats::nls(
      formula= .(ui$nlsFormula),
      data=nlmixr2est::.nlmixrNlsData(),
      start=.(.env$par.ini.list),
      control=.(.nls.control),
      algorithm=.(.ctl$algorithm),
      trace=.(.ctl$trace), # nolint
      model=FALSE,
      lower=.(.env$lower),
      upper=.(.env$upper)
    ))
    .ret <- eval(.ret)
    .ret <- .nlmFinalizeList(.env, .ret, printLine=TRUE,
                             hessianCov=FALSE)
  }
  .ret
}

.nlsGetTheta <- function(nls, ui) {
  .iniDf <- ui$iniDf
  .theta0 <- nls$par
  if (inherits(nls, "nls.lm")) {
    .sd <- nls$sd
  } else {
    .sd <- sd(resid(nls))
  }
  setNames(vapply(seq_along(.iniDf$ntheta), function(t) {
    if (.iniDf$err[t] %in% c("add", "prop", "pow")) {
      .sd
    } else if (.iniDf$fix[t]) {
      .iniDf$est[t]
    } else {
      .theta0[.iniDf$name[t]]
    }
  }, double(1), USE.NAMES=FALSE), .iniDf$name)
}

.nlsControlToFoceiControl <- function(env, assign=TRUE) {
  .nlsControl <- env$nlsControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$nlsControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.nlsControl$sumProd,
                                optExpression=.nlsControl$optExpression,
                                literalFix=.nlsControl$literalFix,
                                literalFixRes=FALSE,
                                scaleTo=0,
                                calcTables=.nlsControl$calcTables,
                                addProp=.nlsControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.nlsControl$compress,
                                ci=.nlsControl$ci,
                                sigdigTable=.nlsControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

.nlsFamilyFit <- function(env, ...) {
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
  .nls <- .collectWarn(.nlsFitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$nls <- .nls[[1]]
  .ret$parHistData <- .ret$nls$parHistData
  .ret$nls$parHistData <- NULL

  .ret$message <- NULL
  if (rxode2::rxGetControl(.ui, "returnNls", FALSE)) {
    return(.ret$nls)
  }
  .ret$cov <- .ret$nls$cov
  if (inherits(.ret$nls, "nls.lm")) {
    .ret$message <- .ret$nls$message
    .ret$cov <- .ret$nls$cov
    .ret$covMethod <- paste0(.ret$nls$covMethod, " (LM)")
    .ret$objective <- -2 * .ret$nls$logLik
  } else {
    .ret$message <- .ret$nls$convInfo$stopMessage
    .ret$covMethod <- "nls"
    .ret$objective <- -2 * as.numeric(logLik(.ret$nls))
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .nlsGetTheta(.ret$nls, .ui)
  #.ret$etaMat <- NULL
  #.ret$etaObf <- NULL
  #.ret$omega <- NULL
  .ret$control <- .control
  .ret$extra <- paste0(" with ", crayon::bold$yellow(.control$algorithm),  " algorithm")
  .nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui)) {
    rm(list="control", envir=.ui)
  }
  .ret$est <- "nls"
  # There is no parameter history for nlse
  .ret$model <- .ui$ebe
  .ret$ofvType <- "nls"
  .nlsControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="nls")
  .env <- .ret$env
  .env$method <- "nls"
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.nls <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'nls', try 'focei'",
                                   .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'nls'",
                                   .var.name=.ui$modelName)
  rxode2::assertRxUiSingleEndpoint(.ui, " for the estimation routine 'nls'",
                                   .var.name=.ui$modelName)
  rxode2::assertRxUiEstimatedResiduals(.ui, " for the estimation routine 'nls'",
                                       .var.name=.ui$modelName)
  rxode2::warnRxBounded(.ui, " which are ignored in 'nls'", .var.name=.ui$modelName)
  # No add+prop or add+pow
  # Single endpoint
  .nlsFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .nlsFamilyFit(env,  ...)
}
attr(nlmixr2Est.nls, "covPresent") <- TRUE
