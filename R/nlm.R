#' nlmixr2 defaults controls for nlm
#'
#' @inheritParams stats::nlm
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @param covMethod allows selection of "r", which uses nlmixr2's
#'   `nlmixr2Hess()` for the hessian calculation or "nlm" which uses
#'   the hessian from `stats::nlm(.., hessian=TRUE)`. When using
#'   `nlmixr2's` hessian for optimization or `nlmixr2's` gradient for
#'   solving this defaults to "nlm" since `stats::optimHess()` assumes
#'   an accurate gradient and is faster than `nlmixr2Hess`
#' @param returnNlm is a logical that allows a return of the `nlm`
#'   object
#' @param solveType tells if `nlm` will use nlmixr2's analytical
#'   gradients when available (finite differences will be used for
#'   event-related parameters like parameters controlling lag time,
#'   duration/rate of infusion, and modeled bioavailability). This can
#'   be:
#'
#'  - `"hessian"` which will use the analytical gradients to create a
#'     Hessian with finite differences.
#'
#' - `"gradient"` which will use the gradient and let `nlm` calculate
#'    the finite difference hessian
#'
#' - `"fun"` where nlm will calculate both the finite difference
#'    gradient and the finite difference Hessian
#'
#'  When using nlmixr2's finite differences, the "ideal" step size for
#'  either central or forward differences are optimized for with the
#'  Shi2021 method which may give more accurate derivatives
#'
#' @param shiErr This represents the epsilon when optimizing the ideal
#'   step size for numeric differentiation using the Shi2021 method
#'
#' @param hessErr This represents the epsilon when optimizing the
#'   Hessian step size using the Shi2021 method.
#'
#' @param shi21maxHess Maximum number of times to optimize the best
#'   step size for the hessian calculation
#'
#' @param gradTo this is the factor that the gradient is scaled to
#'   before optimizing.  This only works with
#'   scaleType="nlmixr2".
#'
#' @return nlm control object
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
#' fit2 <- nlmixr(mod, dsn, est="nlm")
#'
#' print(fit2)
#'
#' # you can also get the nlm output with fit2$nlm
#'
#' fit2$nlm
#'
#' # The nlm control has been modified slightly to include
#' # extra components and name the parameters
#' }
nlmControl <- function(typsize = NULL,
                       fscale = 1, print.level = 0, ndigit = NULL, gradtol = 1e-6,
                       stepmax = NULL,
                       steptol = 1e-6, iterlim = 10000, check.analyticals = FALSE,
                       returnNlm=FALSE,
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

                       rxControl=NULL,
                       optExpression=TRUE, sumProd=FALSE,
                       literalFix=TRUE,
                       literalFixRes=TRUE,
                       addProp = c("combined2", "combined1"),
                       calcTables=TRUE, compress=FALSE,
                       covMethod=c("r", "nlm", ""),
                       adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {
  checkmate::assertNumeric(shiErr, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(hessErr, lower=0, any.missing=FALSE, len=1)

  checkmate::assertIntegerish(shi21maxFD, lower=1, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(shi21maxHess, lower=1, any.missing=FALSE, len=1)

  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFixRes, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertNumeric(stepmax, lower=0, len=1, null.ok=TRUE, any.missing=FALSE)
  checkmate::assertIntegerish(print.level, lower=0, upper=2, any.missing=FALSE)
  checkmate::assertNumeric(ndigit, lower=0, len=1, any.missing=FALSE, null.ok=TRUE)
  checkmate::assertNumeric(gradtol, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(steptol, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(iterlim, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(check.analyticals, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnNlm, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %fin% c("genRxControl"))]
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

  .solveTypeIdx <- c("hessian" = 3L, "grad" = 2L, "fun" = 1L)
  if (checkmate::testIntegerish(solveType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    solveType <- as.integer(solveType)
  } else {
    solveType <- setNames(.solveTypeIdx[match.arg(solveType)], NULL)
  }
  if (missing(covMethod) && any(solveType == 2:3)) {
    covMethod <- "nlm"
  } else {
    covMethod <- match.arg(covMethod)
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

  .ret <- list(covMethod=covMethod,
               typsize = typsize,
               fscale = fscale, print.level = print.level, ndigit=ndigit, gradtol = gradtol,
               stepmax = stepmax,
               steptol = steptol, iterlim = iterlim,
               check.analyticals = check.analyticals,
               optExpression=optExpression,
               literalFix=literalFix,
               literalFixRes=literalFixRes,
               sumProd=sumProd,
               rxControl=rxControl,
               returnNlm=returnNlm,

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

               addProp=match.arg(addProp),
               calcTables=calcTables,
               compress=compress,
               solveType=solveType,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "nlmControl"
  .ret
}

#' @export
rxUiDeparse.nlmControl <- function(object, var) {
  .default <- nlmControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}


#' Get the nlm family control
#'
#' @param env nlm optimization environment
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

#' A surrogate function for nlm to call for ode solving
#'
#' @param pars Parameters that will be estimated
#' @return Predictions
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.nlmixrNlmFunC <- function(pars) {
  .Call(`_nlmixr2est_nlmSolveSwitch`, pars)
}

#' Get the THETA lines from rxode2 UI and assign fixed
#'
#' @param rxui This is the rxode2 ui object
#' @return The theta/eta lines
#' @author Matthew L. Fidler
#' @noRd
.uiGetThetaDropFixed <- function(rxui) {
  .iniDf <- rxui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  .env <- new.env(parent=emptyenv())
  .env$t <- 0
  lapply(.w, function(i) {
    if (.iniDf$fix[i]) {
      eval(str2lang(paste0("quote(", .iniDf$name[i], " <- ", .iniDf$est[i], ")")))
    } else {
      .env$t <- .env$t + 1
      eval(str2lang(paste0("quote(", .iniDf$name[i], " <- THETA[", .env$t, "])")))
    }
  })
}

#'@export
rxUiGet.nlmModel0 <- function(x, ...) {
  .ui <- rxode2::rxUiDecompress(x[[1]])
  nlmixr2global$rxPredLlik <- TRUE
  on.exit(nlmixr2global$rxPredLlik <- FALSE)
  .predDf <- .ui$predDf
  .save <- .predDf
  .predDf[.predDf$distribution == "norm", "distribution"] <- "dnorm"
  assign(".predDfFocei", .predDf, envir=.ui)
  #assign("predDf", .predDf, envir=.ui)
  on.exit(assign("predDf", .save, envir=.ui))
  .ret <- rxode2::rxCombineErrorLines(.ui, errLines=rxGetDistributionFoceiLines(.ui),
                              prefixLines=.uiGetThetaDropFixed(.ui),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE)
  .ret <- .ret[[2]]
  .ret <- as.call(c(quote(`{`),
                    lapply(seq_along(.ret)[-1], function(i) {
                      .ret[[i]]
                    }),
                    list(str2lang("rx_pred_ <- -rx_pred_"))))
  as.call(c(list(quote(`rxModelVars`)), .ret))
}
attr(rxUiGet.nlmModel0, "rstudio") <- quote(rxModelVar({}))

#' Load the nlm model into symengine
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
  .malert("pruning branches ({.code if}/{.code else}) of population log-likelihood model...")
  .ret <- rxode2::.rxPrune(.x, envir = .env,
                           strAssign = rxModelVars(x[[1]])$strAssign)
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
  .p <- .nlmPrune(x)
  .loadSymengine(.p, promoteLinSens = FALSE)
}
attr(rxUiGet.loadPruneNlm, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlmParams <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .w <- which(!.iniDf$fix)
  .env <- new.env(parent=emptyenv())
  .env$t <- 0
  paste0("params(",
         paste(c(vapply(.w, function(i) {
           .env$t <- .env$t + 1
           paste0("THETA[", .env$t, "]")
         }, character(1), USE.NAMES = FALSE), "DV"),
         collapse=","), ")")
}
attr(rxUiGet.nlmParams, "rstudio") <- "params()"

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
    .malert("stabilizing round off errors in population log-likelihood model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "population log-likelihood model")
    .msuccess("done")
  }
  .cmt <-  rxUiGet.foceiCmtPreModel(x, ...)
  .interp <- rxUiGet.interpLinesStr(x, ...)
  if (.interp != "") {
    .cmt <-paste0(.cmt, "\n", .interp)
  }
  list(predOnly=rxode2::rxode2(paste(c(rxUiGet.nlmParams(x, ...), .cmt,
                                       .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")),
       eventTheta=.eventTheta)
}

#' @export
rxUiGet.loadPruneNlmSens <- function(x, ...) {
  .loadSymengine(.nlmPrune(x), promoteLinSens = TRUE)
}
attr(rxUiGet.loadPruneNlmSens, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlmThetaS <- function(x, ...) {
  .s <- rxUiGet.loadPruneNlmSens(x, ...)
  .sensEtaOrTheta(.s, theta=TRUE)
}
attr(rxUiGet.nlmThetaS, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlmHdTheta <- function(x, ...) {
  .s <- rxUiGet.nlmThetaS(x)
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
attr(rxUiGet.nlmHdTheta, "rstudio") <- emptyenv()

#' Finalize nlm rxode2 based on symengine saved info
#'
#' @param .s Symengine/rxode2 object
#' @return Nothing
#' @author Matthew L Fidler
#' @noRd
.rxFinalizeNlm <- function(.s, sum.prod = FALSE,
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
  .s$..nlmS <- paste(c(
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
    .malert("stabilizing round off errors in nlm llik gradient problem...")
    .s$..nlmS <- rxode2::rxSumProdModel(.s$..nlmS)
    .msuccess("done")
    .malert("stabilizing round off errors in nlm llik pred-only problem...")
    .s$..pred.nolhs <- rxode2::rxSumProdModel(.s$..pred.nolhs)
    .msuccess("done")
  }
  if (optExpression) {
    .s$..nlmS <- rxode2::rxOptExpr(.s$..nlmS, "nlm llik gradient")
    .s$..pred.nolhs <- rxode2::rxOptExpr(.s$..pred.nolhs, "nlm pred-only")
  }
}

#' @export
rxUiGet.nlmEnv <- function(x, ...) {
  .s <- rxUiGet.nlmHdTheta(x, ...)
  .s$params <- rxUiGet.nlmParams(x, ...)
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizeNlm(.s, .sumProd, .optExpression)
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
  ## if (.sumProd) {
  ##   .malert("stabilizing round off errors in pred-only model...")
  ##   s$..pred.nolhs <- rxode2::rxSumProdModel(.s$..pred.nolhs)
  ##   .msuccess("done")
  ## }
  ## if (.optExpression) {
  ##   s$..pred.nolhs <- rxode2::rxOptExpr(.s$..pred.nolhs,
  ##                                       ifelse(.getRxPredLlikOption(),
  ##                                              "Llik pred-only model",
  ##                                              "pred-only model"))
  ## }
  ## s$..pred.nolhs <- paste(c(
  ##   paste0("params(", paste(inner$params, collapse = ","), ")"),
  ##   s$..pred.nolhs
  ## ), collapse = "\n")

  .s$.eventTheta <- .eventTheta

  .s
}
attr(rxUiGet.nlmEnv, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlmSensModel <- function(x, ...) {
  .s <- rxUiGet.nlmEnv(x, ...)
  list(thetaGrad=rxode2::rxode2(.s$..nlmS),
       predOnly=rxode2::rxode2(.s$..pred.nolhs),
       eventTheta=.s$.eventTheta)
}

#' @export
rxUiGet.nlmParNameFun <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .env <- new.env(parent=emptyenv())
  .env$i <- 1
  .w <- which(!.iniDf$fix)
  eval(str2lang(
    paste0("function(p) {c(",
           paste(vapply(.w, function(t) {
             .ret <- paste0("'THETA[", .env$i, "]'=p[", .env$i, "]")
             .env$i <- .env$i + 1
             .ret
           }, character(1), USE.NAMES=FALSE), collapse=","), ")}")))
}
attr(rxUiGet.nlmParNameFun, "rstudio") <- function(){c(`THETA[1]`=1, `THETA[2]`=2, `THETA[3]`=3)}

#' @export
rxUiGet.optimParNameFun <- rxUiGet.nlmParNameFun

#' @export
rxUiGet.nlmParIni <- function(x, ...) {
  .ui <- x[[1]]
  .ui$iniDf$est[!.ui$iniDf$fix]
}
attr(rxUiGet.nlmParIni, "rstudio") <- c(1, 2, 3)

#' @export
rxUiGet.optimParIni <- rxUiGet.nlmParIni

#' @export
rxUiGet.nlmParName <- function(x, ...) {
  .ui <- x[[1]]
  .ui$iniDf$name[!.ui$iniDf$fix]
}
attr(rxUiGet.nlmParName, "rstudio") <- c("THETA[1]", "THETA[2]", "THETA[3]")

#' @export
rxUiGet.optimParName <- rxUiGet.nlmParName

#' Setup the data for nlm estimation
#'
#' @param dataSav Formatted Data
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmFitDataSetup <- function(dataSav) {
  .dsAll <- dataSav[dataSav$EVID != 2, ] # Drop EVID=2 for estimation
  if (any(names(.dsAll) == "CENS")) {
    if (!all(.dsAll$CENS == 0)) {
      stop("'nlm' does not work with censored data", call. =FALSE)
    }
  }
  nlmixr2global$nlmEnv$data <- rxode2::etTrans(.dsAll, nlmixr2global$nlmEnv$model)
}

.nlmFitModel <- function(ui, dataSav) {
  .ctl <- ui$control
  class(.ctl) <- NULL
  .p <- setNames(ui$nlmParIni, ui$nlmParName)
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
  .hessian <- .ctl$covMethod == "nlm"
  if (.ctl$solveType == 1L) {
    .mi <-  ui$nlmRxModel
  } else {
    .mi <- ui$nlmSensModel
  }
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl)
  on.exit({.nlmFreeEnv()})
  .ret <- eval(bquote(stats::nlm(
    f=.(.nlmixrNlmFunC),
    p=.(.env$par.ini),
    hessian=.(.hessian),
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
  .nlmFinalizeList(.env, .ret, par="estimate", printLine=TRUE,
                   hessianCov=TRUE)
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
             .iniDf$est[i]
           } else {
             nlm$estimate[.iniDf$name[i]]
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
                                literalFix=.nlmControl$literalFix,
                                literalFixRes=.nlmControl$literalFixRes,
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
  .ret$parHistData <- .ret$nlm$parHistData
  .ret$nlm$parHistData <- NULL
  .ret$message <- NULL
  lapply(.nlm[[2]], function(x){
    warning(x, call.=FALSE)
  })

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
  .ret$objective <- 2 * as.numeric(.ret$nlm$minimum)
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
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'nlm', try 'focei'",
                                   .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'nlm'",
                                   .var.name=.ui$modelName)
  rxode2::warnRxBounded(.ui, " which are ignored in 'nlm'",
                        .var.name=.ui$modelName)
  .nlmFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .nlmFamilyFit(env,  ...)
}
attr(nlmixr2Est.nlm, "covPresent") <- TRUE
