#' nlmixr2 defaults controls for nlm
#'
#' @inheritParams iterPrintParams
#' @inheritParams stats::nlm
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @param covMethod "r" uses nlmixr2's `nlmixr2Hess()` for the hessian, or
#'   "nlm" uses the hessian from `stats::nlm(.., hessian=TRUE)`; defaults to
#'   "nlm" when using nlmixr2's hessian/gradient for solving.
#' @param returnNlm is a logical that allows a return of the `nlm`
#'   object
#' @param solveType controls whether `nlm` uses nlmixr2's analytical
#'   gradients (event-related parameters like lag time/duration/rate/F use
#'   Shi2021 finite differences instead): `"hessian"` builds a Hessian from
#'   the analytical gradient via finite differences, `"gradient"` supplies
#'   the gradient and lets `nlm` compute the finite-difference Hessian, and
#'   `"fun"` lets `nlm` compute both by finite differences.
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
#' @param sensMethod Method used to compute the ODE parameter sensitivities:
#'   `"default"` (the default) defers to the global option
#'   `getOption("nlmixr2est.adjoint")` (itself `"forward"` by default);
#'   `"forward"` uses the classic variational (forward) sensitivity ODEs;
#'   `"adjoint"` uses the in-engine discrete adjoint with the matching adjoint
#'   (`s`) method.
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
                       indTolRelax=TRUE,

                       eventType=c("central", "forward"),
                       shiErr=(.Machine$double.eps)^(1/3),
                       shi21maxFD=20L,

                       optimHessType=c("central", "forward"),
                       hessErr =(.Machine$double.eps)^(1/3),
                       shi21maxHess=20L,

                       censOption=c("gauss", "laplace"),

                       eventSens=c("jump", "fd"),

                       sensMethod=c("default", "forward", "adjoint"),

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
                       calcTables=TRUE, compress=FALSE,
                       covMethod=c("r", "nlm", ""),
                       adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL,
                       boundedTransform=TRUE, ...) {
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
  checkmate::assertLogical(boundedTransform, len=1, any.missing=FALSE)

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
  # censOption: FOCEI-family censored (M2/M3/M4) 2nd-derivative treatment -- "gauss" (historic
  # Gauss-Newton, default) or "laplace" (exact).  Accepted for a uniform interface but INERT for
  # NLM (its finite-difference Hessian already reflects censoring exactly); kept for alignment.
  if (checkmate::testIntegerish(censOption, len=1, lower=0, upper=1, any.missing=FALSE)) {
    censOption <- as.integer(censOption)
  } else {
    censOption <- setNames(c("gauss" = 0L, "laplace" = 1L)[match.arg(censOption)], NULL)
  }

  ## eventSens: "jump" routes dosing-parameter (alag/F/rate/dur) sensitivities
  ## through rxode2's analytic event jumps; "fd" uses the legacy path that misses them.
  eventSens <- match.arg(eventSens)

  ## sensMethod: "forward" builds the ODE parameter sensitivities the classic
  ## (variational) way; "adjoint" solves them with the in-engine discrete
  ## adjoint (matching s-method); "default" defers to
  ## getOption("nlmixr2est.adjoint").
  sensMethod <- match.arg(sensMethod)

  ## eventSens: "jump" routes dosing-parameter (alag/F/rate/dur) sensitivities
  ## through rxode2's analytic event jumps; "fd" uses the legacy path that misses them.
  eventSens <- match.arg(eventSens)

  ## sensMethod: "forward" builds the ODE parameter sensitivities the classic
  ## (variational) way; "adjoint" solves them with the in-engine discrete
  ## adjoint (matching s-method); "default" defers to
  ## getOption("nlmixr2est.adjoint").
  sensMethod <- match.arg(sensMethod)

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
               indTolRelax=indTolRelax,

               eventType=eventType,
               shiErr=shiErr,
               shi21maxFD=as.integer(shi21maxFD),

               optimHessType=optimHessType,
               hessErr=hessErr,
               shi21maxHess=as.integer(shi21maxHess),
               censOption=censOption,

               eventSens=eventSens,
               sensMethod=sensMethod,

               iterPrintControl = .iterPrintControl,
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
               genRxControl=.genRxControl,
               boundedTransform=boundedTransform)
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
  .nlmFamilyControlGeneric(env, nlmixr2est::nlmControl, "nlmControl")
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
  ## Declare the model covariates (ui$allCovs) explicitly after DV.  Referenced
  ## covariates land here anyway (auto-detected after the thetas), so this only
  ## pins the order -- but it also keeps covariates whose only reference is dropped
  ## by log-likelihood pruning (e.g. a plugin's externally-loaded parameter block,
  ## read by a compiled function at a fixed par_ptr index rather than by name) in
  ## the solve parameter layout so they retain a stable par_ptr slot.
  .covs <- .ui$allCovs
  if (is.null(.covs)) .covs <- character(0)
  paste0("params(",
         paste(c(vapply(.w, function(i) {
           .env$t <- .env$t + 1
           paste0("THETA[", .env$t, "]")
         }, character(1), USE.NAMES = FALSE), "DV", .covs),
         collapse=","), ")")
}
attr(rxUiGet.nlmParams, "rstudio") <- "params()"

#' Extract rx_pred_f_ and rx_r_ model lines from symengine environment
#'
#' @param .s symengine environment
#' @return named list with `f_line` and `r_line` character strings
#' @noRd
.nlmGetFRLines <- function(.s) {
  .f_line <- ""
  .r_line <- ""
  if (exists("rx_pred_f_", envir = .s, inherits = FALSE)) {
    .f <- get("rx_pred_f_", envir = .s)
    .f_line <- paste0("rx_pred_f_=", rxode2::rxFromSE(.f))
  }
  if (exists("rx_r_", envir = .s, inherits = FALSE)) {
    .r <- get("rx_r_", envir = .s)
    .r_line <- paste0("rx_r_=", rxode2::rxFromSE(.r))
  }
  list(f_line = .f_line, r_line = .r_line)
}

#' @export
rxUiGet.nlmRxModel <- function(x, ...) {
  .s <- rxUiGet.loadPruneNlm(x, ...)
  # For matExp() models materialize the implied d/dt() from the k_from_to rate
  # constants.  When this fires we must also emit the model LHS (which defines
  # the k_from_to constants and other assignments) ahead of the d/dt() lines so
  # the derivative expressions can resolve them.
  .isMatExp <- isTRUE(.rxInjectMatExpDdt(.s))
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  ## .lhs0 <- .s$..lhs0
  ## if (is.null(.lhs0)) .lhs0 <- ""
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- ""
  .lhs <- character(0)
  if (.isMatExp) {
    .lhs <- .s$..lhs
    if (is.null(.lhs)) .lhs <- character(0)
  }
  # variables referenced by lag()/history functions (eg the AR(1) residual) are
  # not part of rx_pred_ itself; include their definitions so the history
  # reference resolves in the compiled model
  .lagDefs <- character(0)
  if (!is.null(.s$..laggedVars) && length(.s$..laggedVars) > 0L && !is.null(.s$..lhs)) {
    .pat <- paste0("^(", paste0(.s$..laggedVars, collapse = "|"), ")=")
    .lagDefs <- .s$..lhs[grepl(.pat, .s$..lhs)]
  }
  # Add rx_pred_f_ and rx_r_ as lhs outputs for censoring support
  .fr <- .nlmGetFRLines(.s)
  .ret <- paste(c(
    #.s$..stateInfo["state"],
    #.lhs0,
    .lhs,
    .ddt,
    .lagDefs,
    ## DDE non-constant delay() pre-history (base past(state,tau)<-expr)
    rxode2::.rxPastBaseLinesFromEnv(.s),
    .prd,
    .fr$f_line,
    .fr$r_line,
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
    .ret <- rxode2::rxOptExpr(.ret, "population log-likelihood model",
                              parallel = .optExprCores(x[[1]]))
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

#' Resolve the nlm/focei sensitivity method (forward vs adjoint) + s-method
#'
#' Reads `sensMethod` from the control (`"default"` defers to
#' `getOption("nlmixr2est.adjoint")`).  When the adjoint is selected it maps the
#' base ODE method to its in-engine discrete-adjoint (`s`) variant (base code +
#' 200, falling back to `dop853s`) and flags whether that variant needs the
#' stiff analytic Jacobian.
#'
#' @param ui rxode2 UI environment (i.e. `x[[1]]`)
#' @param nParam number of parameters differentiated by the sensitivities,
#'   recorded on the returned list.  Defaults to the estimated THETA count (nlm
#'   family); the focei inner path passes the ETA count instead.
#' @return list with `useAdjoint`; when TRUE also `sMethodInt`, `sMethodName`,
#'   `stiff`, `nParam`, `nState`.
#' @author Matthew L. Fidler
#' @noRd
.nlmAdjointResolve <- function(ui, nParam = NULL) {
  .sensMethod <- rxode2::rxGetControl(ui, "sensMethod", "default")
  ## "default" (the control default when sensMethod is not specified directly)
  ## defers to the global policy option so the package-wide default can be
  ## changed in one place: getOption("nlmixr2est.adjoint").
  if (identical(.sensMethod, "default")) {
    .sensMethod <- getOption("nlmixr2est.adjoint", "forward")
    if (!(.sensMethod %in% c("forward", "adjoint"))) .sensMethod <- "forward"
  }
  ## ODE (d/dt) state count -- the adjoint expansion differentiates these; a
  ## linCmt()/algebraic model reports 0 here (its pseudo-compartments are solved
  ## analytically, not integrated) so the adjoint does not apply.
  .nState <- length(rxode2::rxStateOde(ui))
  # the discrete adjoint only applies to ODE-state sensitivities; models with no
  # ODE states (e.g. linCmt()/algebraic) fall back to the forward path.
  if (!identical(.sensMethod, "adjoint") || .nState == 0L) {
    return(list(useAdjoint = FALSE))
  }
  if (is.null(nParam)) {
    .iniDf <- ui$iniDf
    nParam <- sum(!.iniDf$fix & !is.na(.iniDf$ntheta))
  }
  .sm <- .nlmAdjointSMethod(ui)
  list(useAdjoint = TRUE, sMethodInt = .sm$sMethodInt, sMethodName = .sm$sMethodName,
       stiff = .sm$stiff, nParam = nParam, nState = .nState)
}

#' Map the base ODE method to its in-engine discrete-adjoint (`s`) variant
#'
#' Independent of the forward/adjoint decision: given the control's base method,
#' returns the adjoint `s`-method (base code + 200, falling back to `dop853s`)
#' and whether it needs the stiff analytic Jacobian.  Used both by
#' [.nlmAdjointResolve()] and, on the solve side, wherever an already-built
#' adjoint model (one carrying `rx__adjFX_*`) must be pointed at its s-method.
#'
#' @param ui rxode2 UI environment (i.e. `x[[1]]`)
#' @return list with `sMethodInt`, `sMethodName`, `stiff`.
#' @author Matthew L. Fidler
#' @noRd
.nlmAdjointSMethod <- function(ui) {
  .rxControl <- rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl())
  .baseInt <- suppressWarnings(as.integer(.rxControl$method))
  if (length(.baseInt) != 1L || is.na(.baseInt)) .baseInt <- 2L # liblsoda
  .adjIdx <- rxode2::odeMethodToInt(NULL)
  .adjCodes <- .adjIdx[.adjIdx >= 200L]
  if (.baseInt >= 200L) {
    .up <- .baseInt
  } else {
    .up <- .baseInt + 200L
    if (!(.up %in% .adjCodes)) .up <- 200L # no direct variant -> dop853s
  }
  .sName <- names(.adjIdx)[match(.up, .adjIdx)][1]
  list(sMethodInt = .up, sMethodName = .sName,
       stiff = isTRUE(rxode2::.rxAdjointMethodStiff(.sName)))
}

#' Build the adjoint replacement for the forward variational sensitivity block
#'
#' Runs rxode2's adjoint expansion on the (THETA/ETA-substituted) model and
#' keeps only the pieces that replace the forward `.sens` lines: the
#' `rx__sens_<state>_BY_<param>__` output compartments, the analytic `df()/dy()`
#' Jacobian (stiff `s`-methods only), and the `rx__adj*` sweep lhs.  The base
#' ODE and dosing modifiers are dropped -- the nlm/focei assembly already emits
#' those.  Used for both the nlm THETA gradient and the focei inner ETA
#' sensitivities (same output structure, different `calcSens`).
#'
#' A build-time parity guard asserts the emitted sensitivity-column set matches
#' the forward layout exactly, since the C++ reads those columns by name.
#'
#' @param object THETA/ETA-substituted (pruned) model text
#' @param calcSens the parameter names differentiated by (THETA_n or ETA_n)
#' @param stiff whether the chosen s-method needs the analytic Jacobian
#' @return character vector of model lines to splice in place of `.sens`
#' @author Matthew L. Fidler
#' @noRd
.rxAdjointSensLines <- function(object, calcSens, stiff) {
  .ex <- rxode2::.rxAdjointExpand(object, calcSens, stiff = isTRUE(stiff))
  .lines <- strsplit(.ex$text, "\n", fixed = TRUE)[[1]]
  .keep <- grepl("^d/dt\\(rx__sens_", .lines) |
    grepl("^df\\(", .lines) |
    grepl("^rx__adj", .lines)
  .ret <- .lines[.keep]
  .want <- as.vector(vapply(calcSens, function(.p) {
    paste0("rx__sens_", .ex$st, "_BY_", .p, "__")
  }, character(length(.ex$st))))
  .got <- sub("^d/dt\\((rx__sens_.*__)\\)=.*$", "\\1",
              .ret[grepl("^d/dt\\(rx__sens_", .ret)])
  if (!setequal(.got, .want)) {
    stop("adjoint sensitivity columns do not match the forward layout", call. = FALSE) # nocov
  }
  .ret
}

#' Finalize nlm rxode2 based on symengine saved info
#'
#' @param .s Symengine/rxode2 object
#' @return Nothing
#' @author Matthew L Fidler
#' @noRd
.rxFinalizeNlm <- function(.s, sum.prod = FALSE,
                           optExpression = TRUE, cores = 0L) {
  .rxInjectMatExpDdt(.s)
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
  .lhs <- .s$..lhs
  if (is.null(.lhs)) .lhs <- character(0)
  .sens <- .s$..sens
  if (is.null(.sens)) .sens <- character(0)
  # sensMethod="adjoint": swap the forward variational block for the adjoint one
  # (rx__sens_* output compartments + df/dy + rx__adj* lhs).  The nlm.cpp read
  # is by-name, so this is transparent to the C++ gradient path; the adjoint
  # s-method is selected on the solve side (.nlmFitModel).
  if (!is.null(.s$..adjSens)) .sens <- .s$..adjSens
  # Extract rx_pred_f_ and rx_r_ for censoring support
  .fr <- .nlmGetFRLines(.s)
  .s$..nlmS <- paste(c(
    .s$params,
    .s$..stateInfo["state"],
    .lhs,
    .ddt,
    .sens,
    ## DDE non-constant delay() pre-history: base past(state,tau)<-expr + the
    ## per-sensitivity-compartment histories (analytic nlm gradient/Hessian).
    .s$..pastLines,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .s$..HdTheta,
    .fr$f_line,
    .fr$r_line,
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
    .lhs,
    .ddt,
    ## DDE non-constant delay() pre-history (base past(state,tau)<-expr; the
    ## pred-only model has no sensitivity compartments)
    .s$..pastBaseLines,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .fr$f_line,
    .fr$r_line,
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
    .s$..nlmS <- rxode2::rxOptExpr(.s$..nlmS, "nlm llik gradient", parallel = cores)
    .s$..pred.nolhs <- rxode2::rxOptExpr(.s$..pred.nolhs, "nlm pred-only", parallel = cores)
  }
}

#' @export
rxUiGet.nlmEnv <- function(x, ...) {
  .s <- rxUiGet.nlmHdTheta(x, ...)
  .s$params <- rxUiGet.nlmParams(x, ...)
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  # When sensMethod resolves to "adjoint", build the adjoint sensitivity block
  # against the SAME THETA_n parameter names the forward `.rxSens` uses, so the
  # rx__sens_<state>_BY_THETA_n__ output columns line up name-for-name.
  .adj <- .nlmAdjointResolve(x[[1]])
  if (isTRUE(.adj$useAdjoint) && exists("..maxTheta", .s) && .s$..maxTheta > 0L) {
    .calcSens <- paste0("THETA_", seq_len(.s$..maxTheta), "_")
    .s$..adjSens <- .rxAdjointSensLines(.nlmPrune(x), .calcSens, .adj$stiff)
  }
  .rxFinalizeNlm(.s, .sumProd, .optExpression, .optExprCores(x[[1]]))
  .s$..outer <- NULL
  if (exists("..maxTheta", .s)) {
    .eventTheta <- rep(0L, .s$..maxTheta)
  } else {
    .eventTheta <- integer(0)
  }
  ## eventTheta flags dosing-parameter (alag/F/rate/dur) THETAs; under "fd" nlm
  ## overrides their gradient with finite differences, under "jump" it's left analytic since rxode2 injects the jump directly.
  .eventSens <- rxode2::rxGetControl(x[[1]], "eventSens", "jump")
  if (!identical(.eventSens, "jump")) {
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
  ## "jump" attaches rxode2's analytic event (alag/F/rate/dur) sensitivities to
  ## the thetaGrad model; nlm has no FD fallback so under "fd" the gradient simply misses the jump.
  .eventSens <- rxode2::rxGetControl(x[[1]], "eventSens", "jump")
  ## adjoint carries its own dosing-parameter corrections (rx__adjdF/Dlag/Drate)
  ## in the sweep, so it must NOT also inject the forward variational jump.
  if (!is.null(.s$..adjSens)) .eventSens <- "fd"
  list(thetaGrad=rxode2::rxode2(.s$..nlmS, eventSens=.eventSens),
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
  nlmixr2global$nlmEnv$data <- rxode2::etTrans(.dsAll, nlmixr2global$nlmEnv$model)
}

#' Set up an nlm-family objective for repeated hook-firing evaluation
#'
#' Preprocesses the data and LOADS the nlm population (predOnly) problem into the
#' C++ engine, so that repeated \code{nlmSolveR(theta)} calls evaluate the
#' population objective -- firing any registered likelihood-contribution hook
#' (e.g. a plugin's per-observation cotangent capture) -- WITHOUT re-running the
#' optimizer.  One compiled setup is reused across evaluations.  Intended for
#' extension packages (e.g. nlmixr2nn) that optimize an out-of-band parameter
#' block (network weights injected via a par-loader) and need the exact
#' error-model cotangent from the nlm C++ solve at each iterate.  Free the loaded
#' problem with \code{.nlmFreeEnv()} when done.
#'
#' @param ui rxode2/nlmixr2 model.  Uses \code{ui$control} when present.
#' @param data event data.
#' @param control optional nlm-family control; defaults to \code{nlmControl()}
#'   (or \code{ui$control} if that is an nlm-family control).
#' @return (invisibly) the scaled starting parameter vector to hand to
#'   \code{nlmSolveR()}; the C++ problem is left loaded.
#' @export
#' @keywords internal
#' @author Matthew L. Fidler
nlmObjectiveSetup <- function(ui, data, control = NULL) {
  .ui <- rxode2::rxUiDecompress(ui)
  if (is.null(control)) {
    control <- if (!is.null(.ui$control)) .ui$control else nlmControl()
  }
  .ui$control <- control
  .ctl <- .ui$control
  class(.ctl) <- NULL
  .ret <- new.env(parent = emptyenv())
  .foceiPreProcessData(data, .ret, .ui, .ctl$rxControl)
  .p <- setNames(.ui$nlmParIni, .ui$nlmParName)
  ## solveType 1 / nlmRxModel: the objective-only predOnly model (no thetaGrad).
  ## The hook fires from nlmSolveFid during the objective solve; the caller gets
  ## the weight gradient from its own augmented-sensitivity solve, so no analytic
  ## theta gradient is needed here.
  .env <- .nlmSetupEnv(.p, .ui, .ret$dataSav, .ui$nlmRxModel, .ctl)
  invisible(.env$par.ini)
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
  ## Event ("jump") sensitivities are activated in .nlmSetupEnv and deactivated in .nlmFreeEnv; nothing extra needed here.
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
                                sigdigTable=.nlmControl$sigdigTable,
                                indTolRelax=.nlmControl$indTolRelax,
                                eventSens=.nlmControl$eventSens)
  if (assign) env$control <- .foceiControl
  .foceiControl
}


.nlmFamilyFit <- function(env, ...) {
  .nlmFamilyFitGeneric(
    env, "nlm", .nlmFitModel, .nlmGetTheta,
    objective = function(.fit) 2 * as.numeric(.fit$minimum),
    controlToFocei = .nlmControlToFoceiControl,
    returnFlag = "returnNlm",
    emitFitWarnings = TRUE,
    message = function(.fit) {
      if (.fit$code == 1) {
        "relative gradient is close to zero, current iterate is probably solution"
      } else if (.fit$code == 2) {
        "successive iterates within tolerance, current iterate is probably solution"
      } else if (.fit$code == 3) {
        c("last global step failed to locate a point lower than 'estimate'",
          "either 'estimate' is an approximate local minimum of the function or 'steptol' is too small")
      } else if (.fit$code == 4) {
        "iteration limit exceeded"
      } else if (.fit$code == 5) {
        c("maximum step size 'stepmax' exceeded five consecutive times",
          "either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or 'stepmax' is too small")
      } else {
        ""
      }
    })
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
attr(nlmixr2Est.nlm, "unbounded") <- TRUE
