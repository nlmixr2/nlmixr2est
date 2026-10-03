#' nlmixr2 defaults controls for nls
#'
#' @inheritParams iterPrintParams
#' @inheritParams stats::nls
#' @inheritParams stats::nls.control
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @inheritParams nlmControl
#' @inheritParams minpack.lm::nls.lm.control
#' @param returnNls logical; when TRUE, will return the nls object
#'   instead of the nlmixr object
#' @return nls control object
#'
#' @section Pooled (population-only) estimation:
#'
#' `est="nls"` is a *pooled* estimation method.  It estimates population
#' parameters only; there is no eta-conditional inner problem, and
#' sensitivities are taken with respect to the population parameters alone.
#'
#' A model carrying any random effect is therefore refused before fitting
#' rather than being fit to something other than what it says -- it stops
#' with `can only have population estimates for the estimation routine
#' 'nls', try 'focei'`.
#'
#' A successful `nls` fit consequently has no `$omega`, no empirical Bayes
#' estimates and no shrinkage, and its objective function table carries a
#' single `Pop` row.  Use a mixed effects routine -- `est="focei"`,
#' `est="saem"` and so on -- for any model with between-subject variability.
#'
#' Note that this help page inherits its parameter list from
#' `foceiControl()`, `saemControl()` and `nlmControl()`, because `nls`
#' shares many options with them.  Inherited options that only have meaning
#' for random effects do not apply to `nls`.
#'
#' @export
#' @author Matthew L. Fidler
#' @examples
#' \donttest{
#'
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45
#'     tcl <- log(c(0, 2.7, 100))
#'     tv <- 3.45
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka)
#'     cl <- exp(tcl)
#'     v <- exp(tv)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' # Note that `one.cmt` declares no random effects: `nls` is a pooled
#' # method and refuses a model that has any.
#'
#' # Uses nlsLM from minpack.lm if available
#'
#' fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est = "nls", nlsControl(algorithm = "LM"))
#'
#' # Uses port and respect parameter boundaries
#' fit2 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est = "nls", nlsControl(algorithm = "port"))
#'
#' # You can access the underlying nls object with `$nls`
#' fit2$nls
#' }
nlsControl <- function(
  maxiter = 10000,
  tol = NULL,
  minFactor = 1 / 1024,
  printEval = FALSE,
  warnOnly = FALSE,
  scaleOffset = 0,
  nDcentral = FALSE,
  algorithm = c("LM", "default", "plinear", "port"),
  ############################################
  ## minpack.lm
  ftol = NULL,
  ptol = NULL,
  gtol = 0,
  diag = list(),
  epsfcn = 0,
  factor = 100,
  maxfev = integer(),
  nprint = 0,
  #### nlm C++ style style to give gradients
  solveType = c("grad", "fun"),
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  indTolRelax = TRUE,
  eventType = c("central", "forward"),
  shiErr = (.Machine$double.eps)^(1 / 3),
  shi21maxFD = 20L,
  useColor = NULL,
  printNcol = NULL, #
  print = 1L, #
  normType = c("rescale2", "mean", "rescale", "std", "len", "constant"), #
  scaleType = c("nlmixr2", "norm", "mult", "multAdd"), #
  scaleCmax = 1e5, #
  scaleCmin = 1e-5, #
  scaleC = NULL,
  scaleTo = 1.0,
  gradTo = 1.0,
  ############################################
  trace = FALSE, # nolint
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  returnNls = FALSE,
  addProp = c("combined2", "combined1"),
  eventSens = c("jump", "fd"),
  linCmtSensCarry = c("auto", "none"),
  calcTables = TRUE,
  compress = TRUE,
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  boundedTransform = TRUE,
  ...
) {
  algorithm <- match.arg(algorithm)
  if (algorithm == "LM" && !requireNamespace("minpack.lm", quietly = TRUE)) {
    .malert("to use the LM algorithm you must have minpack.lm installed")
    .malert("changing to default `nls` method")
    algorithm <- "default"
  }
  checkmate::assertIntegerish(stickyRecalcN, any.missing = FALSE, lower = 0, len = 1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(odeRecalcFactor, len = 1, lower = 1, any.missing = FALSE)
  checkmate::assertLogical(indTolRelax, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(shiErr, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(shi21maxFD, lower = 1, any.missing = FALSE, len = 1)

  .eventTypeIdx <- c("central" = 2L, "forward" = 1L)
  if (checkmate::testIntegerish(eventType, len = 1, lower = 1, upper = 6, any.missing = FALSE)) {
    eventType <- as.integer(eventType)
  } else {
    eventType <- setNames(.eventTypeIdx[match.arg(eventType)], NULL)
  }

  solveType <- match.arg(solveType)

  # nls tolerances from sigdig: the LM/nls `tol` is 10^(-sigdig) (1e-4 at the
  # default sigdig=4); the minpack.lm function/parameter tolerances keep their
  # sqrt(eps) default at sigdig=4 and tighten one order per significant digit.  A
  # user value wins, sigdig=NULL keeps the historic defaults.
  if (is.null(ftol)) {
    ftol <- if (!is.null(sigdig)) .sigdigScale(sqrt(.Machine$double.eps), sigdig) else sqrt(.Machine$double.eps)
  }
  if (is.null(ptol)) {
    ptol <- if (!is.null(sigdig)) .sigdigScale(sqrt(.Machine$double.eps), sigdig) else sqrt(.Machine$double.eps)
  }
  checkmate::assertNumeric(ftol, len = 1, any.missing = FALSE, lower = 0)
  checkmate::assertNumeric(ptol, len = 1, any.missing = FALSE, lower = 0)
  checkmate::assertNumeric(gtol, len = 1, any.missing = FALSE, lower = 0)
  checkmate::assertNumeric(epsfcn, len = 1, any.missing = FALSE, lower = 0)
  checkmate::assertNumeric(factor, len = 1, any.missing = FALSE, lower = 1)
  checkmate::assertIntegerish(maxfev, min.len = 0, max.len = 1, any.missing = FALSE)

  checkmate::assertLogical(trace, len = 1, any.missing = FALSE) # nolint
  checkmate::assertLogical(nDcentral, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(scaleOffset, any.missing = FALSE, finite = TRUE)
  checkmate::assertLogical(warnOnly, len = 1, any.missing = FALSE)
  checkmate::assertLogical(printEval, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(maxiter, len = 1, any.missing = FALSE, lower = 1)
  if (is.null(tol)) {
    tol <- if (!is.null(sigdig)) .sigdigOptTol(sigdig) else 1e-05
  }
  checkmate::assertNumeric(tol, len = 1, any.missing = FALSE, lower = 0)
  checkmate::assertNumeric(minFactor, len = 1, any.missing = FALSE, lower = 0)
  checkmate::assertLogical(optExpression, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFix, len = 1, any.missing = FALSE)
  checkmate::assertLogical(sumProd, len = 1, any.missing = FALSE)
  checkmate::assertLogical(returnNls, len = 1, any.missing = FALSE)
  checkmate::assertLogical(calcTables, len = 1, any.missing = FALSE)
  checkmate::assertLogical(compress, len = 1, any.missing = TRUE)
  checkmate::assertLogical(adjObf, len = 1, any.missing = TRUE)
  checkmate::assertLogical(boundedTransform, len = 1, any.missing = FALSE)
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl", "iterPrintControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste(paste0("'", .bad, "'", sep = ""), collapse = ", "), call. = FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      # nls's Levenberg-Marquardt step is sensitive to ODE solver noise, so give it
      # a tighter solve (3 orders below the shared sigdig target) than the optimizer
      rxControl <- .rxControlScaleSigdig(rxode2::rxControl(sigdig = sigdig), sigdig, tighten = 3)
    } else {
      rxControl <- rxode2::rxControl(atol = 1e-4, rtol = 1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {} else if (is.list(rxControl)) {
    rxControl <- .rxControlScaleSigdig(
      do.call(rxode2::rxControl, rxControl),
      sigdig,
      skip = names(rxControl),
      tighten = 3
    )
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'", call = FALSE)
  }
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower = 1, finite = TRUE, any.missing = TRUE, len = 1)
    if (is.null(sigdigTable)) {
      sigdigTable <- round(sigdig)
    }
  }
  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower = 1, len = 1, any.missing = FALSE)

  .iterPrintControl <- .absorbIterPrintControl(
    print = print,
    printNcol = printNcol,
    useColor = useColor,
    iterPrintControl = .xtra$iterPrintControl
  )
  scaleType <- .ctlIdx(scaleType, .scaleTypeIdx, match.arg(scaleType))

  normType <- .ctlIdx(normType, .normTypeIdx, match.arg(normType))
  checkmate::assertNumeric(scaleCmax, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(scaleCmin, lower = 0, any.missing = FALSE, len = 1)
  if (!is.null(scaleC)) {
    checkmate::assertNumeric(scaleC, lower = 0, any.missing = FALSE)
  }
  checkmate::assertNumeric(scaleTo, len = 1, lower = 0, any.missing = FALSE)
  checkmate::assertNumeric(gradTo, len = 1, lower = 0, any.missing = FALSE)

  .ret <- list(
    algorithm = algorithm,
    maxiter = maxiter,
    tol = tol,
    trace = trace, # nolint
    minFactor = minFactor,
    printEval = printEval,
    warnOnly = warnOnly,
    scaleOffset = scaleOffset,
    nDcentral = nDcentral,
    solveType = solveType,
    linCmtSensCarry = match.arg(linCmtSensCarry),
    stickyRecalcN = stickyRecalcN,
    maxOdeRecalc = maxOdeRecalc,
    odeRecalcFactor = odeRecalcFactor,
    indTolRelax = indTolRelax,
    eventType = eventType,
    shiErr = shiErr,
    shi21maxFD = shi21maxFD,
    ftol = ftol,
    ptol = ptol,
    gtol = gtol,
    diag = diag,
    epsfcn = epsfcn,
    factor = factor,
    maxfev = maxfev,
    nprint = nprint,
    iterPrintControl = .iterPrintControl,
    scaleType = scaleType,
    normType = normType,
    scaleCmax = scaleCmax,
    scaleCmin = scaleCmin,
    scaleC = scaleC,
    scaleTo = scaleTo,
    gradTo = gradTo,
    optExpression = optExpression,
    literalFix = literalFix,
    sumProd = sumProd,
    rxControl = rxControl,
    returnNls = returnNls,
    addProp = match.arg(addProp),
    eventSens = match.arg(eventSens),
    calcTables = calcTables,
    compress = compress,
    ci = ci,
    sigdig = sigdig,
    sigdigTable = sigdigTable,
    genRxControl = .genRxControl,
    boundedTransform = boundedTransform
  )
  class(.ret) <- "nlsControl"
  .ret
}

#' @export
rxUiDeparse.nlsControl <- function(object, var) .deparseControl(object, var, nlsControl())

#' Get the nls family control
#'
#' @param env nlme optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlsFamilyControl <- function(env, ...) {
  .nlmFamilyControlGeneric(env, nlmixr2est::nlsControl, "nlsControl")
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.nlsControl <- function(control, env) assign("nlsControl", control, envir = env)


#' @rdname nmObjGetControl
#' @export
nmObjGetControl.nls <- function(x, ...) .nmObjGetControlByClass(x, "nlsControl")

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.nls <- function(control) .getValidCtl(control, "nlsControl", "nls")


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
  do.call(
    rxode2::rxSolve,
    c(
      list(
        object = nlmixr2global$nlsEnv$model,
        params = nlmixr2global$nlsEnv$parFun(...),
        events = nlmixr2global$nlsEnv$data
      ),
      nlmixr2global$nlsEnv$rxControl
    )
  )$rx_pred_
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
  .line <- rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1, .errNum, rxPredLlik = .getRxPredLlikOption())
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
    # First residual error is divided out (nls estimates it directly); add+prop/add+pow unsupported.
    .errType <- as.character(pred1$errType)
    if (.errType == "add") {
      # In these cases you are simply dividing out the additive error
      # Simply force this to be one.
      .lineExtra <- c(.lineExtra, list(quote(rx_r_ ~ 1)))
    } else if (.errType == "prop") {
      #   rx_r_ ~ (rx_pred_f_ * prop.sd)^2
      .f <- pred1$f
      .type <- as.character(pred1$errTypeF)
      .lineExtra <- c(
        .lineExtra,
        list(switch(
          .type,
          untransformed = quote(rx_r_ ~ (rx_pred_f_)^2),
          transformed = quote(rx_r_ ~ (rx_pred_)^2),
          f = bquote(rx_r_ ~ (.(str2lang(.f)))^2),
          none = quote(rx_r_ ~ (rx_pred_f_)^2)
        ))
      )
    } else if (.errType == "pow") {
      .cnd <- pred1$cond
      if (!is.na(pred1$c)) {
        .p2 <- str2lang(pred1$c)
      } else {
        .w <- which(env$iniDf$err %in% c("pow2", "powF2", "powT2") & env$iniDf$condition == .cnd)
        if (length(.w) == 1L) {
          .p2 <- str2lang(env$iniDf$name[.w])
        } else {
          stop("cannot find exponent of power expression", call. = FALSE)
        }
      }
      .f <- pred1$f
      .type <- as.character(pred1$errTypeF)
      .lineExtra <- c(
        .lineExtra,
        list(switch(
          .type,
          untransformed = bquote(rx_r_ ~ (rx_pred_f_)^(2 * .(.p2))),
          transformed = bquote(rx_r_ ~ (rx_pred_)^(2 * .(.p2))),
          f = bquote(rx_r_ ~ (.(str2lang(.f)))^(2 * .(.p2))),
          none = quote(rx_r_ ~ (rx_pred_f_)^(2 * .(.p2)))
        ))
      )
    }
  }
  c(.line, .lineExtra)
}

#' @rdname rxGetDistributionNlsLines
#' @export
rxGetDistributionNlsLines.default <- function(line) {
  stop("only normally related endoints can be used with 'nls', try 'nlm'", call. = FALSE)
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
  # .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  .env <- new.env(parent = emptyenv())
  .env$i <- 0
  .w <- which(!(rxui$iniDf$err %in% c("add", "prop", "pow")))
  lapply(.w, function(i) {
    if (rxui$iniDf$fix[i]) {
      return(eval(parse(text = paste0("quote(", .iniDf$name[i], " <- ", rxui$iniDf$est[i], ")"))))
    }
    if (rxui$iniDf$err[i] %in% c("add", "prop", "pow")) {
      return(NULL)
    }
    .env$i <- .env$i + 1
    eval(parse(text = paste0("quote(", .iniDf$name[i], " <- THETA[", .env$i, "])")))
  })
}

#' @export
rxUiGet.nlsModel0 <- function(x, ...) {
  .f <- x[[1]]
  .nlmFamilyModel0(
    .f,
    rxGetDistributionNlsLines(.f),
    .uiGetNlsTheta(.f),
    quote(rx_pred_ <- (rx_dv_ - rx_pred_) / sqrt(rx_r_))
  )
}
attr(rxUiGet.nlsModel0, "rstudio") <- quote(rxModelVars({}))

# The nls models are built by the nlm-family build stack in R/nlm.R;
# .nlmFamilySpec() lists what differs between the nlm and nls builds.

#' @export
rxUiGet.loadPruneNls <- function(x, ...) {
  .loadSymengine(.nlmFamilyPrune(x, "nls"), promoteLinSens = FALSE)
}
attr(rxUiGet.loadPruneNls, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsRxModel <- function(x, ...) {
  .nlmFamilyRxModel(x, "nls", ...)
}

#' @export
rxUiGet.loadPruneNlsSens <- function(x, ...) {
  .loadSymengine(.nlmFamilyPrune(x, "nls"), promoteLinSens = TRUE)
}
attr(rxUiGet.loadPruneNlsSens, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsThetaS <- function(x, ...) {
  .nlmFamilyThetaS(x, "nls")
}
attr(rxUiGet.nlsThetaS, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsHdTheta <- function(x, ...) {
  .nlmFamilyHdTheta(x, "nls")
}
attr(rxUiGet.nlsHdTheta, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsEnv <- function(x, ...) {
  .nlmFamilyEnv(x, "nls", ...)
}
attr(rxUiGet.nlsEnv, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlsSensModel <- function(x, ...) {
  .nlmFamilySensModel(x, "nls", ...)
}


#' @export
rxUiGet.nlsParStart <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  setNames(
    lapply(.w, function(i) {
      .ui$iniDf$est[i]
    }),
    .ui$iniDf$name[.w]
  )
}

#' @export
rxUiGet.nlsParStartTheta <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  setNames(
    vapply(
      .w,
      function(i) {
        .ui$iniDf$est[i]
      },
      double(1),
      USE.NAMES = FALSE
    ),
    paste0("THETA[", seq_along(.ui$iniDf$name[.w]), "]")
  )
}
attr(rxUiGet.nlsParStartTheta, "rstudio") <- c(`THETA[1]` = 0.1)

#' @export
rxUiGet.nlsParams <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  paste0("params(", paste(c(paste0("THETA[", seq_along(.ui$iniDf$name[.w]), "]"), "DV"), collapse = ", "), ")")
}
attr(rxUiGet.nlsParams, "rstudio") <- "params(THETA[1], DV)"

#' @export
rxUiGet.nlsParLower <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  setNames(
    vapply(
      .w,
      function(i) {
        .ui$iniDf$lower[i]
      },
      double(1),
      USE.NAMES = FALSE
    ),
    .ui$iniDf$name[.w]
  )
}
attr(rxUiGet.nlsParLower, "rstudio") <- c(`ka` = 0.01)

#' @export
rxUiGet.nlsParUpper <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow")))
  setNames(
    vapply(
      .w,
      function(i) {
        .ui$iniDf$upper[i]
      },
      double(1),
      USE.NAMES = FALSE
    ),
    .ui$iniDf$name[.w]
  )
}
attr(rxUiGet.nlsParUpper, "rstudio") <- c(`ka` = 1000)

#' @export
rxUiGet.nlsParNameFun <- function(x, ...) {
  .iniDf <- x[[1]]$iniDf
  # every THETA, with the residual-error and fixed ones at their estimates
  .values <- .iniDf$name
  .w <- .iniDf$err %in% c("add", "prop", "pow") | .iniDf$fix
  .values[.w] <- paste(.iniDf$est[.w])
  .nlmFamilyParNameFun(.nlsFormulaArgs(x)[-1], .values)
}
attr(rxUiGet.nlsParNameFun, "rstudio") <- function() {}

.nlsFormulaArgs <- function(x) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .args <- vapply(
    seq_along(.iniDf$ntheta),
    function(t) {
      if (.iniDf$err[t] %in% c("add", "prop", "pow")) {
        ""
      } else if (.iniDf$fix[t]) {
        ""
      } else {
        .iniDf$name[t]
      }
    },
    character(1),
    USE.NAMES = FALSE
  )
  c("DV", .args[.args != ""])
}

#' @export
rxUiGet.nlsFormula <- function(x, ..., grad = FALSE) {
  .args <- .nlsFormulaArgs(x)
  str2lang(paste0(
    "~nlmixr2est::.nlmixrNlsFunValGrad(",
    paste(.args, collapse = ", "),
    ")"
  ))
}
attr(rxUiGet.nlsFormula, "rstudio") <- quote(~ nlmixr2est::.nlmixrNlsFunValGrad(DV, ka, V, CL))

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
  .cens <- which(tolower(names(dataSav)) == "cens")
  .lim <- which(tolower(names(dataSav)) == "limit")
  .evid <- which(tolower(names(dataSav)) == "evid")
  if (length(.evid) == 1) {
    .wObs <- which(dataSav[[.evid]] == 0 | dataSav[[.evid]] == 2)
  } else {
    .wObs <- seq_len(nrow(dataSav))
  }
  if (length(.cens) == 1L) {
    if (!all(dataSav[.wObs, .cens] == 0)) {
      stop("'nls' does not work with censored data", call. = FALSE)
    }
  }
  if (length(.lim) == 1L) {
    if (any(is.finite(dataSav[.wObs, .lim]))) {
      stop("'nls' does not work with limit data", call. = FALSE)
    }
  }
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl, lower = ui$nlsParLower, upper = ui$nlsParUpper)
  .env$par.ini.list <- setNames(as.list(.env$par.ini), names(ui$nlsParStart))

  if (.ctl$algorithm == "LM") {
    .nls.control <- minpack.lm::nls.lm.control(
      ftol = .ctl$ftol,
      ptol = .ctl$ptol,
      gtol = .ctl$gtol,
      diag = .ctl$diag,
      epsfcn = .ctl$epsfcn,
      factor = .ctl$factor,
      maxfev = .ctl$maxfev,
      maxiter = .ctl$maxiter,
      nprint = .ctl$nprint
    )
    class(.ctl) <- NULL
    if (.ctl$solveType == 11L) {
      .ret <- bquote(minpack.lm::nls.lm(
        par = .(.env$par.ini),
        lower = .(.env$lower),
        upper = .(.env$upper),
        fn = nlmixr2est::.nlmixrNlsFunVal,
        control = .(.nls.control)
      ))
    } else {
      .ret <- bquote(minpack.lm::nls.lm(
        par = .(.env$par.ini),
        lower = .(.env$lower),
        upper = .(.env$upper),
        fn = nlmixr2est::.nlmixrNlsFunVal,
        jac = nlmixr2est::.nlmixrNlsFunGrad,
        control = .(.nls.control)
      ))
    }
    .ret <- eval(.ret)
    .ret <- .nlmFinalizeList(.env, .ret, par = "par", printLine = TRUE, hessianCov = TRUE)
    .ret$sd <- sd(.ret$fvec)
    .ret$logLik <- sum(stats::dnorm(.ret$fvec, log = TRUE))
  } else {
    nlmixr2global$nlsEnv$dataNls <- dataSav[dataSav$EVID == 0, ]
    .nls.control <- stats::nls.control(
      maxiter = .ctl$maxiter,
      tol = .ctl$tol,
      minFactor = .ctl$minFactor,
      printEval = .ctl$printEval,
      warnOnly = .ctl$warnOnly,
      scaleOffset = .ctl$scaleOffset,
      nDcentral = .ctl$nDcentral
    )
    class(.ctl) <- NULL
    .ret <- bquote(stats::nls(
      formula = .(ui$nlsFormula),
      data = nlmixr2est::.nlmixrNlsData(),
      start = .(.env$par.ini.list),
      control = .(.nls.control),
      algorithm = .(.ctl$algorithm),
      trace = .(.ctl$trace), # nolint
      model = FALSE,
      lower = .(.env$lower),
      upper = .(.env$upper)
    ))
    .ret <- eval(.ret)
    .ret <- .nlmFinalizeList(.env, .ret, printLine = TRUE, hessianCov = FALSE)
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
  setNames(
    vapply(
      seq_along(.iniDf$ntheta),
      function(t) {
        if (.iniDf$err[t] %in% c("add", "prop", "pow")) {
          .sd
        } else if (.iniDf$fix[t]) {
          .iniDf$est[t]
        } else {
          .theta0[.iniDf$name[t]]
        }
      },
      double(1),
      USE.NAMES = FALSE
    ),
    .iniDf$name
  )
}

.nlsControlToFoceiControl <- function(env, assign = TRUE) {
  .nlmFamilyControlToFoceiControl(env, "nlsControl", assign, literalFixRes = FALSE)
}

.nlsFamilyFit <- function(env, ...) {
  .nlmFamilyFitGeneric(
    env,
    "nls",
    .nlsFitModel,
    .nlsGetTheta,
    controlToFocei = .nlsControlToFoceiControl,
    returnFlag = "returnNls",
    # objective + cov + covMethod are set per-branch in postSetup (before
    # .nlmFamilyAdjustOutput, whose is.null guards then keep them)
    objective = NULL,
    message = function(.fit) {
      if (inherits(.fit, "nls.lm")) .fit$message else .fit$convInfo$stopMessage
    },
    extra = function(.control) {
      paste0(" with ", crayon::bold$yellow(.control$algorithm), " algorithm")
    },
    postSetup = function(.ret, .ui, .fit) {
      .ret$cov <- .ret$nls$cov
      if (inherits(.ret$nls, "nls.lm")) {
        .ret$covMethod <- paste0(.ret$nls$covMethod, " (LM)")
        .ret$objective <- -2 * .ret$nls$logLik
      } else {
        .ret$covMethod <- "nls"
        .ret$objective <- -2 * as.numeric(stats::logLik(.ret$nls))
      }
      .ret
    }
  )
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.nls <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiNoAutoregressive(.ui, " for the estimation routine 'nls', try 'focei'", .var.name = .ui$modelName)
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'nls', try 'focei'", .var.name = .ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'nls'", .var.name = .ui$modelName)
  rxode2::assertRxUiSingleEndpoint(.ui, " for the estimation routine 'nls'", .var.name = .ui$modelName)
  rxode2::assertRxUiEstimatedResiduals(.ui, " for the estimation routine 'nls'", .var.name = .ui$modelName)
  rxode2::warnRxBounded(.ui, " which are ignored in 'nls'", .var.name = .ui$modelName)
  # No add+prop or add+pow
  # Single endpoint
  .nlsFamilyControl(env, ...)
  on.exit(
    {
      if (exists("control", envir = .ui)) rm("control", envir = .ui)
    },
    add = TRUE
  )
  .nlsFamilyFit(env, ...)
}
attr(nlmixr2Est.nls, "covPresent") <- TRUE
attr(nlmixr2Est.nls, "unbounded") <- TRUE
