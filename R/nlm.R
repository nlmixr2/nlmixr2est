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
#' @param linCmtSensCarry use the linCmt() sensitivity carry for a theta on
#'   a covariate-driven linCmt() parameter ("auto", the default) or keep the
#'   standard gradient ("none"); see `foceiControl()`'s argument of the same
#'   name.
#' @param sensMethod Method used to compute the ODE parameter sensitivities.
#'   `"forward"` uses the classic variational (forward) sensitivity ODEs;
#'   `"default"` is the same thing.
#'
#' @return nlm control object
#' @export
#' @author Matthew L. Fidler
#' @examples
#' \donttest{
#' # A logit regression example with emax model
#'
#' dsn <- data.frame(i = 1:1000)
#' dsn$time <- exp(rnorm(1000))
#' dsn$DV <- rbinom(1000, 1, exp(-1 + dsn$time) / (1 + exp(-1 + dsn$time)))
#'
#' mod <- function() {
#'   ini({
#'     E0 <- 0.5
#'     Em <- 0.5
#'     E50 <- 2
#'     g <- fix(2)
#'   })
#'   model({
#'     v <- E0 + Em * time^g / (E50^g + time^g)
#'     ll(bin) ~ DV * v - log(1 + exp(v))
#'   })
#' }
#'
#' fit2 <- nlmixr(mod, dsn, est = "nlm")
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
nlmControl <- function(
  typsize = NULL,
  fscale = 1,
  print.level = 0,
  ndigit = NULL,
  gradtol = NULL,
  stepmax = NULL,
  steptol = NULL,
  iterlim = 10000,
  check.analyticals = FALSE,
  returnNlm = FALSE,
  solveType = c("hessian", "grad", "fun"),
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  indTolRelax = TRUE,
  eventType = c("central", "forward"),
  shiErr = (.Machine$double.eps)^(1 / 3),
  shi21maxFD = 20L,
  optimHessType = c("central", "forward"),
  hessErr = (.Machine$double.eps)^(1 / 3),
  shi21maxHess = 20L,
  censOption = c("gauss", "laplace"),
  eventSens = c("jump", "fd"),
  linCmtSensCarry = c("auto", "none"),
  sensMethod = c("default", "forward"),
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
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", "nlm", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 3,
  sigdigTable = NULL,
  boundedTransform = TRUE,
  ...
) {
  checkmate::assertNumeric(shiErr, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(hessErr, lower = 0, any.missing = FALSE, len = 1)

  checkmate::assertIntegerish(shi21maxFD, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(shi21maxHess, lower = 1, any.missing = FALSE, len = 1)

  checkmate::assertLogical(optExpression, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFix, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFixRes, len = 1, any.missing = FALSE)
  checkmate::assertLogical(sumProd, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(stepmax, lower = 0, len = 1, null.ok = TRUE, any.missing = FALSE)
  checkmate::assertIntegerish(print.level, lower = 0, upper = 2, any.missing = FALSE)
  checkmate::assertNumeric(ndigit, lower = 0, len = 1, any.missing = FALSE, null.ok = TRUE)
  # nlm gradtol/steptol keyed to the shared sigdig target (10^-sigdig), matching the
  # ODE rtol so nlm converges to the precision the solve supports; a user value wins
  if (is.null(gradtol)) {
    gradtol <- if (!is.null(sigdig)) .sigdigOptTol(sigdig) else 1e-6
  }
  if (is.null(steptol)) {
    steptol <- if (!is.null(sigdig)) .sigdigOptTol(sigdig) else 1e-6
  }
  checkmate::assertNumeric(gradtol, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(steptol, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(iterlim, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assertLogical(check.analyticals, len = 1, any.missing = FALSE)
  checkmate::assertLogical(returnNlm, len = 1, any.missing = FALSE)
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

  checkmate::assertIntegerish(stickyRecalcN, any.missing = FALSE, lower = 0, len = 1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(odeRecalcFactor, len = 1, lower = 1, any.missing = FALSE)
  checkmate::assertLogical(indTolRelax, any.missing = FALSE, len = 1)

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(ndigit)) {
    ndigit <- sigdig
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- .rxControlScaleSigdig(rxode2::rxControl(sigdig = sigdig), sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol = 1e-4, rtol = 1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {} else if (is.list(rxControl)) {
    rxControl <- .rxControlScaleSigdig(do.call(rxode2::rxControl, rxControl), sigdig, skip = names(rxControl))
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

  solveType <- .nlmCtlCode(solveType, c("hessian" = 3L, "grad" = 2L, "fun" = 1L), "solveType")
  if (missing(covMethod) && any(solveType == 2:3)) {
    covMethod <- "nlm"
  } else {
    covMethod <- match.arg(covMethod)
  }

  eventType <- .nlmCtlCode(eventType, c("central" = 2L, "forward" = 1L), "eventType")
  optimHessType <- .nlmCtlCode(optimHessType, c("central" = 2L, "forward" = 1L), "optimHessType")
  # censOption: FOCEI-family censored (M2/M3/M4) 2nd-derivative treatment -- "gauss" (historic
  # Gauss-Newton, default) or "laplace" (exact).  Accepted for a uniform interface but INERT for
  # NLM (its finite-difference Hessian already reflects censoring exactly); kept for alignment.
  if (checkmate::testIntegerish(censOption, len = 1, lower = 0, upper = 1, any.missing = FALSE)) {
    censOption <- as.integer(censOption)
  } else {
    censOption <- setNames(c("gauss" = 0L, "laplace" = 1L)[match.arg(censOption)], NULL)
  }

  ## eventSens: "jump" routes dosing-parameter (alag/F/rate/dur) sensitivities
  ## through rxode2's analytic event jumps; "fd" uses the legacy path that misses them.
  eventSens <- match.arg(eventSens)

  ## sensMethod: forward (variational) ODE parameter sensitivities.  Retained as
  ## a control so an explicit sensMethod="forward" keeps working.
  sensMethod <- match.arg(sensMethod)

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
    covMethod = covMethod,
    typsize = typsize,
    fscale = fscale,
    print.level = print.level,
    ndigit = ndigit,
    gradtol = gradtol,
    stepmax = stepmax,
    steptol = steptol,
    iterlim = iterlim,
    check.analyticals = check.analyticals,
    optExpression = optExpression,
    literalFix = literalFix,
    literalFixRes = literalFixRes,
    sumProd = sumProd,
    rxControl = rxControl,
    returnNlm = returnNlm,
    stickyRecalcN = as.integer(stickyRecalcN),
    maxOdeRecalc = as.integer(maxOdeRecalc),
    odeRecalcFactor = odeRecalcFactor,
    indTolRelax = indTolRelax,
    eventType = eventType,
    shiErr = shiErr,
    shi21maxFD = as.integer(shi21maxFD),
    optimHessType = optimHessType,
    hessErr = hessErr,
    shi21maxHess = as.integer(shi21maxHess),
    censOption = censOption,
    eventSens = eventSens,
    sensMethod = sensMethod,
    iterPrintControl = .iterPrintControl,
    scaleType = scaleType,
    normType = normType,
    scaleCmax = scaleCmax,
    scaleCmin = scaleCmin,
    scaleC = scaleC,
    scaleTo = scaleTo,
    gradTo = gradTo,
    addProp = match.arg(addProp),
    calcTables = calcTables,
    compress = compress,
    solveType = solveType,
    linCmtSensCarry = match.arg(linCmtSensCarry),
    ci = ci,
    sigdig = sigdig,
    sigdigTable = sigdigTable,
    genRxControl = .genRxControl,
    boundedTransform = boundedTransform
  )
  class(.ret) <- "nlmControl"
  .ret
}

#' @export
rxUiDeparse.nlmControl <- function(object, var) .deparseControl(object, var, nlmControl())


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
nmObjHandleControlObject.nlmControl <- function(control, env) assign("nlmControl", control, envir = env)

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.nlm <- function(x, ...) .nmObjGetControlByClass(x, "nlmControl")

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.nlm <- function(control) .getValidCtl(control, "nlmControl", "nlm")

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
  .env <- new.env(parent = emptyenv())
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

#' @export
rxUiGet.nlmModel0 <- function(x, ...) {
  .ui <- rxode2::rxUiDecompress(x[[1]])
  # on.exit() registered BEFORE mutating either flag: an interrupt landing
  # between setting a flag and registering its reset would otherwise leak
  # it TRUE for the rest of the R session, changing how a later fit's model
  # is generated.
  on.exit(nlmixr2global$rxCensNuFix <- FALSE, add = TRUE)
  on.exit(nlmixr2global$rxPredLlik <- FALSE, add = TRUE)
  nlmixr2global$rxPredLlik <- TRUE
  # expose the censoring inputs (a real rx_r_, plus rx_nu_) for the
  # llik-forced endpoint (.fixCensRNuLine, R/focei.R)
  nlmixr2global$rxCensNuFix <- TRUE
  .predDf <- .ui$predDf
  .save <- .predDf
  .predDf[.predDf$distribution == "norm", "distribution"] <- "dnorm"
  assign(".predDfFocei", .predDf, envir = .ui)
  on.exit(assign("predDf", .save, envir = .ui), add = TRUE)
  # rx_r_/rx_nu_ for the llik-forced norm/dnorm/t/cauchy path are already
  # fixed at the source (.fixCensRNuLine, R/focei.R) -- see #979.
  .errLines <- rxGetDistributionFoceiLines(.ui)
  .nlmFamilyModel0(.ui, .errLines, .uiGetThetaDropFixed(.ui), quote(rx_pred_ <- -rx_pred_))
}
attr(rxUiGet.nlmModel0, "rstudio") <- quote(rxModelVar({}))

#' Assemble the `rxModelVars({...})` call an nlm-family model is built from
#'
#' @param ui rxode2 UI object
#' @param errLines,prefixLines error and theta lines for
#'   `rxode2::rxCombineErrorLines()`
#' @param lastLine last line, which turns `rx_pred_` into what the method
#'   minimizes: the -LL for nlm (`rxUiGet.nlmModel0`), the weighted residual
#'   for nls (`rxUiGet.nlsModel0`)
#' @return `rxModelVars({...})` call
#' @noRd
.nlmFamilyModel0 <- function(ui, errLines, prefixLines, lastLine) {
  .ret <- rxode2::rxCombineErrorLines(
    ui,
    errLines = errLines,
    prefixLines = prefixLines,
    paramsLine = NA,
    modelVars = TRUE,
    cmtLines = FALSE,
    dvidLine = FALSE
  )[[2]]
  as.call(list(quote(`rxModelVars`), as.call(c(as.list(.ret), list(lastLine)))))
}

#' What differs between the nlm and nls model builds
#'
#' The nlm-family methods (nlm, nlminb, optim, bobyqa, ...) minimize the
#' population -LL of `rxUiGet.nlmModel0`; nls minimizes the weighted residual
#' of `rxUiGet.nlsModel0`.  Both models go through the one build stack below
#' (`.nlmFamilyPrune()` to `.nlmFamilySensModel()`).  Apart from the theta
#' numbering of the `params()` line and the names in messages, the builds
#' differ only in:
#'
#' - `matExpForcing` (see `.sensEtaOrTheta()`): nlm flattens a matExp() model
#'   with an `indLin()` forcing term to ODEs.  nls keeps the native
#'   matrix-exponential sensitivities; its residual Jacobian from them matches
#'   finite differences.
#' - `censFR`: nlm emits the `rx_pred_f_`/`rx_r_`/`rx_nu_` outputs that the
#'   censoring likelihood in src/nlm.cpp reads.  nls refuses censored and
#'   limit data (`.nlsFitModel()`), and its `rx_pred_` is a residual, not the
#'   -LL that the censoring likelihood replaces.
#' - `lhs`: nlm copies the model lhs into its gradient and pred-only models.
#'   They define the `k_*` rate constants of a flattened matExp() model and the
#'   variables referenced by `lag()`.  nls needs neither: it never flattens a
#'   matExp() model, and a variable referenced by `lag()` has no symbolic
#'   sensitivity, so it cannot enter the residual Jacobian.  (The
#'   objective-only model of both defines such a variable, see
#'   `.nlmFamilyRxModel()`.)
#'
#' @param type `"nlm"` or `"nls"`
#' @return list of the settings for `type`
#' @noRd
.nlmFamilySpec <- function(type) {
  switch(
    type,
    nlm = list(
      model = "population log-likelihood model",
      llik = "nlm llik",
      role = "Nlm",
      loadPrune = rxUiGet.loadPruneNlm,
      params = rxUiGet.nlmParams,
      matExpForcing = FALSE,
      censFR = TRUE,
      lhs = TRUE
    ),
    nls = list(
      model = "nls model",
      llik = "nls",
      role = "Nls",
      loadPrune = rxUiGet.loadPruneNls,
      params = rxUiGet.nlsParams,
      matExpForcing = TRUE,
      censFR = FALSE,
      lhs = FALSE
    )
  )
}

#' Prune the `if`/`else` branches of an nlm-family model
#'
#' @param x rxode2 UI object
#' @param type `"nlm"` or `"nls"` (see `.nlmFamilySpec()`)
#' @return String for loading into symengine
#' @author Matthew L. Fidler
#' @noRd
.nlmFamilyPrune <- function(x, type) {
  .x <- x[[1]]
  .x <- switch(type, nlm = .x$nlmModel0, nls = .x$nlsModel0)[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .malert(paste0("pruning branches ({.code if}/{.code else}) of ", .nlmFamilySpec(type)$model, "..."))
  .ret <- rxode2::.rxPrune(.x, envir = .env, strAssign = rxode2::rxModelVars(x[[1]])$strAssign)
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
  .p <- .nlmFamilyPrune(x, "nlm")
  .loadSymengine(.p, promoteLinSens = FALSE)
}
attr(rxUiGet.loadPruneNlm, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlmParams <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .w <- which(!.iniDf$fix)
  .env <- new.env(parent = emptyenv())
  .env$t <- 0
  ## Declare the model covariates (ui$allCovs) explicitly after DV.  Referenced
  ## covariates land here anyway (auto-detected after the thetas), so this only
  ## pins the order -- but it also keeps covariates whose only reference is dropped
  ## by log-likelihood pruning (e.g. a plugin's externally-loaded parameter block,
  ## read by a compiled function at a fixed par_ptr index rather than by name) in
  ## the solve parameter layout so they retain a stable par_ptr slot.
  .covs <- .ui$allCovs
  if (is.null(.covs)) {
    .covs <- character(0)
  }
  paste0(
    "params(",
    paste(
      c(
        vapply(
          .w,
          function(i) {
            .env$t <- .env$t + 1
            paste0("THETA[", .env$t, "]")
          },
          character(1),
          USE.NAMES = FALSE
        ),
        "DV",
        .covs
      ),
      collapse = ","
    ),
    ")"
  )
}
attr(rxUiGet.nlmParams, "rstudio") <- "params()"

#' Extract rx_pred_f_, rx_r_ and rx_nu_ model lines from symengine environment
#'
#' Like `rx_pred_f_`/`rx_r_`, a plain `rx_nu_ ~ nu` (or `=`) line does not
#' survive `rxOptExpr`/symengine's own dead-code elimination -- nothing else
#' in the compiled model reads it back, so it is pruned along with any other
#' truly-unused intermediate before `rxCombineErrorLines`'s caller ever gets
#' to see it. It has to be pulled straight out of the symengine environment
#' (which still has the full unpruned variable set) and re-spliced into the
#' final text, exactly like `rx_pred_f_`/`rx_r_` already are (#979).
#'
#' @param .s symengine environment
#' @return named list with `f_line`, `r_line` and `nu_line` character strings
#' @noRd
.nlmGetFRLines <- function(.s) {
  .f_line <- ""
  .r_line <- ""
  .nu_line <- ""
  if (exists("rx_pred_f_", envir = .s, inherits = FALSE)) {
    .f <- get("rx_pred_f_", envir = .s)
    .f_line <- paste0("rx_pred_f_=", rxode2::rxFromSE(.f))
  }
  if (exists("rx_r_", envir = .s, inherits = FALSE)) {
    .r <- get("rx_r_", envir = .s)
    .r_line <- paste0("rx_r_=", rxode2::rxFromSE(.r))
  }
  if (exists("rx_nu_", envir = .s, inherits = FALSE)) {
    .nu <- get("rx_nu_", envir = .s)
    .nu_line <- paste0("rx_nu_=", rxode2::rxFromSE(.nu))
  }
  list(f_line = .f_line, r_line = .r_line, nu_line = .nu_line)
}

#' Flag the THETAs that dosing parameters (alag/F/rate/dur) depend on
#'
#' @param s symengine environment
#' @param flag when `FALSE`, return all zeros
#' @return 0/1 integer vector, one element per THETA
#' @noRd
.nlmFamilyEventTheta <- function(s, flag = TRUE) {
  if (exists("..maxTheta", s)) {
    .eventTheta <- rep(0L, s$..maxTheta)
  } else {
    .eventTheta <- integer(0)
  }
  if (!flag) {
    return(.eventTheta)
  }
  for (.v in s$..eventVars) {
    .vars <- as.character(get(.v, envir = s))
    .vars <- rxode2::rxGetModel(paste0("rx_lhs=", rxode2::rxFromSE(.vars)))$params
    for (.v2 in .vars) {
      .reg <- rex::rex(start, "THETA[", capture(any_numbers), "]", end)
      if (regexpr(.reg, .v2) != -1) {
        .num <- as.numeric(sub(.reg, "\\1", .v2))
        .eventTheta[.num] <- 1L
      }
    }
  }
  .eventTheta
}

#' @export
rxUiGet.nlmRxModel <- function(x, ...) {
  .nlmFamilyRxModel(x, "nlm", ...)
}

#' Objective-only (`solveType = "fun"`) model of an nlm-family method
#'
#' @param x rxode2 UI object
#' @param type `"nlm"` or `"nls"` (see `.nlmFamilySpec()`)
#' @param ... passed to the `rxUiGet` methods
#' @return list with the `predOnly` model and the `eventTheta` flags
#' @noRd
.nlmFamilyRxModel <- function(x, type, ...) {
  .spec <- .nlmFamilySpec(type)
  .s <- .spec$loadPrune(x, ...)
  # For matExp() models materialize the implied d/dt() from the k_from_to rate
  # constants.  When this fires we must also emit the model LHS (which defines
  # the k_from_to constants and other assignments) ahead of the d/dt() lines so
  # the derivative expressions can resolve them.
  .isMatExp <- isTRUE(.rxInjectMatExpDdt(.s))
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .ddt <- .s$..ddt
  if (is.null(.ddt)) {
    .ddt <- ""
  }
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
  # rx_pred_f_/rx_r_/rx_nu_ outputs for censoring support
  .fr <- if (.spec$censFR) .nlmGetFRLines(.s) else list()
  .ret <- paste(
    c(
      .lhs,
      .ddt,
      .lagDefs,
      ## DDE non-constant delay() pre-history (base past(state,tau)<-expr)
      rxode2::.rxPastBaseLinesFromEnv(.s),
      .prd,
      .fr$f_line,
      .fr$r_line,
      .fr$nu_line,
      ""
    ),
    collapse = "\n"
  )
  .eventTheta <- .nlmFamilyEventTheta(.s)
  .s$.eventTheta <- .eventTheta
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  if (.sumProd) {
    .malert(paste0("stabilizing round off errors in ", .spec$model, "..."))
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, .spec$model, parallel = .optExprCores(x[[1]]))
    .msuccess("done")
  }
  .cmt <- rxUiGet.foceiCmtPreModel(x, ...)
  # mtime() lines are re-emitted here (#919); see .mtimeLinesStr()
  .cmt <- .addPreModelLines(.cmt, rxUiGet.interpLinesStr(x, ...), .mtimeLinesStr(.s))
  ## no splitBolus() here -- this model solves the pre-split events, so
  ## declaring it would split the doses twice (see .foceiPreProcessData())
  list(
    predOnly = .nlmixr2estRxode2(
      paste(c(.spec$params(x, ...), .cmt, .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse = "\n"),
      paste0("rx", .spec$role, "PredOnly")
    ),
    eventTheta = .eventTheta
  )
}

#' @export
rxUiGet.loadPruneNlmSens <- function(x, ...) {
  .loadSymengine(.nlmFamilyPrune(x, "nlm"), promoteLinSens = TRUE)
}
attr(rxUiGet.loadPruneNlmSens, "rstudio") <- emptyenv()

#' @export
rxUiGet.nlmThetaS <- function(x, ...) {
  .nlmFamilyThetaS(x, "nlm")
}
attr(rxUiGet.nlmThetaS, "rstudio") <- emptyenv()

#' Load an nlm-family model with its theta sensitivities into symengine
#'
#' @param x rxode2 UI object
#' @param type `"nlm"` or `"nls"` (see `.nlmFamilySpec()`)
#' @return symengine environment from `.sensEtaOrTheta()`
#' @noRd
.nlmFamilyThetaS <- function(x, type) {
  .s <- .loadSymengine(.nlmFamilyPrune(x, type), promoteLinSens = TRUE)
  .sensEtaOrTheta(.s, theta = TRUE, rxui = x[[1]], matExpForcing = .nlmFamilySpec(type)$matExpForcing)
}

#' @export
rxUiGet.nlmHdTheta <- function(x, ...) {
  .nlmFamilyHdTheta(x, "nlm")
}
attr(rxUiGet.nlmHdTheta, "rstudio") <- emptyenv()

#' Calculate the d(f)/d(theta) lines of an nlm-family model
#'
#' @param x rxode2 UI object
#' @param type `"nlm"` or `"nls"` (see `.nlmFamilySpec()`)
#' @return symengine environment with `..HdTheta` set
#' @noRd
.nlmFamilyHdTheta <- function(x, type) {
  .s <- .nlmFamilyThetaS(x, type)
  .stateVars <- rxode2stateOde(.s)
  .predMinusDv <- rxode2::rxGetControl(x[[1]], "predMinusDv", TRUE)
  .grd <- rxode2::rxExpandFEta_(
    .stateVars,
    .s$..maxTheta,
    ifelse(.predMinusDv, 1L, 2L),
    isTheta = TRUE
  )
  if (rxode2::.useUtf()) {
    .malert("calculate \u2202(f)/\u2202(\u03B8)")
  } else {
    .malert("calculate d(f)/d(theta)")
  }
  rxode2::rxProgress(dim(.grd)[1])
  on.exit({
    rxode2::rxProgressAbort()
  })
  .zero <- new.env(parent = emptyenv())
  .zero$any <- FALSE
  .zero$all <- TRUE
  # linCmt() sensitivity carry for a theta on a covariate-driven linCmt()
  # parameter (#1003): the naive line for a carry-eligible theta is replaced
  # wholesale; everything else is byte-identical (foceiLinCmtCarryTheta.R)
  .thetaVars <- paste0("THETA_", seq_len(.s$..maxTheta), "_")
  .carry <- .rxCarryThetaPairsForBuild(x, .s, .thetaVars)
  .ret <- apply(.grd, 1, .nlmFamilyHdThetaLine, .s = .s, .carry = .carry, .predMinusDv = .predMinusDv, .zero = .zero)
  if (.zero$all) {
    stop("none of the predictions depend on 'THETA'", call. = FALSE)
  }
  if (.zero$any) {
    warning("some of the predictions do not depend on 'THETA'", call. = FALSE)
  }
  .s$..HdTheta <- .ret
  .s$..linCmtCarryThetaPairs <- if (is.null(.carry)) NULL else .carry$pairs
  .s$..pred.minus.dv <- .predMinusDv
  rxode2::rxProgressStop()
  .s
}

#' One d(f)/d(theta) line of `.nlmFamilyHdTheta()`
#'
#' @param x row of the `rxode2::rxExpandFEta_()` table; its `calc`
#'   expression refers to the symengine environment as `.s`
#' @param .s symengine environment
#' @param .carry linCmt() sensitivity carry from `.rxCarryThetaPairsForBuild()`
#' @param .predMinusDv the `predMinusDv` control setting
#' @param .zero environment whose `any`/`all` record whether any/all of the
#'   derivatives are zero
#' @return the `dfe=expr` model line
#' @noRd
.nlmFamilyHdThetaLine <- function(x, .s, .carry, .predMinusDv, .zero) {
  .l <- x["calc"]
  .l <- eval(parse(text = .l))
  .ret <- paste0(x["dfe"], "=", rxode2::rxFromSE(.l))
  if (!is.null(.carry)) {
    .w <- which(.carry$pairs$eta == sub("^.*_BY_(THETA_[0-9]+)___$", "\\1_", x["dfe"]))
    if (length(.w) == 1L) {
      .ret <- .rxCarryThetaEmit(.carry$pairs, .w, .s, x["dfe"], .carry$fp, .predMinusDv)
    }
  }
  .zErr <- suppressWarnings(try(as.numeric(get(x["dfe"], .s)), silent = TRUE))
  if (identical(.zErr, 0)) {
    .zero$any <- TRUE
  } else if (.zero$all) {
    .zero$all <- FALSE
  }
  rxode2::rxTick()
  .ret
}

#' Finalize nlm-family rxode2 models based on symengine saved info
#'
#' @param .s Symengine/rxode2 object
#' @param interpLines covariate interpolation lines (`locf()`/`nocb()`/...) to
#'   emit; symengine drops them, so they have to be added back here
#' @param type `"nlm"` or `"nls"` (see `.nlmFamilySpec()`)
#' @return Nothing; sets the gradient model (`..nlmS` or `..nlsS`) and the
#'   pred-only model (`..pred.nolhs`) in `.s`
#' @author Matthew L Fidler
#' @noRd
.rxFinalizeNlm <- function(.s, sum.prod = FALSE, optExpression = TRUE, cores = 0L, interpLines = "", type = "nlm") {
  .spec <- .nlmFamilySpec(type)
  interpLines <- interpLines[interpLines != ""]
  # see focei.R's .rxFinalizeInner(): do not re-flatten a matExp-native ..ddt (#860)
  if (!isTRUE(.s$..matExpNative)) {
    .rxInjectMatExpDdt(.s)
  }
  if (isTRUE(.s$..matExpNative)) {
    # see focei.R's .rxFinalizeInner(): rxSumProdModel()/rxOptExpr() do not
    # support "indLin(state) <- expr" (Michaelis-Menten forcing)
    sum.prod <- FALSE
    optExpression <- FALSE
  }
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
  if (is.null(.ddt)) {
    .ddt <- character(0)
  }
  .lhs <- character(0)
  if (.spec$lhs) {
    .lhs <- .s$..lhs
    if (is.null(.lhs)) {
      .lhs <- character(0)
    }
    # matExp-native sensitivities (#860): see focei.R's .rxFinalizeInner()
    .lhs <- .rxDropMatExpNativeLhs(.lhs, .s)
  }
  .sens <- .s$..sens
  if (is.null(.sens)) {
    .sens <- character(0)
  }
  # rx_pred_f_/rx_r_/rx_nu_ outputs for censoring support
  .fr <- if (.spec$censFR) .nlmGetFRLines(.s) else list()
  .grad <- paste(
    c(
      .s$params,
      .s$..stateInfo["state"],
      interpLines,
      .lhs,
      .ddt,
      .sens,
      ## DDE non-constant delay() pre-history: base past(state,tau)<-expr + the
      ## per-sensitivity-compartment histories (analytic gradient/Jacobian).
      .s$..pastLines,
      .yj,
      .lambda,
      .hi,
      .low,
      .prd,
      .s$..HdTheta,
      .fr$f_line,
      .fr$r_line,
      .fr$nu_line,
      .s$..stateInfo["statef"],
      .s$..stateInfo["dvid"],
      ""
    ),
    collapse = "\n"
  )
  .lhs0 <- .s$..lhs0
  if (is.null(.lhs0)) {
    .lhs0 <- ""
  }
  .s$..pred.nolhs <- paste(
    c(
      .s$params,
      .s$..stateInfo["state"],
      interpLines,
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
      .fr$nu_line,
      .s$..stateInfo["statef"],
      .s$..stateInfo["dvid"],
      ""
    ),
    collapse = "\n"
  )
  if (sum.prod) {
    .malert(paste0("stabilizing round off errors in ", .spec$llik, " gradient problem..."))
    .grad <- rxode2::rxSumProdModel(.grad)
    .msuccess("done")
    .malert(paste0("stabilizing round off errors in ", .spec$llik, " pred-only problem..."))
    .s$..pred.nolhs <- rxode2::rxSumProdModel(.s$..pred.nolhs)
    .msuccess("done")
  }
  if (optExpression) {
    .grad <- rxode2::rxOptExpr(.grad, paste0(.spec$llik, " gradient"), parallel = cores)
    .s$..pred.nolhs <- rxode2::rxOptExpr(.s$..pred.nolhs, paste0(type, " pred-only"), parallel = cores)
  }
  # mtime() lines go in AFTER the optimization, which cannot parse them (#919)
  assign(paste0("..", type, "S"), .addMtimeLines(.grad, .s), envir = .s)
  .s$..pred.nolhs <- .addMtimeLines(.s$..pred.nolhs, .s)
}

#' @export
rxUiGet.nlmEnv <- function(x, ...) {
  .nlmFamilyEnv(x, "nlm", ...)
}
attr(rxUiGet.nlmEnv, "rstudio") <- emptyenv()

#' Build the symengine environment holding an nlm-family gradient model
#'
#' @param x rxode2 UI object
#' @param type `"nlm"` or `"nls"` (see `.nlmFamilySpec()`)
#' @param ... passed to the `rxUiGet` methods
#' @return symengine environment; `.rxFinalizeNlm()` describes the models in it
#' @noRd
.nlmFamilyEnv <- function(x, type, ...) {
  .s <- .nlmFamilyHdTheta(x, type)
  .s$params <- .nlmFamilySpec(type)$params(x, ...)
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizeNlm(
    .s,
    .sumProd,
    .optExpression,
    .optExprCores(x[[1]]),
    interpLines = rxUiGet.interpLinesStr(x, ...),
    type = type
  )
  .s$..outer <- NULL
  ## eventTheta flags dosing-parameter (alag/F/rate/dur) THETAs, whose gradient
  ## is then taken by finite differences.  Under eventSens="jump" none are
  ## flagged: rxode2 injects the jump sensitivities analytically.
  .s$.eventTheta <- .nlmFamilyEventTheta(.s, !identical(rxode2::rxGetControl(x[[1]], "eventSens", "jump"), "jump"))
  .s
}

#' @export
rxUiGet.nlmSensModel <- function(x, ...) {
  .nlmFamilySensModel(x, "nlm", ...)
}

#' Gradient and pred-only models of an nlm-family method
#'
#' @param x rxode2 UI object
#' @param type `"nlm"` or `"nls"` (see `.nlmFamilySpec()`)
#' @param ... passed to the `rxUiGet` methods
#' @return list with the `thetaGrad` and `predOnly` models and the
#'   `eventTheta` flags
#' @noRd
.nlmFamilySensModel <- function(x, type, ...) {
  .role <- .nlmFamilySpec(type)$role
  .s <- .nlmFamilyEnv(x, type, ...)
  ## "jump" attaches rxode2's analytic event (alag/F/rate/dur) sensitivities to
  ## the thetaGrad model; under "fd" the flagged eventTheta are finite-differenced.
  .eventSens <- rxode2::rxGetControl(x[[1]], "eventSens", "jump")
  list(
    thetaGrad = .nlmixr2estRxode2(
      get(paste0("..", type, "S"), envir = .s),
      paste0("rx", .role, "Grad"),
      eventSens = .eventSens
    ),
    predOnly = .nlmixr2estRxode2(.s$..pred.nolhs, paste0("rx", .role, "Pred")),
    eventTheta = .s$.eventTheta
  )
}

#' Build a function mapping parameter values to named `THETA[#]` values
#'
#' @param args arguments of the function
#' @param values `THETA[1]`, `THETA[2]`, ... values, as R code in terms of
#'   `args`
#' @return `function(<args>) {c('THETA[1]'=<values[1]>, ...)}`
#' @noRd
.nlmFamilyParNameFun <- function(args, values) {
  eval(str2lang(paste0(
    "function(",
    paste(args, collapse = ", "),
    ") {c(",
    paste(sprintf("'THETA[%d]'=%s", seq_along(values), values), collapse = ","),
    ")}"
  )))
}

#' @export
rxUiGet.nlmParNameFun <- function(x, ...) {
  .nlmFamilyParNameFun("p", sprintf("p[%d]", seq_len(sum(!x[[1]]$iniDf$fix))))
}
attr(rxUiGet.nlmParNameFun, "rstudio") <- function() {
  c(`THETA[1]` = 1, `THETA[2]` = 2, `THETA[3]` = 3)
}

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
nlmObjectiveSetup <- function(ui, data, control = NULL, gradient = FALSE, scale = c("control", "natural")) {
  scale <- match.arg(scale)
  ## assertRxUi accepts a model function as well as a ui; .copyUi (not
  ## rxUiDecompress) then isolates it, because decompressing an already-
  ## decompressed ui hands back the SAME environment and assigning $control below
  ## would permanently rebind the caller's model to this control.
  .ui <- rxode2::.copyUi(rxode2::assertRxUi(ui))
  if (is.null(control)) {
    control <- if (!is.null(.ui$control)) .ui$control else nlmControl()
  }
  .ui$control <- control
  .ctl <- .ui$control
  class(.ctl) <- NULL
  if (gradient) {
    ## the C API (#953) hands out value + analytic d/d(theta): load the
    ## sensitivity model with the gradient solve type
    .ctl$solveType <- 2L
  }
  if (identical(scale, "natural")) {
    ## identity scale (scaleTypeNone + normTypeConstant), so the
    ## evaluated theta IS the model's theta -- what a sampler needs
    .ctl$scaleType <- 5L
    .ctl$normType <- 6L
  }
  .ret <- new.env(parent = emptyenv())
  .foceiPreProcessData(data, .ret, .ui, .ctl$rxControl)
  .p <- setNames(.ui$nlmParIni, .ui$nlmParName)
  ## gradient=FALSE -- solveType 1 / nlmRxModel: the objective-only predOnly
  ## model (no thetaGrad).  The hook fires from nlmSolveFid during the
  ## objective solve; the caller gets the weight gradient from its own
  ## augmented-sensitivity solve, so no analytic theta gradient is needed.
  .mi <- if (gradient) .ui$nlmSensModel else .ui$nlmRxModel
  .env <- .nlmSetupEnv(.p, .ui, .ret$dataSav, .mi, .ctl)
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
    stop("'typsize' needs to match the number of estimated parameters (or equal 1)", call. = FALSE)
  }
  .stepmax <- .ctl$stepmax
  if (is.null(.stepmax)) {
    .stepmax <- max(1000 * sqrt(sum((.p / .typsize)^2)), 1000)
  }
  .hessian <- .ctl$covMethod == "nlm"
  if (.ctl$solveType == 1L) {
    .mi <- ui$nlmRxModel
  } else {
    .mi <- ui$nlmSensModel
  }
  ## Event ("jump") sensitivities are activated in .nlmSetupEnv and deactivated in
  ## .nlmFreeEnv; nothing extra needed here.
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl)
  on.exit({
    .nlmFreeEnv()
  })
  .ret <- eval(bquote(stats::nlm(
    f = .(.nlmixrNlmFunC),
    p = .(.env$par.ini),
    hessian = .(.hessian),
    typsize = .(.typsize),
    fscale = .(.ctl$fscale),
    print.level = .(.ctl$print.level),
    ndigit = .(.ctl$ndigit),
    gradtol = .(.ctl$gradtol),
    stepmax = .(.stepmax),
    steptol = .(.ctl$steptol),
    iterlim = .(.ctl$iterlim),
    check.analyticals = .(.ctl$check.analyticals)
  )))
  .nlmFinalizeList(.env, .ret, par = "estimate", printLine = TRUE, hessianCov = TRUE)
}
.nlmControlToFoceiControl <- function(env, assign = TRUE) {
  .nlmFamilyControlToFoceiControl(env, "nlmControl", assign)
}


.nlmFamilyFit <- function(env, ...) {
  .nlmFamilyFitGeneric(
    env,
    "nlm",
    .nlmFitModel,
    "estimate",
    objective = "minimum",
    controlToFocei = .nlmControlToFoceiControl,
    returnFlag = "returnNlm",
    emitFitWarnings = TRUE,
    message = function(.fit) {
      if (.fit$code == 1) {
        "relative gradient is close to zero, current iterate is probably solution"
      } else if (.fit$code == 2) {
        "successive iterates within tolerance, current iterate is probably solution"
      } else if (.fit$code == 3) {
        c(
          "last global step failed to locate a point lower than 'estimate'",
          "either 'estimate' is an approximate local minimum of the function or 'steptol' is too small"
        )
      } else if (.fit$code == 4) {
        "iteration limit exceeded"
      } else if (.fit$code == 5) {
        c(
          "maximum step size 'stepmax' exceeded five consecutive times",
          "either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or 'stepmax' is too small" # nolint: line_length_linter.
        )
      } else {
        ""
      }
    }
  )
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.nlm <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'nlm', try 'focei'", .var.name = .ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'nlm'", .var.name = .ui$modelName)
  rxode2::warnRxBounded(.ui, " which are ignored in 'nlm'", .var.name = .ui$modelName)
  .nlmFamilyControl(env, ...)
  on.exit(
    {
      if (exists("control", envir = .ui)) rm("control", envir = .ui)
    },
    add = TRUE
  )
  .nlmFamilyFit(env, ...)
}
attr(nlmixr2Est.nlm, "covPresent") <- TRUE
attr(nlmixr2Est.nlm, "unbounded") <- TRUE
