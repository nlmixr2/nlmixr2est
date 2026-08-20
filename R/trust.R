#' Control for the trust estimation method in nlmixr2
#'
#' `est="trust"` is a trust-region Newton method for the population theta
#' vector, backed by the `RcppTrust` package's thread-safe `trust_solve_c()`.
#' Unlike every other nlm-family method (`nlm`/`nlminb`/`bobyqa`/`newuoa`/
#' `uobyqa`/`n1qn1`/`lbfgsb3c`/`optim`), whose optimization loop lives in R
#' (one `.Call()` per iteration), `trust`'s entire loop runs inside a single
#' `.Call()` -- `RcppTrust` needs no R API, so there is no per-iteration R
#' round-trip. It also supplies a full analytic-or-Shi21-finite-difference
#' gradient AND a full finite-difference-of-the-gradient Hessian every
#' iteration (there is no analytic outer-theta Hessian in this package), so
#' each outer iteration costs roughly `ntheta` extra full population-gradient
#' solves on top of the gradient solve itself -- state this cost plainly, it
#' is the price of true Newton-trust behavior rather than a derivative-free
#' or quasi-Newton method. `trust` is unbounded, like `n1qn1`/`nlm`: box
#' constraints are handled upstream by `preProcessBoundedTransform.R`, not by
#' the optimizer itself.
#'
#' @inheritParams iterPrintParams
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @inheritParams nlmControl
#'
#' @param rinit Initial trust-region radius, in the SAME scaled-parameter
#'   space every nlm-family method (bobyqa included) optimizes in. `NULL`
#'   (default) derives it from the scaled starting vector using
#'   `minqa::bobyqa()`'s own default-`rhobeg` formula,
#'   `min(0.95, 0.2*max(abs(par.ini)))`, so swapping `est="bobyqa"` for
#'   `est="trust"` on the same model starts from a comparable trust-region
#'   size.
#' @param rmax Maximum trust-region radius, same scaled-parameter space.
#'   `NULL` (default) derives it as `8 * rinit` -- the same growth-ceiling
#'   multiplier already used for `foceiControl(innerOpt="trust")`'s per-
#'   subject eta trust region.
#' @param iterlim Maximum number of `trust_solve_c()` iterations.
#' @param fterm,mterm `trust_solve_c()`'s function-value and predicted-
#'   decrease convergence tolerances. `NULL` (default) uses
#'   `10^(-sigdig)`, matching `bobyqaControl()`'s `rhoend` and
#'   `foceiControl()`'s inner `epsilon`. `mterm` defaults to `fterm` when
#'   not given separately.
#' @param optimHessType Finite-difference type for the per-iteration outer
#'   Hessian: `1L` forward (default), `2L` central. Exposed explicitly
#'   (unlike bobyqa/n1qn1, where this is only ever set implicitly) because
#'   `trust` recomputes this Hessian every outer iteration, not once
#'   post-fit, making its accuracy/cost tradeoff far more consequential; a
#'   stiff or noisy model may benefit from `optimHessType=2`.
#' @param shi21maxHess Maximum Shi (2021) adaptive-step-size iterations for
#'   the per-iteration outer Hessian.
#' @param hessErr Target relative error for the per-iteration outer Hessian's
#'   Shi (2021) step-size search. `NULL` (default) uses
#'   `(.Machine$double.eps)^(1/3)`, the same fallback `.nlmSetupEnv()` itself
#'   applies for every nlm-family method.
#' @param returnTrust return the raw `nlmTrustFit()` output list instead of
#'   the nlmixr2 fit.
#' @param covMethod Method for calculating the covariance. `"r"` (the
#'   default) reuses the LAST outer iteration's already-computed Hessian
#'   (skipping `nlmixr2est`'s own post-fit finite-difference Hessian
#'   recompute, since `trust` already has one in hand); `""` skips the
#'   covariance step.
#' @return trust control structure
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
#' fit <- nlmixr(mod, dsn, est="trust")
#'
#' print(fit)
#'
#' # you can also get the raw trust output with fit$trust
#'
#' fit$trust
#' }
trustControl <- function(rinit=NULL, rmax=NULL, iterlim=1000L,
                         fterm=NULL, mterm=NULL,
                         optimHessType=1L, shi21maxHess=20L, hessErr=NULL,

                         returnTrust=FALSE,
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
                         calcTables=TRUE, compress=FALSE,
                         covMethod=c("r", ""),
                         adjObf=TRUE, ci=0.95, sigdig=3, sigdigTable=NULL,
                         boundedTransform=TRUE, ...) {

  # trust_solve_c()'s fterm/mterm from sigdig, matching bobyqaControl()'s
  # rhoend and foceiControl()'s inner epsilon; a user-supplied value wins,
  # sigdig=NULL keeps the historic .Machine-derived fallback.
  if (is.null(fterm)) {
    fterm <- if (!is.null(sigdig)) .sigdigOptTol(sigdig) else 1.4901161193847656e-08
  }
  if (is.null(mterm)) mterm <- fterm
  checkmate::assertNumeric(fterm, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(mterm, len=1, any.missing=FALSE, lower=0)
  checkmate::assertIntegerish(iterlim, len=1, any.missing=FALSE, lower=1)

  # rinit/rmax cannot be numerically resolved here (they need the scaled
  # starting vector, only known once .trustFitModel() calls .nlmSetupEnv());
  # validated the same way foceiControl()'s trustRinit/trustRmax are.
  checkmate::assertNumeric(rinit, lower=0, finite=TRUE, null.ok=TRUE, len=1)
  if (!is.null(rinit) && rinit <= 0) {
    stop("'rinit' must be > 0", call.=FALSE)
  }
  checkmate::assertNumeric(rmax, lower=0, finite=TRUE, null.ok=TRUE, len=1)
  if (!is.null(rmax) && rmax <= 0) {
    stop("'rmax' must be > 0", call.=FALSE)
  }
  if (!is.null(rinit) && !is.null(rmax) && rinit > rmax) {
    stop("'rinit' cannot be larger than 'rmax'", call.=FALSE)
  }

  if (is.null(hessErr)) hessErr <- (.Machine$double.eps)^(1/3)
  checkmate::assertNumeric(hessErr, len=1, any.missing=FALSE, lower=0)
  checkmate::assertIntegerish(optimHessType, len=1, any.missing=FALSE, lower=1, upper=2)
  checkmate::assertIntegerish(shi21maxHess, len=1, any.missing=FALSE, lower=1)

  checkmate::assertLogical(returnTrust, len=1, any.missing=FALSE)
  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFixRes, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
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
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- .rxControlScaleSigdig(rxode2::rxControl(sigdig=sigdig), sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol=1e-4, rtol=1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- .rxControlScaleSigdig(do.call(rxode2::rxControl, rxControl), sigdig, skip = names(rxControl))
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
    rinit=rinit,
    rmax=rmax,
    iterlim=as.integer(iterlim),
    fterm=fterm,
    mterm=mterm,
    optimHessType=as.integer(optimHessType),
    shi21maxHess=as.integer(shi21maxHess),
    hessErr=hessErr,

    returnTrust=returnTrust,
    covMethod=match.arg(covMethod),
    optExpression=optExpression,
    literalFix=literalFix,
    literalFixRes=literalFixRes,
    sumProd=sumProd,
    rxControl=rxControl,
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
    calcTables=calcTables,
    compress=compress,
    ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
    genRxControl=.genRxControl,
    boundedTransform=boundedTransform)
  class(.ret) <- "trustControl"
  .ret
}

#' @export
rxUiDeparse.trustControl <- function(object, var) {
  .default <- trustControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' Get the trust family control
#'
#' @param env trust optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.trustFamilyControl <- function(env, ...) {
  .nlmFamilyControlGeneric(env, nlmixr2est::trustControl, "trustControl")
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.trustControl <- function(control, env) {
  assign("trustControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.trust <- function(x, ...) {
  .env <- x[[1]]
  if (exists("trustControl", .env)) {
    .control <- get("trustControl", .env)
    if (inherits(.control, "trustControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "trustControl")) return(.control)
  }
  stop("cannot find trust related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.trust <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- trustControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("trustControl", .ctl)
  if (!inherits(.ctl, "trustControl")) {
    .minfo("invalid control for `est=\"trust\"`, using default")
    .ctl <- trustControl()
  } else {
    .ctl <- do.call(trustControl, .ctl)
  }
  .ctl
}

.trustControlToFoceiControl <- function(env, assign=TRUE) {
  .trustControl <- env$trustControl
  .foceiControl <- foceiControl(rxControl=.trustControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.trustControl$sumProd,
                                optExpression=.trustControl$optExpression,
                                literalFix=.trustControl$literalFix,
                                literalFixRes=.trustControl$literalFixRes,
                                scaleTo=0,
                                calcTables=.trustControl$calcTables,
                                addProp=.trustControl$addProp,
                                interaction=0L,
                                compress=.trustControl$compress,
                                ci=.trustControl$ci,
                                sigdigTable=.trustControl$sigdigTable,
                                indTolRelax=.trustControl$indTolRelax,
                                eventSens=.trustControl$eventSens)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' Warn when a trust result's Newton decrement contradicts trust_solve_c()'s
#' own converged flag
#'
#' Split out of `.trustFitModel()` so the warning logic is testable against a
#' synthetic `tres` list without needing a full model fit -- the false-
#' positive-convergence case this guards against (RcppTrust reports
#' `converged=TRUE` off its own internal step-size tolerance while the true
#' Newton step at that point is still large, verified in `nlmTrustFit()`,
#' `src/nlm.cpp`) is a genuine numerical edge case, not reliably reproducible
#' on demand from a real fit.
#'
#' @param tres the list returned by `nlmTrustFit()`
#' @return `NULL`, called for the `warning()` side effect
#' @author Matthew L. Fidler
#' @noRd
.trustWarnUnderConverged <- function(tres) {
  if (isTRUE(tres$underConverged)) {
    warning("trust exited without a verified stationary point (the Newton step ",
            "at the reported optimum still exceeds tolerance); consider a ",
            "larger 'rmax'/'iterlim', or est=\"bobyqa\"/\"nlm\"", call.=FALSE)
  }
  invisible(NULL)
}

.trustFitModel <- function(ui, dataSav) {
  .ctl <- ui$control
  class(.ctl) <- NULL
  .p <- setNames(ui$nlmParIni, ui$nlmParName)
  .mi <- ui$nlmSensModel # trust always needs the gradient/Hessian model
  .env <- .nlmSetupEnv(.p, ui, dataSav, .mi, .ctl,
                       lower=ui$optimParLower, upper=ui$optimParUpper)
  on.exit({.nlmFreeEnv()})

  # rinit/rmax need the scaled starting vector -- not available inside
  # trustControl() itself. rinit mirrors minqa::bobyqa()'s own default-rhobeg
  # formula; rmax reuses the 8x growth-ceiling convention already established
  # for foceiControl(innerOpt="trust")'s per-subject eta trust region.
  .rinit <- .ctl$rinit
  if (is.null(.rinit)) .rinit <- min(0.95, 0.2*max(abs(.env$par.ini)))
  .rmax <- .ctl$rmax
  if (is.null(.rmax)) .rmax <- 8 * .rinit

  .tctl <- list(rinit=.rinit, rmax=.rmax, iterlim=.ctl$iterlim,
               fterm=.ctl$fterm, mterm=.ctl$mterm)
  .tres <- nlmTrustFit(.env$par.ini, .tctl)
  .trustWarnUnderConverged(.tres)
  if (isTRUE(.ctl$returnTrust)) return(.tres)

  .ret <- list(par=setNames(.tres$par, NULL),
              fval=.tres$value,
              hessian=matrix(.tres$hessian, nrow=length(.p), ncol=length(.p)),
              convergence=if (isTRUE(.tres$converged)) 0L else 1L,
              iterations=.tres$iterations,
              message=if (isTRUE(.tres$converged)) "converged" else "did not fully converge")
  .nlmFinalizeList(.env, .ret, par="par", printLine=TRUE,
                   hessianCov=TRUE)
}

#' Get the full theta for the trust method
#'
#' @param nlm enhanced trust return
#' @param ui ui object
#' @return named theta matrix
#' @author Matthew L. Fidler
#' @noRd
.trustGetTheta <- function(nlm, ui) {
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

.trustFamilyFit <- function(env, ...) {
  .nlmFamilyFitGeneric(
    env, "trust", .trustFitModel, .trustGetTheta,
    objective = function(.fit) 2 * as.numeric(.fit$fval),
    controlToFocei = .trustControlToFoceiControl,
    returnFlag = "returnTrust")
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.trust <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'trust', try 'focei'",
                                   .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'trust'",
                                   .var.name=.ui$modelName)
  rxode2::warnRxBounded(.ui, " which are ignored in 'trust'",
                        .var.name=.ui$modelName)
  .trustFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .trustFamilyFit(env, ...)
}
attr(nlmixr2Est.trust, "covPresent") <- TRUE
attr(nlmixr2Est.trust, "unbounded") <- TRUE
