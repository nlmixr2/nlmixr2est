#' Control for advi (automatic differentiation variational inference) in nlmixr2
#'
#' Variational-inference NLME estimation following Kucukelbir et al. (2017): the
#' latent variables are transformed to an unconstrained real coordinate space, a
#' Gaussian variational family is posited there, and the ELBO is maximized by
#' stochastic gradient ascent using the reparameterization trick.  The gradient
#' of the log-joint is obtained from the FOCEi forward sensitivities (inner
#' per-subject eta gradient and the outer population sensitivity contraction),
#' not from automatic differentiation.  The whole optimization loop runs in C++.
#'
#' @inheritParams saemControl
#' @inheritParams foceiControl
#'
#' @param seed Random seed for the ADVI optimization (reparameterization
#'   sampling); default 42.  The Monte-Carlo gradient is stochastic, so a fixed
#'   seed makes every fit reproducible.  Reparameterization noise is drawn from a
#'   counter-based stream keyed by the global iteration index, so a shorter run
#'   is a bit-for-bit prefix of a longer one and results are independent of the
#'   number of cores.
#' @param iters Total number of ADVI (stochastic gradient ascent) iterations.
#' @param nMc Number of Monte-Carlo samples used to approximate the ELBO
#'   gradient at each iteration (the paper's `M`; typically 1-10).
#' @param adviFamily Variational family in the unconstrained space.
#'   `"fullRank"` (default) uses a block full-rank Gaussian: a dense
#'   `neta x neta` Cholesky factor per subject plus a dense block over the
#'   population vector (mean-field across blocks).  `"meanField"` uses a fully
#'   factorized (diagonal) Gaussian.  Mean-field is faster but is known to
#'   underestimate marginal variances.
#' @param pointEstimate When `TRUE` (default) run a variational-EM hybrid: the
#'   variational posterior covers the per-subject etas only, and the population
#'   parameters (thetas / omega / residual error) are point estimates maximized
#'   by the ADVI gradient (stochastic maximum likelihood); output semantics match
#'   FOCEi/SAEM.  When `FALSE` run full Bayes: the variational posterior also
#'   covers the unconstrained population vector, with flat priors.
#' @param optim Stochastic optimizer.  `"advi"` (default) uses the paper's
#'   adaptive step-size sequence (Eqs 10-11); `"adam"` uses Adam.
#' @param adaptEta When `TRUE` (default) adaptively choose the step-size scale
#'   `eta` by a short search over `etaCandidates` before the main loop; when
#'   `FALSE` use a fixed `eta` (the first `etaCandidates` entry).
#' @param etaCandidates Candidate step-size scales searched when `adaptEta` is
#'   `TRUE` (the paper searches `c(0.01, 0.1, 1, 10, 100)`).
#' @param tau Stabilizing constant `tau > 0` in the step-size denominator
#'   (paper Eq 10); the step-size is insensitive to it.
#' @param alpha Weighting `alpha` in (0, 1) of new vs old gradient information in
#'   the step-size memory recursion (paper Eq 11).
#' @param tol Convergence tolerance on the relative change in the ELBO; the loop
#'   may stop early once the change stays below this.  `0` disables early
#'   stopping (run all `iters`).
#' @param likelihood Inner likelihood used for the per-subject objective and
#'   gradient, run through the FOCEi inner interface: `"focei"` (default),
#'   `"foce"`, `"focep"`, or `"laplace"`.
#' @param returnAdvi When `TRUE` return the raw ADVI optimization object instead
#'   of the nlmixr2 fit.
#'
#' @return advi control structure (class `adviControl`)
#' @export
#' @author Matthew L. Fidler
adviControl <- function(seed = 42L,
                        iters = 300L,
                        nMc = 1L,
                        adviFamily = c("fullRank", "meanField"),
                        pointEstimate = TRUE,
                        optim = c("advi", "adam"),
                        adaptEta = TRUE,
                        etaCandidates = c(0.01, 0.1, 1, 10, 100),
                        tau = 1.0,
                        alpha = 0.1,
                        tol = 1e-4,
                        likelihood = c("focei", "foce", "focep", "laplace"),
                        returnAdvi = FALSE,

                        print = 1L,
                        useColor = NULL,
                        printNcol = NULL,

                        covMethod = c("advi", "analytic", "r,s", "r", "s", ""),
                        optExpression = TRUE,
                        sumProd = FALSE,
                        literalFix = TRUE,
                        literalFixRes = TRUE,
                        addProp = c("combined2", "combined1"),
                        calcTables = TRUE,
                        compress = FALSE,
                        adjObf = TRUE,
                        ci = 0.95,
                        sigdig = NULL,
                        sigdigTable = NULL,

                        stickyRecalcN = 4,
                        maxOdeRecalc = 5,
                        odeRecalcFactor = 10^(0.5),
                        indTolRelax = TRUE,
                        eventSens = c("jump", "fd"),
                        rxControl = NULL,
                        ...) {

  checkmate::assertIntegerish(seed, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(iters, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(nMc, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertLogical(pointEstimate, len = 1, any.missing = FALSE)
  checkmate::assertLogical(adaptEta, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(etaCandidates, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 1)
  checkmate::assertNumeric(tau, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(alpha, lower = 0, upper = 1, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(tol, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assertLogical(returnAdvi, len = 1, any.missing = FALSE)
  checkmate::assertLogical(optExpression, len = 1, any.missing = FALSE)
  checkmate::assertLogical(sumProd, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFix, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFixRes, len = 1, any.missing = FALSE)
  checkmate::assertLogical(calcTables, len = 1, any.missing = FALSE)
  checkmate::assertLogical(compress, len = 1, any.missing = TRUE)
  checkmate::assertLogical(adjObf, len = 1, any.missing = TRUE)
  checkmate::assertIntegerish(stickyRecalcN, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(odeRecalcFactor, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assertLogical(indTolRelax, len = 1, any.missing = FALSE)
  adviFamily <- match.arg(adviFamily)
  optim <- match.arg(optim)
  likelihood <- match.arg(likelihood)
  # match.arg cannot match ""; treat it (skip covariance) like foceiControl does
  if (length(covMethod) == 1L && covMethod == "") {
    covMethod <- ""
  } else {
    covMethod <- match.arg(covMethod)
  }
  addProp <- match.arg(addProp)
  eventSens <- match.arg(eventSens)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl", "iterPrintControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste(paste0("'", .bad, "'"), collapse = ", "),
         call. = FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- rxode2::rxControl(sigdig = sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol = 1e-4, rtol = 1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'",
         call. = FALSE)
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

  .iterPrintControl <- .absorbIterPrintControl(print = print,
                                               printNcol = printNcol,
                                               useColor = useColor,
                                               iterPrintControl = .xtra$iterPrintControl)

  .ret <- list(seed = as.integer(seed),
               iters = as.integer(iters),
               nMc = as.integer(nMc),
               adviFamily = adviFamily,
               pointEstimate = pointEstimate,
               optim = optim,
               adaptEta = adaptEta,
               etaCandidates = as.numeric(etaCandidates),
               tau = tau,
               alpha = alpha,
               tol = tol,
               likelihood = likelihood,
               returnAdvi = returnAdvi,
               covMethod = covMethod,
               optExpression = optExpression,
               sumProd = sumProd,
               literalFix = literalFix,
               literalFixRes = literalFixRes,
               addProp = addProp,
               calcTables = calcTables,
               compress = compress,
               adjObf = adjObf,
               ci = ci,
               sigdig = sigdig,
               sigdigTable = sigdigTable,
               stickyRecalcN = as.integer(stickyRecalcN),
               maxOdeRecalc = as.integer(maxOdeRecalc),
               odeRecalcFactor = odeRecalcFactor,
               indTolRelax = indTolRelax,
               eventSens = eventSens,
               iterPrintControl = .iterPrintControl,
               rxControl = rxControl,
               genRxControl = .genRxControl)
  class(.ret) <- "adviControl"
  .ret
}

#' @export
rxUiDeparse.adviControl <- function(object, var) {
  .default <- adviControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.adviControl <- function(control, env) {
  assign("adviControl", control, envir = env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.advi <- function(x, ...) {
  .env <- x[[1]]
  if (exists("adviControl", .env)) {
    .control <- get("adviControl", .env)
    if (inherits(.control, "adviControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "adviControl")) return(.control)
  }
  stop("cannot find advi related control object", call. = FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.advi <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- adviControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("adviControl", .ctl)
  if (!inherits(.ctl, "adviControl")) {
    .minfo("invalid control for `est=\"advi\"`, using default")
    .ctl <- adviControl()
  } else {
    .ctl <- do.call(adviControl, .ctl)
  }
  .ctl
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.advi <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'advi'", .var.name = .ui$modelName)
  ## absorb the validated control (set by getValidNlmixrControl before dispatch)
  if (exists("control", envir = env) && inherits(env$control, "adviControl")) {
    assign("adviControl", env$control, envir = env)
  } else {
    assign("adviControl", adviControl(), envir = env)
  }
  ## Seed the ENTIRE estimation ONCE here and restore the caller's global RNG
  ## state afterward; the counter-based reparameterization stream inside the loop
  ## makes a shorter run a bit-for-bit prefix of a longer one.
  rxode2::rxWithSeed(env$adviControl$seed, {
    .adviFitModel(env)
  })
}
attr(nlmixr2Est.advi, "covPresent") <- TRUE
## ADVI optimizes in the unconstrained real coordinate space
attr(nlmixr2Est.advi, "unbounded") <- TRUE
attr(nlmixr2Est.advi, "iov") <- TRUE
