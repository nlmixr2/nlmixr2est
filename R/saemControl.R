#' Control Options for SAEM
#'
#' @param seed Random Seed for SAEM step.  (Needs to be set for
#'     reproducibility.)  By default this is 99.
#'
#' @param nBurn Number of iterations in the first phase, ie the  MCMC/Stochastic Approximation
#'     steps. This is equivalent to Monolix's \code{K_0} or \code{K_b}.
#'
#' @param nEm Number of iterations in the Expectation-Maximization
#'     (EM) Step. This is equivalent to Monolix's \code{K_1}.
#'
#' @param nmc Number of Markov Chains. By default this is 3.  When
#'     you increase the number of chains the numerical integration by
#'     MC method will be more accurate at the cost of more
#'     computation.  In Monolix this is equivalent to \code{L}.
#'
#' @param nu This is a vector of 3 integers. They represent the
#'     numbers of transitions of the three different kernels used in
#'     the Hasting-Metropolis algorithm.  The default value is \code{c(2,2,2)},
#'     representing 40 for each transition initially (each value is
#'     multiplied by 20).
#'
#'     The first value represents the initial number of multi-variate
#'     Gibbs samples are taken from a normal distribution.
#'
#'     The second value represents the number of uni-variate, or multi-
#'     dimensional random walk Gibbs samples are taken.
#'
#'     The third value represents the number of bootstrap/reshuffling or
#'     uni-dimensional random samples are taken.
#'
#' @inheritParams iterPrintParams
#'
#' @param trace An integer indicating if you want to trace(1) the
#'     SAEM algorithm process.  Useful for debugging, but not for
#'     typical fitting.
#'
#' @param covMethod Method for calculating covariance.  In this
#'     discussion, R is the Hessian matrix of the objective
#'     function. The S matrix is the sum of each individual's
#'     gradient cross-product (evaluated at the individual empirical
#'     Bayes estimates).
#'
#'  "\code{linFim}" Use the Linearized Fisher Information Matrix to calculate the covariance.
#'
#'  "\code{fim}" Use the Fisher Information Matrix accumulated during SAEM
#'  estimation to calculate the covariance.  Like \code{sa} it inverts the observed
#'  information to a full theta + \code{Omega} diagonal + residual covariance, but
#'  uses the (noisier) estimation-phase matrix rather than a dedicated cov phase.
#'
#'  "\code{sa}" Use the stochastic-approximation Fisher Information Matrix.  After
#'  estimation, a dedicated covariance phase (\code{nSaCov} iterations) holds the
#'  parameters at the converged estimate and keeps resimulating the individual
#'  parameters, Monte-Carlo averaging the Louis observed-information integrand into a
#'  converged FIM decoupled from the cooling schedule (the approach used by Monolix;
#'  Kuhn & Lavielle 2005).  Always includes every estimated population parameter
#'  (theta, the \code{Omega} diagonal variances, and residual).
#'
#'  "\code{r,s}" Uses the sandwich matrix to calculate the covariance, that is: \eqn{R^-1 \times S \times R^-1}
#'
#'  "\code{r}" Uses the Hessian matrix to calculate the covariance as \eqn{2\times R^-1}
#'
#'  "\code{s}" Uses the crossproduct matrix to calculate the covariance as \eqn{4\times S^-1}
#'
#'  "" Does not calculate the covariance step.
#'
#' @param covFull Boolean (default \code{TRUE}) indicating the covariance
#'   should include every estimated population parameter -- the structural and
#'   residual thetas plus the \code{Omega} variance/covariance elements -- named
#'   \code{om.<eta>} / \code{cov.<eta>.<eta>}.  When \code{FALSE} the legacy
#'   structural-theta-only covariance is reported.  Ignored by
#'   \code{covMethod="sa"}, which is always full.
#'
#' @param nSaCov Number of iterations in the dedicated stochastic-approximation
#'   covariance phase used by \code{covMethod="sa"} (default \code{500}).  These
#'   iterations run at the converged estimate (parameters frozen) and only
#'   resimulate the individual parameters to build the observed Fisher
#'   information; a larger value gives a less noisy covariance.  Ignored by other
#'   covariance methods.
#'
#' @param logLik boolean indicating that log-likelihood should be
#'     calculate by Gaussian quadrature.
#'
#' @param nnodesGq number of nodes to use for the Gaussian
#'     quadrature when computing the likelihood with this method
#'     (defaults to 1, equivalent to the Laplacian likelihood)
#'
#' @param nsdGq span (in SD) over which to integrate when computing
#'     the likelihood by Gaussian quadrature. Defaults to 3 (eg 3
#'     times the SD)
#'
#' @param adjObf is a boolean to indicate if the objective function
#'     should be adjusted to be closer to NONMEM's default objective
#'     function.  By default this is \code{TRUE}
#'
#' @param tol This is the tolerance for the regression models used
#'   for complex residual errors (ie add+prop etc)
#'
#' @param itmax This is the maximum number of iterations for the
#'   regression models used for complex residual errors.  The number
#'   of iterations is itmax*number of parameters
#'
#' @param type indicates the type of optimization for the residuals; Can be one of c("nelder-mead", "newuoa")
#' @param powRange This indicates the range that powers can take for residual errors;  By default this is 10 indicating the range is c(-10, 10)
#' @param lambdaRange This indicates the range that Box-Cox and Yeo-Johnson parameters are constrained to be;  The default is 3 indicating the range c(-3,3)
#'
#' @param perSa This is the percent of the time the `nBurn`
#'   iterations in phase runs runs a simulated annealing.
#'
#' @param perNoCor This is the percentage of the MCMC phase of the SAEM
#'   algorithm where the variance/covariance matrix has no
#'   correlations. By default this is 0.75 or 75% of the
#'   Monte-carlo iteration.
#'
#'
#' @param perFixOmega This is the percentage of the `nBurn` phase
#'   where the omega values are unfixed to allow better exploration
#'   of the likelihood surface.  After this time, the omegas are
#'   fixed during optimization.
#'
#' @param perFixResid This is the percentage of the `nBurn` phase
#'   where the residual components are unfixed to allow better
#'   exploration of the likelihood surface.
#'
#' @param muRefCov This controls if mu-referenced covariates in `saem`
#'   are handled differently than non mu-referenced covariates.  When
#'   `TRUE`, mu-referenced covariates have special handling.  When
#'   `FALSE` mu-referenced covariates are treated the same as any
#'   other input parameter.
#'
#' @param muRefCovAlg This controls if algebraic expressions that can
#'   be mu-referenced are treated as mu-referenced covariates by:
#'
#'   1. Creating a internal data-variable `nlmixrMuDerCov#` for each
#'      algebraic mu-referenced expression
#'
#'   2. Change the algebraic expression to `nlmixrMuDerCov# * mu_cov_theta`
#'
#'   3. Use the internal mu-referenced covariate for saem
#'
#'   4. After optimization is completed, replace `model({})` with old
#'   `model({})` expression
#'
#'   5. Remove `nlmixrMuDerCov#` from nlmix2 output
#'
#' In general, these covariates should be more accurate since it
#' changes the system to a linear compartment model.  Therefore, by default this is `TRUE`.
#'
#' @param handleUninformativeEtas boolean that tells nlmixr2's saem to
#'   calculate uninformative etas and handle them specially (default
#'   is `TRUE`).
#'
#' @param mixProbMethod For mixture models (`mix()`, more than one
#'   component), stabilizes the mixing-probability estimate against
#'   collapsing onto a single component (the responsibility used to
#'   update it is itself weighted by the current mixing probability,
#'   which can create a runaway feedback loop). Two options:
#'
#'   * `"regularized"` (default): blend `mixProbPriorN` pseudo-subjects,
#'     distributed per the initial mixing probability, into the
#'     responsibility average each iteration (Dirichlet/MAP-EM-style).
#'     Prevents collapse even in difficult cases, at the cost of some
#'     bias toward the initial guess; may need larger `nBurn`/`nEm`.
#'
#'   * `"annealed"`: give the mixing-probability update its own decaying
#'     step-size schedule (`mixProbStepExp`) instead of the
#'     full-replacement step used during `nBurn`. Lower bias, but does
#'     not by itself fix a systematic (non-noise-driven) collapse.
#'
#' @param mixProbStepExp Only used when `mixProbMethod="annealed"`. Decay
#'   exponent for the mixing-probability step size
#'   (`1/iteration^mixProbStepExp`), applied from iteration 1. Default 1;
#'   smaller values decay more slowly.
#'
#' @param mixProbPriorN Only used when `mixProbMethod="regularized"`.
#'   Number of pseudo-subjects blended into the responsibility average
#'   each iteration. Larger values are more robust to collapse but bias
#'   the estimate more and need more `nBurn`/`nEm`. Default 20.
#'
#' @param mixSampleMethod For mixture models with per-component etas
#'   (split-ETA, e.g. `cl <- mix(tcl1 + eta.cl1, p1, tcl2 + eta.cl2)`),
#'   controls the MCMC/sufficient-statistic architecture for the
#'   individual random effects, independent of `mixProbMethod`. BSV
#'   (`$omega`) for split components is unreliable under `"parallel"`
#'   regardless of `mixProbMethod`.
#'
#'   * `"parallel"` (default): one full MCMC chain per component per
#'     subject per iteration, blended post hoc by responsibility. Mirrors
#'     NONMEM's `$MIX` and correctly estimates BSV shared across
#'     components, but cannot cleanly separate per-component BSV for
#'     split-ETA models (each "wrong-hypothesis" chain still explores its
#'     non-owned column(s) as unconstrained prior noise).
#'
#'   * `"msaem"` (experimental): the MSAEM algorithm (Lavielle & Mbogning
#'     2014), as used by Monolix. Simulates one random-effects trajectory
#'     per subject per iteration (label marginalized out via a
#'     closed-form responsibility) instead of parallel per-component
#'     chains, so no post-hoc blending is needed. Not compute-matched to
#'     `"parallel"` at equal `nmc` -- set `nmc` to roughly `nMix` times
#'     its default for a fair comparison. Uses a model-aware stratified
#'     initialization for split-ETA components that reliably achieves
#'     full theta/fixed-effect separation. Split-ETA BSV recovery is
#'     improved (two numerical bugs fixed: an `IGamma2_phi1` blowup that
#'     locked variance to exactly zero, and an inverted responsibility
#'     sign) but still not reliable -- it often settles at a safety-floor
#'     value rather than the true variance. Prefer `"parallel"` unless
#'     specifically evaluating this method.
#'
#' @param ... Other arguments to control SAEM.
#'
#' @inheritParams rxode2::rxSolve
#' @inheritParams foceiControl
#' @return List of options to be used in \code{\link{nlmixr2}} fit for
#'     SAEM.
#' @author Wenping Wang & Matthew L. Fidler
#' @references
#' Kuhn E, Lavielle M (2005). "Maximum likelihood estimation in nonlinear mixed
#' effects models." Computational Statistics & Data Analysis, 49(4), 1020-1038.
#' \doi{10.1016/j.csda.2004.07.002}
#'
#' Jiang L, Roy A, Balasubramanian K, Davis D, Drusvyatskiy D, Na S (2025).
#' "Online Covariance Estimation in Nonsmooth Stochastic Approximation."
#' arXiv:2502.05305. \doi{10.48550/arXiv.2502.05305}
#' @family Estimation control
#' @export
saemControl <- function(seed = 99,
                        nBurn = 200,
                        nEm = 300,
                        nmc = 3,
                        nu = c(2, 2, 2),
                        print = 1L,
                        trace = 0, # nolint
                        covMethod = c("linFim", "fim", "sa", "r,s", "r", "s", ""),
                        covFull = TRUE,
                        nSaCov = 500L,
                        calcTables = TRUE,
                        logLik = FALSE,
                        nnodesGq = 3,
                        nsdGq = 1.6,
                        optExpression = TRUE,
                        literalFix=FALSE,
                        adjObf = TRUE,
                        sumProd = FALSE,
                        addProp = c("combined2", "combined1"),
                        tol = 1e-6,
                        itmax = 30,
                        type = c("nelder-mead", "newuoa"),
                        powRange = 10,
                        lambdaRange = 3,
                        odeRecalcFactor=10^(0.5),
                        maxOdeRecalc=5L,
                        indTolRelax=TRUE,
                        perSa=0.75,
                        perNoCor=0.75,
                        perFixOmega=0.1,
                        perFixResid=0.1,
                        compress=TRUE,
                        rxControl=NULL,
                        sigdig=NULL,
                        sigdigTable=NULL,
                        ci=0.95,
                        muRefCov=TRUE,
                        muRefCovAlg=TRUE,
                        handleUninformativeEtas=TRUE,
                        iovXform = c("sd", "var", "logsd", "logvar"),
                        boundedTransform = TRUE,
                        eventSens = c("jump", "fd"),
                        mixProbMethod = c("regularized", "annealed"),
                        mixProbStepExp = 1,
                        mixProbPriorN = 20,
                        mixSampleMethod = c("parallel", "msaem"),
                        ...) {
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl", "mcmc",
                             "DEBUG", "iterPrintControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  iovXform <- match.arg(iovXform)
  checkmate::assertIntegerish(seed, any.missing=FALSE, min.len=1)
  if (!is.null(.xtra$mcmc)) {
    #mcmc = list(niter = c(nBurn, nEm), nmc = nmc, nu = nu),
    checkmate::assertIntegerish(.xtra$mcmc$niter, len=2, lower=0, any.missing=FALSE, .var.name="mcmc$niter")
    nBurn <- .xtra$mcmc$niter[1]
    nEm   <- .xtra$mcmc$niter[2]
    checkmate::assertIntegerish(.xtra$mcmc$nmc, len=1, lower=1, any.missing=FALSE, .var.name="mcmc$nmc")
    nmc <- .xtra$mcmc$nmc
    checkmate::assertIntegerish(.xtra$mcmc$nu, len=3, lower=1, any.missing=FALSE, .var.name="mcmc$nu")
  }
  checkmate::assertIntegerish(nBurn, any.missing=FALSE, len=1, lower=0)
  checkmate::assertIntegerish(nEm, any.missing=FALSE, len=1, lower=0)
  checkmate::assertIntegerish(nmc, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(nu, any.missing=FALSE, len=3, lower=1)
  # `print` can be either a scalar print-frequency or a pre-built
  # iterPrintControl object; .absorbIterPrintControl validates either form
  # and returns the canonical iterPrintControl list.  list(...)$iterPrintControl
  # catches the round-trip case where the previous saemControl()'s return
  # value is passed back through do.call(saemControl, .ctl).
  .iterPrintControl <- .absorbIterPrintControl(print = print,
                                               iterPrintControl = .xtra$iterPrintControl)
  if (!is.null(.xtra$DEBUG)) {
    trace <- .xtra$DEBUG # nolint
  }
  checkmate::assertIntegerish(trace, any.missing=FALSE, lower=0, upper=1, len=1) # nolint
  checkmate::assertLogical(calcTables, any.missing=FALSE, len=1)
  checkmate::assertLogical(logLik, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(nnodesGq, any.missing=FALSE, lower=1, len=1)
  checkmate::assertNumeric(nsdGq, any.missing=FALSE, lower=1, len=1, finite=TRUE)
  checkmate::assertLogical(optExpression, any.missing=FALSE, len=1)
  checkmate::assertLogical(literalFix, any.missing=FALSE, len=1)
  checkmate::assertLogical(adjObf, any.missing=FALSE, len=1)
  checkmate::assertLogical(sumProd, any.missing=FALSE, len=1)
  checkmate::assertNumeric(tol, any.missing=FALSE, len=1, finite=TRUE)
  checkmate::assertIntegerish(itmax, any.missing=FALSE, len=1, lower=1)
  checkmate::assertNumeric(powRange, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(lambdaRange, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(odeRecalcFactor, any.missing=FALSE, lower=0, len=1, finite=TRUE)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, lower=0, len=1)
  checkmate::assertLogical(indTolRelax, any.missing=FALSE, len=1)
  checkmate::assertNumeric(perSa, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertNumeric(perNoCor, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertNumeric(perFixOmega, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertNumeric(perFixResid, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertLogical(muRefCov, any.missing=FALSE, len=1)
  checkmate::assertLogical(muRefCovAlg, any.missing=FALSE, len=1)
  checkmate::assertLogical(handleUninformativeEtas, any.missing=FALSE, len=1)
  checkmate::assertLogical(boundedTransform, any.missing=FALSE, len=1)
  eventSens <- match.arg(eventSens)
  mixProbMethod <- match.arg(mixProbMethod)
  checkmate::assertNumeric(mixProbStepExp, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(mixProbPriorN, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  mixSampleMethod <- match.arg(mixSampleMethod)

  type <- match.arg(type)
  if (inherits(addProp, "numeric")) {
    if (addProp == 1) {
      addProp <- "combined1"
    } else if (addProp == 2) {
      addProp <- "combined2"
    } else {
      stop("addProp must be 1, 2, \"combined1\" or \"combined2\"", call.=FALSE)
    }
  } else {
    addProp <- match.arg(addProp)
  }
  checkmate::assertLogical(compress, any.missing=FALSE, len=1)
  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
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
  .env <- nlmixr2global$nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  if (is.null(rxControl)) {
    rxControl <- rxode2::rxControl(sigdig=sigdig, envir=.env)
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
    rxControl$envir <- .env
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'", call=FALSE)
  }

  if (checkmate::testIntegerish(covMethod, lower=0, len=1, any.missing=FALSE)) {
    .covMethod <- covMethod
  } else {
    .covMethod <- match.arg(covMethod)
  }

  checkmate::assertLogical(covFull, len=1, any.missing=FALSE)

  .ret <- list(
    mcmc = list(niter = c(nBurn, nEm), nmc = nmc, nu = nu),
    rxControl = rxControl,
    seed = seed,
    iterPrintControl = .iterPrintControl,
    DEBUG = trace, # nolint
    optExpression = optExpression,
    literalFix=literalFix,
    sumProd = sumProd,
    nnodesGq = nnodesGq,
    nsdGq = nsdGq,
    adjObf = adjObf,
    addProp = addProp,
    itmax = itmax,
    tol = tol,
    type = type,
    powRange = powRange,
    lambdaRange = lambdaRange,
    odeRecalcFactor=odeRecalcFactor,
    maxOdeRecalc=maxOdeRecalc,
    indTolRelax=indTolRelax,
    perSa=perSa,
    perNoCor=perNoCor,
    perFixOmega=perFixOmega,
    perFixResid=perFixResid,
    compress=compress,
    genRxControl=.genRxControl,
    sigdigTable=sigdigTable,
    ci=ci,
    covMethod=.covMethod,
    covFull=covFull,
    nSaCov=as.integer(nSaCov),
    logLik=logLik,
    calcTables=calcTables,
    muRefCov=muRefCov,
    muRefCovAlg=muRefCovAlg,
    handleUninformativeEtas=handleUninformativeEtas,
    iovXform=iovXform,
    boundedTransform=boundedTransform,
    eventSens=eventSens,
    mixProbMethod=mixProbMethod,
    mixProbStepExp=mixProbStepExp,
    mixProbPriorN=mixProbPriorN,
    mixSampleMethod=mixSampleMethod
  )
  class(.ret) <- "saemControl"
  .ret
}

.saemDeparseExtra <- function(default, name, value) {
  if (name == "mcmc") {
    .ret <- character(0)
    if (!identical(default$mcmc$niter, value$niter)) {
      if (default$mcmc$niter[1] != value$niter[1]) {
        .ret <- c(.ret, paste0("nBurn=", value$niter[1]))
      }
      if (default$mcmc$niter[2] != value$niter[2]) {
        .ret <- c(.ret, paste0("nEm=", value$niter[2]))
      }
    }
    if (default$mcmc$nmc != value$nmc) {
      .ret <- c(.ret, paste0("nmc=", value$nmc))
    }
    if (!identical(default$mcmc$nu, value$nu)) {
      .ret <- c(.ret, paste0("nu=", deparse1(value$nu)))
    }
    return(paste0(.ret, collapse=","))
  }
  NA_character_
}

#' @export
rxUiDeparse.saemControl <- function(object, var) {
  .default <- saemControl()
  .w <- .deparseDifferent(.default, object, c("genRxControl", "DEBUG"))
  .deparseFinal(.default, object, .w, var, fun=.saemDeparseExtra)
}
