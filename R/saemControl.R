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
#' @param print The number it iterations that are completed before
#'     anything is printed to the console.  By default, this is 1.
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
#'  "\code{fim}" Use the SAEM-calculated Fisher Information Matrix to calculate the covariance.
#'
#'  "\code{r,s}" Uses the sandwich matrix to calculate the covariance, that is: \eqn{R^-1 \times S \times R^-1}
#'
#'  "\code{r}" Uses the Hessian matrix to calculate the covariance as \eqn{2\times R^-1}
#'
#'  "\code{s}" Uses the crossproduct matrix to calculate the covariance as \eqn{4\times S^-1}
#'
#'  "" Does not calculate the covariance step.
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
#' @param ... Other arguments to control SAEM.
#'
#' @inheritParams rxode2::rxSolve
#' @inheritParams foceiControl
#' @return List of options to be used in \code{\link{nlmixr2}} fit for
#'     SAEM.
#' @author Wenping Wang & Matthew L. Fidler
#' @family Estimation control
#' @export
saemControl <- function(seed = 99,
                        nBurn = 200,
                        nEm = 300,
                        nmc = 3,
                        nu = c(2, 2, 2),
                        print = 1,
                        trace = 0, # nolint
                        covMethod = c("linFim", "fim", "r,s", "r", "s", ""),
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
                        ...) {
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %fin% c("genRxControl", "mcmc",
                             "DEBUG"))]
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
  checkmate::assertIntegerish(print, any.missing=FALSE, lower=0, len=1)
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
  checkmate::assertNumeric(perSa, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertNumeric(perNoCor, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertNumeric(perFixOmega, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertNumeric(perFixResid, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertLogical(muRefCov, any.missing=FALSE, len=1)
  checkmate::assertLogical(muRefCovAlg, any.missing=FALSE, len=1)
  checkmate::assertLogical(handleUninformativeEtas, any.missing=FALSE, len=1)

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

  .ret <- list(
    mcmc = list(niter = c(nBurn, nEm), nmc = nmc, nu = nu),
    rxControl = rxControl,
    seed = seed,
    print = print,
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
    perSa=perSa,
    perNoCor=perNoCor,
    perFixOmega=perFixOmega,
    perFixResid=perFixResid,
    compress=compress,
    genRxControl=.genRxControl,
    sigdigTable=sigdigTable,
    ci=ci,
    covMethod=.covMethod,
    logLik=logLik,
    calcTables=calcTables,
    muRefCov=muRefCov,
    muRefCovAlg=muRefCovAlg,
    handleUninformativeEtas=handleUninformativeEtas,
    iovXform=iovXform
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
