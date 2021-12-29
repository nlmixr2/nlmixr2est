##' Control Options for SAEM
##'
##' @param seed Random Seed for SAEM step.  (Needs to be set for
##'     reproducibility.)  By default this is 99.
##'
##' @param nBurn Number of iterations in the Stochastic Approximation
##'     (SA), or burn-in step. This is equivalent to Monolix's \code{K_0} or \code{K_b}.
##'
##' @param nEm Number of iterations in the Expectation-Maximization
##'     (EM) Step. This is equivalent to Monolix's \code{K_1}.
##'
##' @param nmc Number of Markov Chains. By default this is 3.  When
##'     you increase the number of chains the numerical integration by
##'     MC method will be more accurate at the cost of more
##'     computation.  In Monolix this is equivalent to \code{L}
##'
##' @param nu This is a vector of 3 integers. They represent the
##'     numbers of transitions of the three different kernels used in
##'     the Hasting-Metropolis algorithm.  The default value is \code{c(2,2,2)},
##'     representing 40 for each transition initially (each value is
##'     multiplied by 20).
##'
##'     The first value represents the initial number of multi-variate
##'     Gibbs samples are taken from a normal distribution.
##'
##'     The second value represents the number of uni-variate, or multi-
##'     dimensional random walk Gibbs samples are taken.
##'
##'     The third value represents the number of bootstrap/reshuffling or
##'     uni-dimensional random samples are taken.
##'
##' @param print The number it iterations that are completed before
##'     anything is printed to the console.  By default, this is 1.
##'
##' @param covMethod  Method for calculating covariance.  In this
##'     discussion, R is the Hessian matrix of the objective
##'     function. The S matrix is the sum of each individual's
##'     gradient cross-product (evaluated at the individual empirical
##'     Bayes estimates).
##'
##'  "\code{linFim}" Use the Linearized Fisher Information Matrix to calculate the covariance.
##'
##'  "\code{fim}" Use the SAEM-calculated Fisher Information Matrix to calculate the covariance.
##'
##'  "\code{r,s}" Uses the sandwich matrix to calculate the covariance, that is: \eqn{R^-1 \times S \times R^-1}
##'
##'  "\code{r}" Uses the Hessian matrix to calculate the covariance as \eqn{2\times R^-1}
##'
##'  "\code{s}" Uses the crossproduct matrix to calculate the covariance as \eqn{4\times S^-1}
##'
##'  "" Does not calculate the covariance step.
##'
##' @param logLik boolean indicating that log-likelihood should be
##'     calculate by Gaussian quadrature.
##'
##' @param trace An integer indicating if you want to trace(1) the
##'     SAEM algorithm process.  Useful for debugging, but not for
##'     typical fitting.
##'
##' @param nnodes.gq number of nodes to use for the Gaussian
##'     quadrature when computing the likelihood with this method
##'     (defaults to 1, equivalent to the Laplaclian likelihood)
##'
##' @param nsd.gq span (in SD) over which to integrate when computing
##'     the likelihood by Gaussian quadrature. Defaults to 3 (eg 3
##'     times the SD)
##'
##' @param adjObf is a boolean to indicate if the objective function
##'     should be adjusted to be closer to NONMEM's default objective
##'     function.  By default this is \code{TRUE}
##'
##' @param tol This is the tolerance for the regression models used
##'   for complex residual errors (ie add+prop etc)
##'
##' @param itmax This is the maximum number of iterations for the
##'   regression models used for complex residual errors.  The number
##'   of iterations is itmax*number of parameters
##'
##' @param ... Other arguments to control SAEM.
##'
##' @inheritParams rxode2::rxSolve
##' @inheritParams foceiControl
##' @inheritParams configsaem
##' @inheritParams rxode2::rxSEinner
##' @inheritParams rxode2::rxGenSaem
##' @return List of options to be used in \code{\link{nlmixr2}} fit for
##'     SAEM.
##' @author Wenping Wang & Matthew L. Fidler
##' @export
saemControl <- function(seed = 99,
                        nBurn = 200, nEm = 300,
                        nmc = 3,
                        nu = c(2, 2, 2),
                        atol = 1e-06,
                        rtol = 1e-04,
                        method = "liblsoda",
                        transitAbs = FALSE,
                        print = 1,
                        trace = 0,
                        covMethod = c("linFim", "fim", "r,s", "r", "s", ""),
                        calcTables = TRUE,
                        logLik = FALSE,
                        nnodes.gq = 3,
                        nsd.gq = 1.6,
                        optExpression = FALSE,
                        maxsteps = 100000L,
                        adjObf = TRUE,
                        sum.prod = FALSE,
                        addProp = c("combined2", "combined1"),
                        singleOde = TRUE,
                        tol = 1e-6,
                        itmax = 30,
                        type = c("nelder-mead", "newuoa"),
                        powRange = 10,
                        lambdaRange = 3,
                        loadSymengine=FALSE,
                        odeRecalcFactor=10^(0.5),
                        maxOdeRecalc=5L,
                        ...) {
  type <- match.arg(type)
  .xtra <- list(...)
  .rm <- c()
  if (missing(transitAbs) && !is.null(.xtra$transit_abs)) {
    transitAbs <- .xtra$transit_abs
    .rm <- c(.rm, "transit_abs")
  }
  if (missing(nBurn) && !is.null(.xtra$n.burn)) {
    nBurn <- .xtra$n.burn
    .rm <- c(.rm, "n.burn")
  }
  if (missing(nEm) && !is.null(.xtra$n.em)) {
    nEm <- .xtra$n.em
    .rm <- c(.rm, "n.em")
  }
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
  .ret <- list(
    mcmc = list(niter = c(nBurn, nEm), nmc = nmc, nu = nu),
    ODEopt = rxode2::rxControl(
      atol = atol, rtol = rtol, method = method,
      transitAbs = transitAbs, maxsteps = maxsteps, ...
    ),
    seed = seed,
    print = print,
    DEBUG = trace,
    optExpression = optExpression,
    sum.prod = sum.prod,
    nnodes.gq = nnodes.gq,
    nsd.gq = nsd.gq,
    adjObf = adjObf,
    addProp = addProp,
    singleOde = singleOde,
    itmax = itmax,
    tol = tol,
    type = type,
    powRange = powRange,
    lambdaRange = lambdaRange,
    loadSymengine=loadSymengine,
    odeRecalcFactor=odeRecalcFactor,
    maxOdeRecalc=maxOdeRecalc,
    ...
  )
  if (length(.rm) > 0) {
    .ret <- .ret[!(names(.ret) %in% .rm)]
  }
  .ret[["covMethod"]] <- match.arg(covMethod)
  .ret[["logLik"]] <- logLik
  .ret[["calcTables"]] <- calcTables
  class(.ret) <- "saemControl"
  .ret
}
#'  Determine if the parameter is a mu-referenced covariate
#'
#' @param expr Expression to check
#' @param muRefCovariateDataFrame Mu Ref data frame
#' @return A boolaen that tells if the expression is a mu-ref covariate
#' @author Matthew L. Fidler
#' @noRd
.saemDropParametersIsMuRefCovariate <- function(expr, muRefCovariateDataFrame) {
  if (length(expr) == 3) {
    if (identical(expr[[1]], quote(`*`))) {
      if (length(expr[[2]]) == 1 &&
            length(expr[[3]]) == 1) {
        .cov1 <- as.character(expr[[2]])
        .cov2 <- as.character(expr[[3]])
        .w <- which(muRefCovariateDataFrame$covariate == .cov1 &
                      muRefCovariateDataFrame$covariateParameter == .cov2)
        if (length(.w) == 1) return(TRUE)
        .w <- which(muRefCovariateDataFrame$covariate == .cov2 &
                      muRefCovariateDataFrame$covariateParameter == .cov1)
        if (length(.w) == 1) return(TRUE)
      }
    }
  }
  FALSE
}

.saemDropParameters <- function(line, muRefDataFrame, muRefCovariateDataFrame) {
  f <- function(x) {
    if (is.name(x) || is.atomic(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`+`))) {
        if (.saemDropParametersIsMuRefCovariate(x[[2]], muRefCovariateDataFrame)) {
          return(f(x[[3]]))
        }
        if (.saemDropParametersIsMuRefCovariate(x[[3]], muRefCovariateDataFrame)) {
          return(f(x[[2]]))
        }
        if (length(x[[2]]) == 1) {
          .char <- as.character(x[[2]])
          if (.char %in% muRefDataFrame$eta) {
            return(f(x[[3]]))
          }
        }
        if (length(x[[3]]) == 1) {
          .char <- as.character(x[[3]])
          if (.char %in% muRefDataFrame$eta) {
            return(f(x[[2]]))
          }
        }
      }
      as.call(lapply(x, f))
    } else {
      return(x)
    }
  }
  f(line)
}

.saemGenModel <- function(ui, timeVaryingCovariates=character(0)) {
  .muRefCovariateDataFrame <- ui$muRefCovariateDataFrame
  if (length(timeVaryingCovariates) > 0) {
    # Drop time-varying covariates
    .muRefCovariateDataFrame <- .muRefCovariateDataFrame[!(.muRefCovariateDataFrame$covariate %in% timeVaryingCovariates), ]
  }
}


#' @rdname nlmixr2Est
#' @export
nlmixr2Est.saem <- function(env, ...) {

}
