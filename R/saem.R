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
#' Drop mu-referenced parameters
#'
#' @param line Line to change
#' @param muRefDataFrame Mu-referenced data frame
#' @param muRefCovariateDataFrame Mu Referenced Covariates
#' @return Remove mu-referenced etas and covariates
#' @author Matthew L. Fidler
#' @noRd
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
#' Drop mu referenced etas and covariates
#'
#' @param ui rxode2 ui
#' @return model line expression with mu referenced information dropped.
#' @author Matthew L. Fidler
#' @noRd
.saemDropMuRefFromModel <- function(ui) {
  if (exists("muRefFinal", ui)) {
    .muRefFinal <- ui$muRefFinal
  } else {
    .muRefFinal <- ui$muRefCovariateDataFrame
  }
  .muRefDataFrame <- ui$muRefDataFrame
  lapply(ui$lstExpr, function(line){
    .saemDropParameters(line, .muRefDataFrame, .muRefFinal)
  })
}
#' Change a character expression into a quoted name
#'
#' @param chr Character expression
#' @return Quote name
#' @author Matthew L. Fidler
#' @noRd
.enQuote <- function(chr) {
  eval(parse(text = paste0("quote(", chr, ")")))
}

#' This is a S3 method for getting the distribution lines for a base rxode2 saem problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the estimation of saem
#' @author Matthew Fidler
#' @keywords internal
#' @export
nmGetDistributionSaemLines <- function(line) {
  UseMethod("nmGetDistributionSaemLines")
}
#' Creates a saem line object from a predDf line
#'
#' @param x rxode2 ui object
#' @param line Line number for saem error line object
#' @param len Number of prediction statements
#' @return nmGetDistributionSaemLines object
#' @author Matthew L. Fidler
#' @noRd
.createSaemLineObject <- function(x, line) {
  .predDf <- get("predDf", x)
  if (line > nrow(.predDf)) {
    return(NULL)
  }
  .predLine <- .predDf[line, ]
  .ret <- list(x, .predLine)
  class(.ret) <- c(paste(.predLine$distribution), "nmGetDistributionSaemLines")
  .ret
}

#' @rdname nmGetDistributionSaemLines
#' @export
nmGetDistributionSaemLines.rxUi <- function(line) {
  .predDf <- get("predDf", line)
  lapply(seq_along(.predDf$cond), function(c){
    .mod <- .createSaemLineObject(line, c)
    nmGetDistributionSaemLines(.mod)
  })
}

#' @rdname nmGetDistributionSaemLines
#' @export
nmGetDistributionSaemLines.norm <- function(line) {
  .rx <- line[[1]]
  .pred1 <- line[[2]]
  if (.pred1[["linCmt"]]) {
    .var <- quote(linCmt())
  } else {
    .var <- .enQuote(.pred1[["var"]])
  }
  return(list(bquote(rx_pred_ <- .(.var))))
}

#' @export
nmGetDistributionSaemLines.t <- function(line) {
  stop("t isn't supported yet")
}

#' @export
nmGetDistributionSaemLines.default  <- function(line) {
  stop("Distribution not supported")
}

#' @export
rxUiGet.saemModel0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=nmGetDistributionSaemLines(.f),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE,
                              lstExpr=.saemDropMuRefFromModel(.f))
}
#attr(rxUiGet.saemModel0, "desc") <- "saem initial model"

#' Load the saem model into symengine
#'
#' @param x rxode2 UI object
#' @return String for loading into symengine
#' @author Matthew L. Fidler
#' @noRd
.saemPrune <- function(x) {
  .x <- x[[1]]
  .x <- .x$saemModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
    .malert("pruning branches ({.code if}/{.code else}) of saem model...")
  .ret <- rxode2::.rxPrune(.x, envir = .env)
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
rxUiGet.loadPruneSaem <- function(x, ...) {
  .loadSymengine(.saemPrune(x), promoteLinSens = FALSE)
}
#attr(rxUiGet.loadPruneSaem, "desc") <- "load the saem model into symengine"

#' @export
rxUiGet.saemParamsToEstimate <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  c(.iniDf$name[!is.na(.iniDf$ntheta) & is.na(.iniDf$err)], .ui$nonMuEtas)
}
attr(rxUiGet.saemParamsToEstimate, "desc") <- "Get the parameters to estimate"

#' @export
rxUiGet.saemParams <- function(x, ...) {
  .ui <- x[[1]]
  .par <- c(rxUiGet.saemParamsToEstimate(x, ...), .ui$covariates)
  paste0("params(", paste(.par, collapse=","), ")")
}
attr(rxUiGet.saemParams, "desc") <- "Get the params() for a saem model"

#' @export
rxUiGet.saemModel <- function(x, ...) {
  .s <- rxUiGet.loadPruneSaem(x, ...)
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
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  if (.sumProd) {
    .malert("stabilizing round off errors in saem model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "saem model")
     .msuccess("done")
  }
  paste(c(rxUiGet.saemParams(x, ...),
          .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")
}


#' @export
rxUiGet.saemInParsAndMuRefCovariates <- function(x, ...) {
  .ui <- x[[1]]
  # mu ref final removes time varying covariates
  if (exists("muRefFinal", .ui)) {
    .muRefFinal <- .ui$muRefFinal
  } else {
    .muRefFinal <- .ui$muRefCovariateDataFrame
  }
  .cov <- .ui$covariates
  .muCov <- unique(.muRefFinal$covariate)
  .cov <- .cov[!(.cov %in% .muCov)]
  if (length(.ui$predDf$cond) > 1) {
    .cov <- unique(c("CMT", .cov))
  }
  list(inPars=.cov, covars=.muCov)
}
#attr(rxUiGet.saemInParsAndMuRefCovariates, "desc") <- "Get inPars and covars for saem"

#' @export
rxUiGet.saemInPars <- function(x, ...) {
  .ret <- rxUiGet.saemInParsAndMuRefCovariates(x, ...)
  .ret$inPars
}
#attr(rxUiGet.saemInPars, "desc") <- "get inPars"

#' @export
rxUiGet.saemCovars <- function(x, ...) {
  .ret <- rxUiGet.saemInParsAndMuRefCovariates(x, ...)
  .ret$covars
}
#attr(rxUiGet.saemInPars, "desc") <- "get saemn mu-referenced non-time varying covariates"

#' @export
rxUiGet.saemFunction <- function(x, ...) {
  # This function depends on the number of time varying covariates in the data
  .ui <- x[[1]]
  .mod <- rxode2::rxode2(rxUiGet.saemModel(x, ...))
  .fnPred <- bquote(function(a, b, c) {
    rxode2::rxLoad(.(.mod))
    rxode2::rxLock(.(.mod))
    rxode2::rxAllowUnload(FALSE)
    on.exit({
      rxode2::rxUnlock(.(.mod))
      rxode2::rxAllowUnload(TRUE)
      rxode2::rxSolveFree()
    })
    .Call(`_nlmixr2_saem_do_pred`, a, b, c)
  })
  .fn <- bquote(function(a, b, c) {
    rxode2::rxLoad(.(.mod))
    rxode2::rxLock(.(.mod))
    on.exit({
      rxode2::rxUnlock(.(.mod))
      rxode2::rxAllowUnload(TRUE)
      rxode2::rxSolveFree()
    })
    if (missing(b) && missing(c)) {
      .ret <- .Call(`_nlmixr2_saem_fit`, a, PACKAGE = "nlmixr2")
      attr(.ret, "dopred") <- .(.fnPred)
      return(.ret)
    } else {
      .curFn <- .(.fnPred)
      return(.curFn(a, b, c))
    }
  })
  .inPars <- rxUiGet.saemInPars(x, ...)
  .param <- rxode2::rxParam(.mod)
  .estParam <- rxUiGet.saemParamsToEstimate(x, ...)
  .parmUpdate <- vapply(.param, function(x) {
    if (x %in% .estParam) {
      return(1L)
    } else {
      return(0L)
    }
  }, integer(1), USE.NAMES=FALSE)
  .nendpnt <- length(.ui$predDf$cond)
  .fn <- eval(.fn)
  attr(.fn, "form") <- "ode" ## Not sure this is necessary any more
  attr(.fn, "neq") <- length(rxode2::rxState(.mod))
  attr(.fn, "nlhs") <- length(rxode2::rxLhs(.mod))
  attr(.fn, "nrhs") <- sum(.parmUpdate)
  attr(.fn, "paramUpdate") <- .parmUpdate
  attr(.fn, "rx") <- .mod
  attr(.fn, "inPars") <- .inPars
  attr(.fn, "nendpnt") <- .nendpnt
  .fn
}


## message(fit$saem.rx1)
## ka=exp(tka);
## cl=exp(tcl);
## v=exp(tv);
## nlmixr_lincmt_pred=linCmtA(rx__PTR__,t,0,1,1,cl,v,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,ka,0.0,1.0,0.0,0.0);
## nlmixr_pred=nlmixr_lincmt_pred;


## message(fit$saem.rx1)
##     ka = exp(tka)
##     cl = exp(tcl)
##     v = exp(tv)
##     d/dt(depot) = -ka * depot
##     d/dt(center) = ka * depot - cl/v * center
##     cp = center/v
##     nlmixr_pred <- cp
## cmt(cp);


## message(fit$saem.rx1)
##     ktr = exp(tktr)
##     ka = exp(tka)
##     cl = exp(tcl)
##     v = exp(tv)
##     emax = expit(temax)
##     ec50 = exp(tec50)
##     kout = exp(tkout)
##     e0 = exp(te0)
##     DCP = center/v
##     PD = 1 - emax * DCP/(ec50 + DCP)
##     effect(0) = e0
##     kin = e0 * kout
##     d/dt(depot) = -ktr * depot
##     d/dt(gut) = ktr * depot - ka * gut
##     d/dt(center) = ka * gut - cl/v * center
##     d/dt(effect) = kin * PD - kout * effect
##     cp = center/v
##     if (CMT == 5) {
##         nlmixr_pred <- cp
##     }
##     if (CMT == 4) {
##         nlmixr_pred <- effect
##     }
## cmt(cp);

## cmt(effect);


## predDf has var for the variable that the prediction is defined and cmt shows the compartment that is used for this model

.saemGenModel <- function(ui, timeVaryingCovariates=character(0)) {
  .muRefCovariateDataFrame <- ui$muRefCovariateDataFrame
  if (length(timeVaryingCovariates) > 0) {
    # Drop time-varying covariates
    .muRefCovariateDataFrame <- .muRefCovariateDataFrame[!(.muRefCovariateDataFrame$covariate %in% timeVaryingCovariates), ]
  }
  assign("muRefFinal", .muRefCovariateDataFrame, ui)
  on.exit(rm(list="muRefFinal", ui))

}


#' @rdname nlmixr2Est
#' @export
nlmixr2Est.saem <- function(env, ...) {

}
