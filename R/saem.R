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
#attr(rxUiGet.saemParamsToEstimate, "desc") <- "Get the parameters to estimate"

#' @export
rxUiGet.saemThetaName <- rxUiGet.saemParamsToEstimate
#attr(rxUiGet.saemParamsToEstimate, "desc") <- "Get the parameters to estimate"

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

#' @export
rxUiGet.saemFixed <- function(x, ...) {
  .ui <- x[[1]]
  .df <- .ui$iniDf
  .dft <- .df[!is.na(.df$ntheta), ]
  .fixError <- .dft[!is.na(.dft$err), ]
  if (any(.fixError$fix)) {
    stop("Residuals cannot be fixed in SAEM.")
  }
  .dft <- .dft[is.na(.dft$err), ]
  .dft <- setNames(.dft$fix, paste(.dft$name))
  .extra <- .ui$nonMuEtas
  .extra <- setNames(rep(TRUE, length(.extra)), .ui$nonMuEtas)
  c(.dft, .extra)
}
#attr(rxUiGet.saemFixed, "desc") <- "Get the saem fixed parameters"

#' @export
rxUiGet.saemEtaTrans <- function(x, ...) {
  .ui <- x[[1]]
  .etas <- .ui$iniDf[!is.na(.ui$iniDf$neta1), ]
  .etas <- .etas$name[.etas$neta1 == .etas$neta2]
  .thetas <- rxUiGet.saemParamsToEstimate(x, ...)
  .muRefDataFrame <- .ui$muRefDataFrame
  vapply(.etas, function(eta) {
    .w <- which(eta == .muRefDataFrame$eta)
    if (length(.w) == 1L) {
      .muTheta <- .muRefDataFrame$theta[.w]
      .w <- which(.muTheta == .thetas)
      if (length(.w) == 1L) return(.w)
    }
    .w <- which(eta == .thetas)
    if (length(.w) == 1L) return(.w)
    return(NA_integer_)
  }, integer(1), USE.NAMES=FALSE)
}
#attr(rxUiGet.saemEtaTrans, "desc") <- "Get the saem eta to theta translation"
#' @export
rxUiGet.saemOmegaTrans <- function(x, ...) {
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  .o <- order(.etaTrans)
  .etaTrans2 <- .etaTrans
  .c <- 1
  for (i in .o) {
    .etaTrans2[i] <- .c
    .c <- .c + 1
  }
  .etaTrans2
}
#attr(rxUiGet.saemOmegaTrans, "desc") <- "Get the saem omega to UI omega translation"


#' @export
rxUiGet.saemModelOmega <- function(x, ...) {
  .ui <- x[[1]]
  .thetas <- rxUiGet.saemParamsToEstimate(x, ...)
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  .dm <- length(.thetas)
  .mat <- matrix(rep(0, .dm * .dm), .dm)
  .iniDf <- .ui$iniDf
  .etd <- .iniDf[which(!is.na(.iniDf$neta1)), ]
  for (i in seq_along(.etd$neta1)) {
    .mat[.etaTrans[.etd$neta1[i]], .etaTrans[.etd$neta2[i]]] <-
      .mat[.etaTrans[.etd$neta2[i]], .etaTrans[.etd$neta1[i]]] <- 1
  }
  .mat
}
#attr(rxUiGet.saemModelOmega, "desc") <- "Get the saem model omega"

#' @export
rxUiGet.saemLow <- function(x, ...) {
  .ui <- x[[1]]
  .ui$predDf$trLow
}
#attr(rxUiGet.saemLow, "desc") <- "Get the saem error transformation lower boundary"

#' @export
rxUiGet.saemHi <- function(x, ...) {
  .ui <- x[[1]]
  .ui$predDf$trHi
}
#attr(rxUiGet.saemHi, "desc") <- "Get the saem error transformation higher boundary"

#' @export
rxUiGet.saemPropT <- function(x, ...) {
  .ui <- x[[1]]
  as.integer((.ui$predDf$errTypeF=="transformed")*1L)
}
#attr(rxUiGet.saemPropT, "desc") <- "Get the saem transformation type for the function"

#' @export
rxUiGet.saemYj <- function(x, ...) {
  .ui <- x[[1]]
  as.integer(.ui$predDf$transform) - 1
}
#attr(rxUiGet.saemYj, "desc") <- "Get the saem transformation type"

#' @export
rxUiGet.saemResMod <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  vapply(seq_along(.predDf$errType),
         function(i) {
           .errType <- as.integer(.predDf$errType[i])
           .hasLambda <- !is.na(.predDf$lambda[i])
           if (.hasLambda) {
             return(.errType + 5L)
           } else {
             return(.errType)
           }
         }, integer(1), USE.NAMES=FALSE)
}
#attr(rxUiGet.saemResMod, "desc") <- "saem res.mod component"

#' @export
rxUiGet.saemResNames <- function(x, ...) {
  .ui <- x[[1]]
  .err <- .ui$iniDf
  .w <- which(sapply(.err$err, function(x) any(x == c("add", "norm", "dnorm", "dlnorm", "lnorm", "logn", "dlogn"))))
  .ret <- c()
  if (length(.w) == 1) {
    if (!is.na(.err$est[.w])) {
      .ret[length(.ret) + 1] <- paste(.err$name[.w])
    }
  }
  .w <- c(which(.err$err == "prop"), which(.err$err == "propT"))
  if (length(.w) == 1) {
    .ret[length(.ret) + 1] <- paste(.err$name[w])
  }
  return(.ret)
}
#attr(rxUiGet.saemResNames, "desc") <- "Get error names for SAEM"

rxUiGet.saemParHistNames <- function(x, ...) {
  #join_cols(join_cols(Plambda, Gamma2_phi1.diag()), vcsig2).t();
}

#' @export
rxUiGet.saemAres <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .ini <- .ui$iniDf
  .ini <- .ini[!is.na(.ini$err), ]
  return(vapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(vapply(.tmp$err, function(x) {
      x %in% c(
        "add", "norm", "dnorm", "dpois",
        "pois", "dbinom", "binom", "dbern", "bern",
        "lnorm", "dlnorm", "logn", "dlogn")
    }, logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      return(10)
    }
  }, numeric(1), USE.NAMES=FALSE))
}
#attr(rxUiGet.saemAres, "desc") <- "ares initial estimates for saem"

#' @export
rxUiGet.saemBres <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .ini <- .ui$iniDf
  .ini <- .ini[!is.na(.ini$err), ]
  return(vapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(vapply(.tmp$err,
                       function(x) (any(x == "prop") || any(x == "propT")),
                       logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      .w <- which(vapply(.tmp$err,
                         function(x) (any(x == "pow") || any(x == "powT")),
                         logical(1), USE.NAMES=FALSE))
      if (length(.w) == 1) {
        return(.tmp$est[.w])
      } else {
        return(1)
      }
    }
  }, numeric(1), USE.NAMES=FALSE))
}
#attr(rxUiGet.saemBres, "desc") <- "bres initial estimates for saem"

#' @export
rxUiGet.saemCres <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .ini <- .ui$iniDf
  .ini <- .ini[!is.na(.ini$err), ]
  return(vapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(vapply(.tmp$err, function(x) (any(x == "pow2") || any(x == "powT2")),
                       logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      return(1)
    }
  }, numeric(1), USE.NAMES=FALSE))
}
#attr(rxUiGet.saemCres, "desc") <- "cres initial estimates for saem"
#' @export
rxUiGet.saemLres <- function(x, ...) {
  .ui <- x[[1]]
 .predDf <- .ui$predDf
  .ini <- .ui$iniDf
  .ini <- .ini[!is.na(.ini$err), ]
  return(vapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .boxCox <- which(.tmp$err == "boxCox")
    if (length(.boxCox) == 1L) {
      return(.tmp$est[.boxCox])
    }
    .yeoJohnson <- which(.tmp$err == "yeoJohnson")
    if (length(.yeoJohnson) == 1L) {
      return(.tmp$est[.yeoJohnson])
    }
    return(1.0)
  }, numeric(1), USE.NAMES=FALSE))
}
#attr(rxUiGet.saemLres, "desc") <- "lres (lambda) initial estimates for saem"

#' @export
rxUiGet.saemDistribution <- function(x, ...) {
  .ui <- x[[1]]
  .df <- .ui$iniDf$err
  .df <- paste(.df[which(!is.na(.df))])
  if (any(.df %in% c("dpois", "pois"))) {
    return("poisson")
  }
  if (any(.df %in% c("dbern", "bern", "dbinom", "binom"))) {
    if (.df %in% c("dbinom", "binom")) {
      .df <- obj$ini
      .w <- which(.df$err %in% c("dbinom", "binom"))
      if (length(.w) != 1L) stop("Distribution unsupported by SAEM")
      if (!is.na(.df$name[.w])) stop("Distribution unsupported by SAEM")
    }
    return("binomial")
  }
  if (any(.df %in% c("dlnorm", "lnorm", "logn", "dlogn", "logitNorm", "probitNorm"))) {
    return("normal")
  }
  if (any(.df %in% c("dnorm", "norm", "prop", "propT", "add", "pow", "powT", "pow2", "powT2"))) {
    return("normal")
  }

  stop("Distribution unsupported by SAEM")
}


#' @export
rxUiGet.saemLogEta <- function(x, ...) {
  .ui <- x[[1]]
  .thetas <- rxUiGet.saemParamsToEstimate(x, ...)
  .ce <- .ui$muRefCurEval
  vapply(.thetas, function(x) {
    .w <- which(.ce$parameter == x)
    if (length(.w) == 1L) return(.ce$curEval[.w] == "exp")
    FALSE
  }, logical(1))
}
#attr(rxUiGet.saemLogEta, "desc") <- "saem's log.eta for saem"

#' @export
rxUiGet.saemModelList <- function(x, ...) {
  .ui <- x[[1]]
  .mod <- list(saem_mod = .ui$saemFunction)
  .covars <- rxUiGet.saemCovars(x, ...)
  if (length(.covars) > 0) {
    .mod$covars <- .covars
  }
  .mod$res.mod <- rxUiGet.saemResMod(x, ...)
  .mod$log.eta <- rxUiGet.saemLogEta(x, ...)
  .mod$ares    <- rxUiGet.saemAres(x, ...)
  .mod$bres    <- rxUiGet.saemBres(x, ...)
  .mod$omega   <- rxUiGet.saemModelOmega(x, ...)
  .mod
}
#attr(rxUiGet.saemModelList "desc") <- "saem's log.eta for saem"

#' @export
rxUiGet.saemInitTheta <- function(x, ...) {
  .logEta <- rxUiGet.saemLogEta(x, ...)
  .names <- names(.logEta)
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .est <- .iniDf[!is.na(.iniDf$ntheta) & is.na(.iniDf$err), "est"]
  .etaNames <- .iniDf[is.na(.iniDf$ntheta), ]
  .etaNames <- .iniDf[.iniDf$neta1 == .iniDf$neta2, "name"]
  .fixed <- rxUiGet.saemFixed(x, ...)
  .n <- vapply(.fixed, function(x) ifelse(x, "FIXED", ""),
               character(1), USE.NAMES=FALSE)

  setNames(vapply(seq_along(.logEta),
                  function(i){
                    .isEta <- any(.names[i] %in% .etaNames)
                    if (.logEta[i]) {
                      if (.isEta) {
                        return(1)
                      } else {
                        return(exp(.est[i]))
                      }
                    } else {
                      if (.isEta) {
                        return(0)
                      } else {
                        return(.est[i])
                      }
                    }
                  }, numeric(1), USE.NAMES=FALSE), .n)
}
#attr(rxUiGet.saemInitTheta, "desc") <- "initialization for saem's theta"

#' @export
rxUiGet.saemInitOmega <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .eta <- .iniDf[is.na(.iniDf$ntheta), ]
  .eta <- .eta[.eta$neta1 == .eta$neta2, ]
  .eta <- setNames(.eta$est, .eta$name)
  .pars <- rxUiGet.saemParamsToEstimate(x, ...)
  .ret <- rep(1.0, length(.pars))
  .etaTrans <- rxUiGet.saemEtaTrans(x, ...)
  for (i in seq_along(.etaTrans)) {
    .ret[.etaTrans[i]] <- .eta[i]
  }
  setNames(.ret, .pars)
}
#attr(rxUiGet.saemInitOmega, "desc") <- "initialization for saem's omega"

#' @export
rxUiGet.saemInit <- function(x, ...) {
  list(theta=rxUiGet.saemInitTheta(x, ...),
       omega=rxUiGet.saemInitOmega(x, ...))
}
#attr(rxUiGet.saemInit, "desc") <- "initialization for saem's theta and omega"

#' @export
rxUiGet.saemResName <- function(x, ...) {
  .ui <- x[[1]]
  .w <- which(vapply(.ui$iniDf$err,
                     function(x) any(x == c("add", "norm", "dnorm", "dlnorm", "lnorm", "logn", "dlogn")),
                     logical(1),
                     USE.NAMES=FALSE))
  .ret <- NULL
  if (length(.w) == 1) {
    if (!is.na(.ui$iniDf$est[.w])) {
      .ret[length(.ret) + 1] <- paste(.ui$iniDf$name[.w])
    }
  }
  .w <- c(which(.ui$iniDf$err == "prop"), which(.ui$iniDf$err == "propT"))
  if (length(.w) == 1) {
    .ret[length(.ret) + 1] <- paste(obj$name[.w])
  }
  .ret
}
#attr(rxUiGet.saemResName, "desc") <- "Residual Names for saem"

#' Fit a UI model with saem
#'
#' @param ui rxode2 ui
#' @param data nlmixr data
#' @param timeVaryingCovariates Time varying covarites in the data
#' @return lower level saem fit
#' @author Matthew L. Fidler
#' @noRd
.saemFitModel <- function(ui, data, timeVaryingCovariates=character(0)) {
  .muRefCovariateDataFrame <- ui$muRefCovariateDataFrame
  if (length(timeVaryingCovariates) > 0) {
    # Drop time-varying covariates
    .muRefCovariateDataFrame <- .muRefCovariateDataFrame[!(.muRefCovariateDataFrame$covariate %in% timeVaryingCovariates), ]
  }
  assign("muRefFinal", .muRefCovariateDataFrame, ui)
  on.exit(rm(list="muRefFinal", envir=ui))
  .model <- ui$saemModelList
  .inits <- ui$saemInit
  .cfg <- .configsaem(model=.model,
                      data=data,
                      inits=.inits,
                      mcmc=rxode2::rxGetControl(ui, "mcmc", list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2))),
                      ODEopt=rxode2::rxGetControl(ui, "ODEopt", rxode2::rxControl()),
                      distribution="normal",
                      addProp="combined2",
                      fixed=.fixed,
                      seed=rxode2::rxGetControl(ui, "seed", 99),
                      DEBUG=rxode2::rxGetControl(ui, "DEBUG", 0),
                      tol=rxode2::rxGetControl(ui, "tol", 1e-6),
                      itmax=rxode2::rxGetControl(ui, "itmax", 30),
                      type=rxode2::rxGetControl(ui, "type", "nelder-mead"),
                      lambdaRange=rxode2::rxGetControl(ui, "lambdaRange", 3),
                      powRange=rxode2::rxGetControl(ui, "powRange", 10),
                      odeRecalcFactor=rxode2::rxGetControl(ui, "odeRecalcFactor", 10^0.5),
                      maxOdeRecalc=rxode2::rxGetControl(ui, "maxOdeRecalc", 10^0.5))
  .print <- rxode2::rxGetControl(ui, "print", 1)
  if (inherits(.print, "numeric")) {
    .cfg$print <- as.integer(.print)
  }
  .cfg$cres <- ui$saemCres
  .cfg$yj <- ui$saemYj
  .cfg$lres <- ui$saemLres
  .cfg$low <- ui$saemLow
  .cfg$hi <- ui$saemHi
  .cfg$propT <- ui$saemPropT
  .model$saem_mod(.cfg)
}

.saemFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- saemControl()
  }
  if (!inherits(.control, "saemControl")){
    .control <- do.call(nlmixr2::saemControl, .control)
  }
  assign("control", .control, envir=.ui)
}
#' Get SAEM theta
#'
#' @param env Environment that has ui and saem in it
#' @return Nothining, environment is assigned a theta
#' @author Matthew L. Fidler
#' @noRd
.getSaemTheta <- function(env) {
  .ui <- env$ui
  .saem <- env$saem
  .iniDf <- .ui$iniDf
  .predDf <- .ui$predDf
  .thetaNames <- .iniDf[!is.na(.iniDf$ntheta), "name"]
  .theta <- setNames(rep(NA_real_, length(.thetaNames)), .thetaNames)
  .saemThetaNames <- .ui$saemParamsToEstimate
  .thetaSaem <- setNames(as.vector(fixed.effects(env$saem)), .saemThetaNames)
  .resMat <- .saem$resMat
  for (n in .thetaNames) {
    if (n %in% .saemThetaNames) {
      .theta[n] <- .thetaSaem[n]
    }
  }
  for (i in seq_along(.predDf$cond)) {
    .x <- paste(.predDf$cond[i])
    .tmp <- .iniDf[which(.iniDf$condition == .x), ]
    .w <- which(vapply(.tmp$err,
                       function(x) any(x == c("prop", "propT", "pow", "powT")),
                       logical(1),
                       USE.NAMES=FALSE))
    if (length(.w) == 1) {
      .theta[paste(.tmp$name[.w])] <- .resMat[i, 2]
    }
    .w <- which(vapply(.tmp$err,
                       function(x) any(x == c("pow2", "powT2")),
                       logical(1),
                       USE.NAMES=FALSE))
    if (length(.w) == 1) {
      .theta[paste(.tmp$name[.w])] <- .resMat[i, 3]
    }
    .w <- which(vapply(seq_along(.tmp$err),
                       function(x) {
                         .x <- .tmp$err[x]
                         if (any(.x == c(
                           "add", "norm", "dnorm", "lnorm", "dlnorm",
                           "dlogn", "logn"))) {
                           if (!is.na(.tmp$est[x])) {
                             return(TRUE)
                           }
                         }
                         return(FALSE)
                       },
                       logical(1),
                       USE.NAMES=FALSE))
    if (length(.w) == 1) {
      .theta[paste(.tmp$name[.w])] <- .resMat[i, 1]
    }
    .w <- which(vapply(.tmp$err, function(x) {
      any(x == c("boxCox", "yeoJohnson"))
    }, logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1) {
      .theta[paste(.tmp$name[.w])] <- .resMat[i, 4]
    }
  }
  env$theta <- .theta
  invisible()
}

#' Get SAEM omega
#'
#' @param env Environment that has ui and saem in it
#' @return Nothing, environment is assigned the omega
#' @author Matthew L. Fidler
#' @noRd
.getSaemOmega <- function(env) {
  ## Reorder based on translation
  .saem <- env$saem
  .ui <- env$ui
  .etaTrans <- .ui$saemOmegaTrans
  ## saem eta ->  ui eta
  .df <- .ui$iniDf
  .eta <- .df[!is.na(.df$neta1), ]
  .etaNames <- .eta[.eta$neta1 == .eta$neta2, "name"]
  .neta <- length(.etaNames)
  .len <- length(.etaNames)
  .ome <- matrix(rep(0, .len * .len), .len, .len, dimnames=list(.etaNames, .etaNames))
  .curOme <- .saem$Gamma2_phi1
  for (i in seq_along(.eta$name)) {
    .e1 <- .eta$neta1[i]
    .e2 <- .eta$neta2[i]
    .o1 <- .etaTrans[.e1]
    .o2 <- .etaTrans[.e2]
    .ome[.e1, .e2] <- .curOme[.o1, .o2]
    .ome[.e2, .e1] <- .curOme[.o2, .o1]
  }
  env$omega <- .ome
  invisible()
}

.saemAddParHist <- function(env) {
  .saem <- env$saem
  .ui <- env$ui

  .m <- .saem$par_hist
  if (ncol(.m) > length(.allThetaNames)) {
    .m <- .m[, seq_along(.allThetaNames)]
  }
  env$parHistStacked <- data.frame(
    val = as.vector(.m),
    par = rep(.allThetaNames, each = nrow(.m)),
    iter = rep(1:nrow(.m), ncol(.m))
  )
  env$parHist <- data.frame(iter = rep(1:nrow(.m)), as.data.frame(.m))
}

.saemCalcCov <- function(env) {
  .ui <- env$ui
  .saem <- env$saem
  .covMethod <- rxode2::rxGetControl(.ui, "covMethod", "linFim")
  .calcCov <- .covMethod == "linFim"
  if (.covMethod == "") {
    .cov <- NULL
    .addCov <- FALSE
  } else {
    .tn <- .ui$saemParamsToEstimate[!.ui$saemFixed]
    .nth <- length(.tn)

    .ini <- .ui$iniDf
    .ini <- .ini[is.na(.ini$err), ]
    .ini <- .ini[!is.na(.ini$ntheta), ]
    .ini <- .ini[!.ini$fix, ]
    .ini <- paste(.ini$name)
    .calcCovTime <- proc.time()
    if (.calcCov) {
      .covm <- .saem$Ha[1:.nth, 1:.nth]
      .covm <- try(calc.COV(.saem))
      .doIt <- !inherits(.covm, "try-error")
      if (.doIt && dim(.covm)[1] != .nth) .doIt <- FALSE
      if (.doIt) {
        .tmp <- try(chol(.covm), silent = TRUE)
        .addCov <- TRUE
        .sqrtm <- FALSE
        if (inherits(.tmp, "try-error")) {
          .tmp <- .covm
          .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
          if (inherits(.tmp, "try-error")) {
            .calcCov <- FALSE
            .covm <- .saem$Ha[1:.nth, 1:.nth]
            .tmp <- try(chol(.covm), silent = TRUE)
            .addCov <- TRUE
            .sqrtm <- FALSE
            if (inherits(.tmp, "try-error")) {
              .tmp <- .saem$Ha[1:.nth, 1:.nth]
              .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
              if (inherits(.tmp, "try-error")) {
                .addCov <- FALSE
              } else {
                .sqrtm <- TRUE
              }
            } else {
              .tmp <- .saem$Ha[1:.nth, 1:.nth]
            }
          } else {
            .sqrtm <- TRUE
          }
        } else {
          .tmp <- .covm
        }
      } else {
        .tmp <- .saem$Ha[1:.nth, 1:.nth]
        .tmp <- try(chol(.covm), silent = TRUE)
        .calcCov <- FALSE
        .addCov <- TRUE
        .sqrtm <- FALSE
        if (inherits(.tmp, "try-error")) {
          .tmp <- .saem$Ha[1:.nth, 1:.nth]
          .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
          if (inherits(.tmp, "try-error")) {
            .addCov <- FALSE
          } else {
            .sqrtm <- TRUE
          }
        } else {
          .tmp <- .saem$Ha[1:.nth, 1:.nth]
          .calcCov <- FALSE
        }
      }
    } else {
      .tmp <- try(chol(.covm), silent = TRUE)
      .addCov <- TRUE
      .sqrtm <- FALSE
      if (inherits(.tmp, "try-error")) {
        .tmp <- .saem$Ha[1:.nth, 1:.nth]
        .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
        if (inherits(.tmp, "try-error")) {
          .addCov <- FALSE
        } else {
          .sqrtm <- TRUE
        }
      } else {
        .tmp <- .saem$Ha[1:.nth, 1:.nth]
        .calcCov <- FALSE
      }
    }
    if (.addCov) {
      if (!.calcCov) {
        .cov <- rxode2::rxInv(.tmp)
      } else {
        .cov <- .tmp
      }
      attr(.cov, "dimnames") <- list(.tn, .tn)
      .cov <- .cov[.ini, .ini, drop = FALSE]
    }
    .calcCovTime <- proc.time() - .calcCovTime
    .calcCovTime <- .calcCovTime["elapsed"]
  }
  if (.addCov) {
    env$cov <- .cov
    env$.calcCovTime <- .calcCovTime
  }
}

#' Fit the saem family of models
#'
#' @param env Environment from nlmixr2Est
#' @param ... Other arguments
#' @return fit environment with $saem, $saemControl, $dataSav, $origData, $ui
#' @author Matthew L. Fidler
#' @noRd
.saemFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  .foceiPreProcessData(.data, .ret, .ui)
  .et <- rxode2::etTrans(.ret$dataSav, .ui$mv0)
  .nTv <- attr(class(.et), ".rxode2.lst")$nTv
  if (is.null(.nTv)) .nTv <- 0
  .tv <- character(0)
  if (.nTv != 0) {
    .tv <- names(.et)[-seq(1, 6)]
  }
  .ret$saem <- .saemFitModel(.ui, .ret$dataSav, timeVaryingCovariates=.tv)
  .ret$saemControl <- .control
  .ret$ui <- .ui
  .saemCalcCov(.ret)
  .getSaemTheta(.ret)
  .getSaemOmega(.ret)
  .nlmixr2FitUpdateParams(.ret)
   if (exists("control", .ui)) {
    rm(list="control", envir=.ui)
  }
  .ret
}


#' @rdname nlmixr2Est
#' @export
nlmixr2Est.saem <- function(env, ...) {
  .saemFamilyControl(env, ...)
  .saemFamilyFit(env,  ...)
}
