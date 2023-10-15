#' nlmixr2 optim defaults
#'
#'
#' @inheritParams stats::nlm
#' @inheritParams foceiControl
#' @inheritParams saemControl
#' @param covMethod allows selection of "r", which uses nlmixr2's
#'   `nlmixr2Hess()` for the hessian calculation or "optim" which uses
#'   the hessian from `stats::optim(.., hessian=TRUE)`
#'
#' @param trace Non-negative integer. If positive, tracing information
#'   on the progress of the optimization is produced. Higher values
#'   may produce more tracing information: for method `"L-BFGS-B"`,
#'   there are six levels of tracing. See `optim()` for more
#'   information
#'
#' @param fnscale An overall scaling to be applied to the value of `fn`
#' and `gr` during optimization. If negative, turns the problem
#' into a maximization problem. Optimization is performed on
#' `fn(par)/fnscale`
#'
#' @param parscale A vector of scaling values for the parameters.
#'   Optimization is performed on `par/parscale` and these should be
#'   comparable in the sense that a unit change in any element
#'   produces about a unit change in the scaled value.  Not used (nor
#'   needed) for `method = "Brent"`
#'
#' @param ndeps A vector of step sizes for the finite-difference
#'   approximation to the gradient, on `par/parscale` scale.  Defaults
#'   to `1e-3`
#'
#' @param maxit The maximum number of iterations. Defaults to `100`
#'   for the derivative-based methods, and `500` for `"Nelder-Mead"`.
#'
#' @param abstol The absolute convergence tolerance. Only useful for
#'    non-negative functions, as a tolerance for reaching zero.
#'
#' @param reltol Relative convergence tolerance.  The algorithm stops
#'   if it is unable to reduce the value by a factor of `reltol *
#'   (abs(val) + reltol)` at a step
#'
#' @param alpha Reflection factor for the `"Nelder-Mead"` method.
#'
#' @param beta Contraction factor for the `"Nelder-Mead"` method
#'
#' @param gamma Expansion  factor for the `"Nelder-Mead"` method
#'
#' @return optimControl object for nlmixr2
#' @export
#' @author Matthew L. Fidler
#' @examples
optimControl <- function(method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                         trace=0,
                         fnscale=1.0,
                         parscale=1.0,
                         ndeps=1e-3,
                         maxit=10000,
                         abstol=1e-8,
                         reltol=1e-8,
                         alpha=1.0,
                         beta=0.5,
                         gamma=2.0,
                         REPORT=NULL,
                         warn.1d.NelderMead=TRUE,
                         type=NULL,
                         lmm=5,
                         factr=1e7,
                         pgtol=0,
                         temp=10,
                         tmax=10,
                         rxControl=NULL,
                         optExpression=TRUE, sumProd=FALSE,
                         returnOptim=FALSE,
                         addProp = c("combined2", "combined1"),
                         calcTables=TRUE, compress=TRUE,
                         covMethod=c("r", "optim", ""),
                         adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {
  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnOptim, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)
  checkmate::assertIntegerish(trace, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(fnscale, len=1, any.missing=FALSE)
  checkmate::assertNumeric(parscale, any.missing=FALSE)
  checkmate::assertNumeric(ndeps, len=1, lower=0, an.missing=FALSE)
  checkmate::assertIntegerish(maxit, len=1, any.missing=FALSE, lower=1)
  checkmate::assertNumeric(abstol, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(reltol, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(alpha, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(beta, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(gamma, len=1, lower=0, any.missing=FALSE)
  ## alpha=1.0,
  ## beta=0.5,
  ## gamma=2.0,
  ## REPORT=NULL,
  ## warn.1d.NelderMead=TRUE,
  ## type=NULL,
  ## lmm=5,
  ## factr=1e7,
  ## pgtol=0,
  ## temp=10,
  ## tmax=10,

  covMethod <- match.arg(covMethod)
  covMethod <- match.arg(method)
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
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

  .ret <- list(method=method,
               covMethod=covMethod,
               trace=trace,
               fnscale=fnscale,
               parscale=parscale,
               ndeps=ndeps,
               maxit=maxit,
               abstol=abstol,
               reltol=reltol,
               alpha=alpha,
               beta=beta,
               gamma=gamma,
               optExpression=optExpression,
               sumProd=sumProd,
               rxControl=rxControl,
               returnOptim=returnOptim, addProp=addProp, calcTables=calcTables,
               compress=compress,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "optimControl"
  .ret
}

#' Get the optim family control
#'
#' @param env optim optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.optimFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2est::optimControl()
  }
  if (!inherits(.control, "optimControl")){
    .control <- do.call(nlmixr2est::optimControl, .control)
  }
  assign("control", .control, envir=.ui)
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.optimControl <- function(control, env) {
  assign("optimControl", control, envir=env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.optim <- function(x, ...) {
  .env <- x[[1]]
  if (exists("optimControl", .env)) {
    .control <- get("optimControl", .env)
    if (inherits(.control, "optimControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "optimControl")) return(.control)
  }
  stop("cannot find optim related control object", call.=FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.optim <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- optimControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("optimControl", .ctl)
  if (!inherits(.ctl, "optimControl")) {
    .minfo("invalid control for `est=\"optim\"`, using default")
    .ctl <- optimControl()
  } else {
    .ctl <- do.call(optimControl, .ctl)
  }
  .ctl
}

.optimEnv <- new.env(parent=emptyenv())

#' @export
rxUiGet.loadPruneOptim <- function(x, ...) {
  .loadSymengine(.nlmPrune(x), promoteLinSens = FALSE)
}

#' @export
rxUiGet.nlmRxModel <- function(x, ...) {
  .s <- rxUiGet.loadPruneOptim(x, ...)
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
    .malert("stabilizing round off errors in optim model...")
    .ret <- rxode2::rxSumProdModel(.ret)
    .msuccess("done")
  }
  if (.optExpression) {
    .ret <- rxode2::rxOptExpr(.ret, "optim model")
    .msuccess("done")
  }
  paste(c(rxUiGet.foceiParams(x, ...), rxUiGet.foceiCmtPreModel(x, ...),
          .ret, .foceiToCmtLinesAndDvid(x[[1]])), collapse="\n")
}

#' Setup the data for nlm estimation
#'
#' @param dataSav Formatted Data
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.optimFitDataSetup <- function(dataSav) {
  .dsAll <- dataSav[dataSav$EVID != 2, ] # Drop EVID=2 for estimation
  if (any(names(.dsAll) == "CENS")) {
    if (!all(.dsAll$CENS == 0)) {
      stop("'optim' does not work with censored data", call. =FALSE)
    }
  }
  .nlmEnv$data <- rxode2::etTrans(.dsAll, .nlmEnv$model)
}
