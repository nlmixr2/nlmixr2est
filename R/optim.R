#' nlmixr2 optim defaults
#'
#'
#' @inheritParams stats::optim
#' @inheritParams foceiControl
#' @inheritParams saemControl
#'
#' @param solveType tells if `optim` will use nlmixr2's analytical
#'   gradients when available (finite differences will be used for
#'   event-related parameters like parameters controlling lag time,
#'   duration/rate of infusion, and modeled bioavailability). This can
#'   be:
#'
#' - `"gradient"` which will use the gradient and let `nlm` calculate
#'    the finite difference hessian
#'
#' - `"fun"` where nlm will calculate both the finite difference
#'    gradient and the finite difference Hessian
#'
#'  When using nlmixr2's finite differences, the "ideal" step size for
#'  either central or forward differences are optimized for with the
#'  Shi2021 method which may give more accurate derivatives
#'
#' These are only applied in the gradient based methods: "BFGS", "CG",
#' "L-BFGS-B"
#'
#' @param returnOptim logical; when TRUE this will return the optim
#'   list instead of the nlmixr2 fit object
#'
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
#' @param REPORT The frequency of reports for the `"BFGS"`,
#'   `"L-BFGS-B"` and `"SANN"` methods if `control$trace` is
#'   positive. Defaults to every 10 iterations for `"BFGS"` and
#'   `"L-BFGS-B"`, or every 100 temperatures for `"SANN"`
#'
#' @param warn.1d.NelderMead a logical indicating if the (default)
#'   `"Nelder-Mead"` method should signal a warning when used for
#'   one-dimensional minimization.  As the warning is sometimes
#'   inappropriate, you can suppress it by setting this option to
#'   `FALSE`
#'
#' @param type for the conjugate-gradients method.  Takes value `1`
#'   for the Fletcher-Reeves update, `2` for Polak-Ribiere and `3` for
#'   Beale-Sorenson.
#'
#' @param lmm is an integer giving the number of BFGS updates retained
#'   in the `"L-BFGS-B"` method, It defaults to `5`
#'
#' @param factr controls the convergence of the `"L-BFGS-B"` method.
#'   Convergence occurs when the reduction in the objective is within
#'   this factor of the machine tolerance. Default is `1e7`, that is a
#'   tolerance of about `1e-8`.
#'
#' @param pgtol helps control the convergence of the ‘"L-BFGS-B"’
#'   method.  It is a tolerance on the projected gradient in the
#'   current search direction. This defaults to zero, when the check
#'   is suppressed
#'
#' @param temp controls the `"SANN"` method. It is the starting
#'   temperature for the cooling schedule. Defaults to `10`.
#'
#' @param tmax is the number of function evaluations at each
#'   temperature for the `"SANN"` method. Defaults to `10`.
#'
#' @return optimControl object for nlmixr2
#' @export
#' @author Matthew L. Fidler
#' @examples
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
#' fit2 <- nlmixr(mod, dsn, est="optim", optimControl(method="BFGS"))
#' fit2
#' }
optimControl <- function(method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                         trace=10,
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
                         stickyRecalcN=4,
                         maxOdeRecalc=5,
                         odeRecalcFactor=10^(0.5),
                         eventType=c("central", "forward"),
                         shiErr=(.Machine$double.eps)^(1/3),
                         shi21maxFD=20L,
                         solveType=c("grad", "fun"),
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
  checkmate::assertIntegerish(trace, len=1, any.missing=FALSE, lower=0) #nolint
  checkmate::assertNumeric(fnscale, len=1, any.missing=FALSE)
  checkmate::assertNumeric(parscale, any.missing=FALSE)
  checkmate::assertNumeric(ndeps, lower=0, any.missing=FALSE)
  checkmate::assertIntegerish(maxit, len=1, any.missing=FALSE, lower=1)
  checkmate::assertNumeric(abstol, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(reltol, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(alpha, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(beta, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(gamma, len=1, lower=0, any.missing=FALSE)
  checkmate::assertIntegerish(REPORT, len=1, lower=0, any.missing=FALSE, null.ok=TRUE)
  checkmate::assertLogical(warn.1d.NelderMead, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(type, len=1, lower=1, upper=3, any.missing=FALSE, null.ok=TRUE)
  checkmate::assertIntegerish(lmm, len=1, lower=1, any.missing=FALSE)
  checkmate::assertNumeric(factr, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(pgtol, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(temp, len=1, lower=0, any.missing=FALSE)
  checkmate::assertIntegerish(tmax, len=1, lower=0, any.missing=FALSE)

  .solveTypeIdx <- c("hessian" = 3L, "grad" = 2L, "fun" = 1L)
  if (checkmate::testIntegerish(solveType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    solveType <- as.integer(solveType)
  } else {
    solveType <- setNames(.solveTypeIdx[match.arg(solveType)], NULL)
  }
  method <- match.arg(method)
  if (missing(covMethod) && any(solveType == 2:3) &&
        method %in% c("BFGS", "CG", "L-BFGS-B")) {
    covMethod <- "optim"
  } else {
    covMethod <- match.arg(covMethod)
  }

  .eventTypeIdx <- c("central" =2L, "forward"=1L)
  if (checkmate::testIntegerish(eventType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    eventType <- as.integer(eventType)
  } else {
    eventType <- setNames(.eventTypeIdx[match.arg(eventType)], NULL)
  }

  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)




  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% "genRxControl")]
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
               trace=trace, # nolint
               fnscale=fnscale,
               parscale=parscale,
               ndeps=ndeps,
               maxit=maxit,
               abstol=abstol,
               reltol=reltol,
               alpha=alpha,
               beta=beta,
               gamma=gamma,
               REPORT=REPORT,
               warn.1d.NelderMead=warn.1d.NelderMead,
               type=type,
               lmm=lmm,
               factr=factr,
               pgtol=pgtol,
               temp=temp,
               tmax=tmax,
               optExpression=optExpression,
               sumProd=sumProd,
               solveType=solveType,
               stickyRecalcN=as.integer(stickyRecalcN),
               maxOdeRecalc=as.integer(maxOdeRecalc),
               odeRecalcFactor=odeRecalcFactor,
               eventType=eventType,
               shiErr=shiErr,
               shi21maxFD=as.integer(shi21maxFD),
               rxControl=rxControl,
               returnOptim=returnOptim, addProp=addProp, calcTables=calcTables,
               compress=compress,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "optimControl"
  .ret
}


#' A surrogate function for optim to call for ode solving
#'
#' @param pars Parameters that will be estimated
#' @return Predictions
#' @details
#' This is an internal function and should not be called directly.
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.nlmixrOptimFunC <- function(pars) {
  .Call(`_nlmixr2est_optimFunC`, pars, FALSE)
}
#' @rdname dot-nlmixrOptimFunC
#' @export
.nlmixrOptimGradC <- function(pars) {
  .Call(`_nlmixr2est_optimFunC`, pars, TRUE)
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

#' Setup the data for optim estimation
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

#' @export
rxUiGet.optimParLower <- function(x, ...) {
  .ui <- x[[1]]
  .ui$iniDf$lower[!.ui$iniDf$fix]
}

#' @export
rxUiGet.optimParUpper <- function(x, ...) {
  .ui <- x[[1]]
  .ui$iniDf$upper[!.ui$iniDf$fix]
}

.optimFitModel <- function(ui, dataSav) {
  # Use nlmEnv and function for DRY principle
  .ctl <- ui$control
  .keep <- c("trace", "fnscale", "parscale", "ndeps", "maxit",
             "abstol", "reltol", "alpha", "beta", "gamma",
             "REPORT", "warn.1d.NelderMead", "type", "lmm",
             "factr", "pgtol", "temp", "tmax")
  if (is.null(.ctl$REPORT)) {
    .keep <- .keep[.keep != "REPORT"]
  }
  if (is.null(.ctl$type)) {
    .keep <- .keep[.keep != "type"]
  }
  .oCtl <- setNames(lapply(.keep, function(x) {.ctl[[x]]}), .keep)
  class(.ctl) <- NULL
  .par <- ui$nlmParIni
  if (length(.oCtl$parscale) == 1L) {
    .oCtl$parscale <- rep(.oCtl$parscale, length(.par))
  } else if (!(length(.oCtl$parscale) == length(.par))) {
    stop("'parscale' should match the number of parmeters (currently ",
         length(.oCtl$parscale), " should be ", length(.par), ")",
         call.=FALSE)
  }
  if (length(.oCtl$ndeps) == 1L) {
    .oCtl$ndeps <- rep(.oCtl$ndeps, length(.par))
  } else if (!(length(.oCtl$ndeps) == length(.par))) {
    stop("'ndeps' should match the number of parmeters (currently ",
         length(.oCtl$ndeps), " should be ", length(.par), ")",
         call.=FALSE)
  }

  # Fill in options for hessian which isn't supported in this method
  .ctl$optimHessType <- 2L
  .ctl$hessErr <- (.Machine$double.eps)^(1/3)
  .ctl$shi21maxHess <- 20L

  if (.ctl$method %in% c("BFGS", "CG", "L-BFGS-B") &&
        .ctl$solveType == 2L) {
    # support gradient
    .f <- ui$nlmSensModel
    .env <- new.env(parent=emptyenv())
    .env$rxControl <- .ctl$rxControl
    .env$predOnly <- .f$predOnly
    .nlmEnv$model <- .env$thetaGrad <- .f$thetaGrad
    .env$param <- setNames(.par, sprintf("THETA[%d]", seq_along(.par)))
    .nlmFitDataSetup(dataSav)
    .env$needFD <- .f$eventTheta
    .env$control <- .ctl
    .env$data <- .nlmEnv$data
    .Call(`_nlmixr2est_nlmSetup`, .env)
    on.exit({
      .Call(`_nlmixr2est_nlmFree`)
      rxode2::rxSolveFree()
    })
    .ret <- bquote(stats::optim(
      par=.(.par),
      fn=.(nlmixr2est::.nlmixrOptimFunC),
      gr=.(nlmixr2est::.nlmixrOptimGradC),
      method=.(.ctl$method),
      control=.(.oCtl),
      lower=.(ui$optimParLower),
      upper=.(ui$optimParUpper),
      hessian=.(.ctl$covMethod == "optim")))
    .ret <- eval(.ret)
  } else {
    # don't support gradient
    .f <- ui$nlmRxModel
    .nlmEnv$model <- .predOnly <- .f$predOnly
    .ctl$solveType <- 1L
    .env <- new.env(parent=emptyenv())
    .env$rxControl <- .ctl$rxControl
    .env$predOnly <- .predOnly
    .env$param <- setNames(.par, sprintf("THETA[%d]", seq_along(.par)))
    .nlmFitDataSetup(dataSav)
    .env$needFD <- .f$eventTheta
    .env$control <- .ctl
    .env$data <- .nlmEnv$data
    .Call(`_nlmixr2est_nlmSetup`, .env)
    on.exit({
      .Call(`_nlmixr2est_nlmFree`)
      rxode2::rxSolveFree()
    })
    .ret <- bquote(stats::optim(
      par=.(.par),
      fn=.(nlmixr2est::.nlmixrOptimFunC),
      method=.(.ctl$method),
      control=.(.oCtl),
      lower=.(ui$optimParLower),
      upper=.(ui$optimParUpper),
      hessian=.(.ctl$covMethod == "optim")))
    .ret <- eval(.ret)
  }
  # be nice and name items
  .name <- ui$nlmParName
  names(.ret$par) <- .name
  if (any(.ctl$covMethod == c("r", "optim"))) {
    .malert("calculating covariance")
    if (.ctl$covMethod != "nlm") {
      .p <- setNames(.ret$par, NULL)
      .hess <- nlmixr2Hess(.p, nlmixr2est::.nlmixrOptimFunC)
      .ret$hessian <- .hess
    }
    dimnames(.ret$hessian) <- list(.name, .name)
    .hess <- .ret$hessian

    # r matrix
    .r <- 0.5 * .hess
    .ch <- try(cholSE(.r), silent = TRUE)
    .covType <- "r"
    if (inherits(.ch, "try-error")) {
      .r2 <- .r %*% .r
      .r2 <- try(sqrtm(.r2), silent=TRUE)
      .covType <- "|r|"
      if (!inherits(.r2, "try-error")) {
        .ch <- try(cholSE(.r), silent=TRUE)
        if (inherits(.ch, "try-error")) {
          .r2 <- .ch # switch to nearPD
        }
      }
      if (inherits(.r2, "try-error")) {
        .covType <- "r+"
        .r2 <- try(nmNearPD(.r), silent=TRUE)
        if (!inherits(.r2, "try-error")) {
          .ch <- try(cholSE(.r), silent=TRUE)
        }
      } else {
        .ch <- try(cholSE(.r), silent=TRUE)
      }
    }
    if (!inherits(.ch, "try-error")) {
      .rinv <- rxode2::rxInv(.ch)
      .rinv <- .rinv %*% t(.rinv)
      .cov <- 2*.rinv
      dimnames(.cov) <- list(.name, .name)
      dimnames(.rinv) <- list(.name, .name)
      .ret$covMethod <- .covType
      if (.ctl$covMethod == "optim") {
        .ret$covMethod <- paste0(.covType, " (optim)")
      } else {
        .ret$covMethod <- .covType
      }
      .ret$cov <- .cov

    } else {
      .ret$covMethod <- "failed"
    }
    dimnames(.r) <- list(.name, .name)
    .ret$r <- .r
    .msuccess("done")
  }
  .ret
}

#' Get the full theta for nlm methods
#'
#' @param optim enhanced nlm return
#' @param ui ui object
#' @return named theta matrix
#' @author Matthew L. Fidler
#' @noRd
.optimGetTheta <- function(nlm, ui) {
  .iniDf <- ui$iniDf
  setNames(vapply(seq_along(.iniDf$name),
                  function(i) {
                    if (.iniDf$fix[i]) {
                      return(.iniDf$est[i])
                    } else {
                      return(nlm$par[.iniDf$name[i]])
                    }
                  }, double(1), USE.NAMES=FALSE),
           .iniDf$name)
}

.optimControlToFoceiControl <- function(env, assign=TRUE) {
  .optimControl <- env$optimControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$optimControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.optimControl$sumProd,
                                optExpression=.optimControl$optExpression,
                                scaleTo=0,
                                calcTables=.optimControl$calcTables,
                                addProp=.optimControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.optimControl$compress,
                                ci=.optimControl$ci,
                                sigdigTable=.optimControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

.optimFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  # The environment needs:
  # - table for table options
  # - $origData -- Original Data
  # - $dataSav -- Processed data from .foceiPreProcessData
  # - $idLvl -- Level information for ID factor added
  # - $covLvl -- Level information for items to convert to factor
  # - $ui for ui fullTheta Full theta information
  # - $etaObf data frame with ID, etas and OBJI
  # - $cov For covariance
  # - $covMethod for the method of calculating the covariance
  # - $adjObf Should the objective function value be adjusted
  # - $objective objective function value
  # - $extra Extra print information
  # - $method Estimation method (for printing)
  # - $omega Omega matrix
  # - $theta Is a theta data frame
  # - $model a list of model information for table generation.  Needs a `predOnly` model
  # - $message Message for display
  # - $est estimation method
  # - $ofvType (optional) tells the type of ofv is currently being used
  # When running the focei problem to create the nlmixr object, you also need a
  #  foceiControl object
  .ret$table <- env$table
  .foceiPreProcessData(.data, .ret, .ui, .control$rxControl)
  .optim <- .collectWarn(.optimFitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$optim <- .optim[[1]]
  .ret$message <- .ret$optim$message
  if (rxode2::rxGetControl(.ui, "returnNlm", FALSE)) {
    return(.ret$nlm)
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .optimGetTheta(.ret$optim, .ui)
  .ret$cov <- .ret$optim$cov
  .ret$covMethod <- .ret$optim$covMethod
  #.ret$etaMat <- NULL
  #.ret$etaObf <- NULL
  #.ret$omega <- NULL
  .ret$control <- .control
  .ret$extra <- paste0(" with ", crayon::bold$yellow(.control$method),  " method")
  .nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui)) {
    rm(list="control", envir=.ui)
  }
  .ret$est <- "optim"
  # There is no parameter history for nlme
  .ret$objective <- 2 * as.numeric(.ret$optim$value)
  .ret$model <- .ui$ebe
  .ret$ofvType <- "optim"
  .optimControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="optim")
  .env <- .ret$env
  .env$method <- "optim"
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.optim <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'optim', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'optim'", .var.name=.ui$modelName)
  .optimFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .optimFamilyFit(env,  ...)
}
