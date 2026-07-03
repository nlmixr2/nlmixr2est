.regFloat1 <- rex::rex(
  or(
    group(some_of("0":"9"), ".", any_of("0":"9")),
    group(any_of("0":"9"), ".", some_of("0":"9"))
  ),
  maybe(group(one_of("E", "e"), maybe(one_of("+", "-")), some_of("0":"9")))
)
.regFloat2 <- rex::rex(some_of("0":"9"), one_of("E", "e"), maybe(one_of("-", "+")), some_of("0":"9"))
.regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
.regNum <- rex::rex(maybe("-"), or(.regDecimalint, .regFloat1, .regFloat2))


use.utf <- function() {
  opt <- getOption("cli.unicode", NULL)
  if (!is.null(opt)) {
    isTRUE(opt)
  } else {
    l10n_info()$`UTF-8` && !is.latex()
  }
 }

is.latex <- function() {
  if (!("knitr" %in% loadedNamespaces())) {
    return(FALSE)
  }
  get("is_latex_output", asNamespace("knitr"))()
}

#' Get the maxfun control for minqa optimizers
#'
#' @param control control to update based on foceiControl()
#' @return control with maxfun updated based on maxOuterIterations
#' @noRd
#' @author Matthew L. Fidler
.controlMaxfun <- function(control) {
  if (!is.null(control$maxOuterIterations)) {
    control$maxfun <- control$maxOuterIterations
  }
  control
}

.uobyqa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- .controlMaxfun(control)
  if (is.null(.ctl$npt)) .ctl$npt <- length(par) * 2 + 1
  .ctl$iprint <- 0L
  .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")]
  .ret <- minqa::uobyqa(par, fn,
                        control = .ctl,
                        lower = lower,
                        upper = upper)
  .ret$x <- .ret$par
  .ret$message <- .ret$msg
  .ret$convergence <- .ret$ierr
  .ret$value <- .ret$fval
  .ret
}

.newuoa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- .controlMaxfun(control)
  if (is.null(.ctl$npt)) .ctl$npt <- length(par) * 2 + 1
  .ctl$iprint <- 0L
  .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")]
  .ret <- minqa::newuoa(par, fn,
                        control = .ctl,
                        lower = lower,
                        upper = upper)
  .ret$x <- .ret$par
  .ret$message <- .ret$msg
  .ret$convergence <- .ret$ierr
  .ret$value <- .ret$fval
  .ret
}

.bobyqa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- .controlMaxfun(control)
  if (is.null(.ctl$npt)) .ctl$npt <- length(par) * 2 + 1
  .ctl$iprint <- 0L
  .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")]
  .ret <- minqa::bobyqa(par, fn,
                        control = .ctl,
                        lower = lower,
                        upper = upper
                        )
  .ret$x <- .ret$par
  .ret$message <- .ret$msg
  .ret$convergence <- .ret$ierr
  .ret$value <- .ret$fval
  .ret
}

#' Get the maxit control
#'
#' @param control control to update based on foceiControl()
#' @return control with maxfun updated based on maxOuterIterations
#' @noRd
#' @author Matthew L. Fidler
.controlMaxit <- function(control) {
  if (!is.null(control$maxOuterIterations)) {
    control$maxit <- control$maxOuterIterations
  }
  control
}

.lbfgsb3c <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  control <- .controlMaxit(control)
  .w <- which(names(control) %in% c("trace", "factr", "pgtol", "abstol", "reltol", "lmm", "maxit", "iprint"))
  .control <- control[.w]
  .ret <- lbfgsb3c::lbfgsb3c(par = as.vector(par), fn = fn, gr = gr, lower = lower, upper = upper, control = .control)
  .ret$x <- .ret$par
  .ret
}


.lbfgsbO <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  control <- .controlMaxit(control)
  .control <- control[names(control) %in% c("trace", "factr", "pgtol", "abstol", "reltol", "lmm", "maxit", "iprint")]
  .w <- which(sapply(.control, is.null))
  .control <- .control[-.w]
  .ret <- optim(
    par = par, fn = fn, gr = gr, method = "L-BFGS-B",
    lower = lower, upper = upper,
    control = .control, hessian = FALSE
  )
  .ret$x <- .ret$par
  .ret
}

.optimize <- function(par, fn, gr, lower=-Inf, upper=Inf, control=list(), ...) {
  # focei assumes par is the initial estimate (ignored in Brent's method)
  # fn  is the function to calculate the objective function
  # gr  is the function to calculate the gradient (ignored)
  # lower is the lower bound, in this case it must be length 1
  # upper is the upper bound, in this case it must be length 1
  .lower <- rxode2::expit(lower)
  if (is.na(.lower)) .lower <- 0.0
  .upper <- rxode2::expit(upper)
  if (is.na(.upper)) .upper <- 1.0
  f <- function(x) {
    fn(rxode2::logit(x))
  }
  .ret <- stats::optimize(f, c(.lower, .upper), tol=control$abstol)
  f <- fn
  .range <- rxode2::logit(.ret$minimum) + c(-4,4)*control$abstol
  .range[1] <- max(lower, .range[1])
  .range[2] <- min(upper, .range[2])
  .ret <- stats::optimize(f, .range, tol=control$abstol)
  .ret$x <- .ret$minimum
  .ret$message <- "stats::optimize for 1 dimensional optimization"
  .ret$convergence <- 0L
  .ret$value <- .ret$objective
  .ret
}


#' Get the maxit control
#'
#' @param control control to update based on foceiControl()
#' @return control with iter.max updated based on maxOuterIterations
#' @noRd
#' @author Matthew L. Fidler
.controlIterMax <- function(control) {
  if (!is.null(control$maxOuterIterations)) {
    control$iter.max <- control$maxOuterIterations
  }
  control
}

.nlminb <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- .controlIterMax(control)
  .ctl <- .ctl[names(.ctl) %in% c(
    "eval.max", "iter.max", "trace", "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol",
    "scale.inti", "diff.g"
  )]
  .ctl$trace <- 0
  .ret <- stats::nlminb(
    start = par, objective = fn, gradient = gr, hessian = NULL, control = .ctl,
    lower = lower, upper = upper
  )
  .ret$x <- .ret$par
  ## .ret$message   already there.
  ## .ret$convergence already there.
  .ret
}

.nloptr <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ..., nloptrAlgoritm = "NLOPT_LD_MMA") {
  rxode2::rxReq("nloptr")
  .ctl <- list(
    algorithm = nloptrAlgoritm,
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn, eval_grad_f = gr,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  .ret
}

.bobyqaNLopt <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- list(
    algorithm = "NLOPT_LN_BOBYQA",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  .ret
}

.slsqp <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .nloptr(par, fn, gr, lower, upper, control, ..., nloptrAlgoritm = "NLOPT_LD_SLSQP")
}

.lbfgsbLG <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctlLocal <- list(
    algorithm = "NLOPT_LD_LBFGS",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ctl <- opts <- list(
    "algorithm" = "NLOPT_LD_AUGLAG",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    maxeval = control$maxOuterIterations,
    "local_opts" = .ctlLocal,
    "print_level" = 0
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn, eval_grad_f = gr,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  .ret
}

.rxode2stateOdeNoOutput <- function(x) {
  setdiff(rxode2stateOde(x), "output")
}

.rxInjectMatExpDdt <- function(s) {
  .mv <- rxode2::rxModelVars(s)
  if (!is.list(.mv$indLin) || length(.mv$indLin) != 4L) {
    return(invisible(FALSE))
  }
  .states <- .rxode2stateOdeNoOutput(s)
  if (length(.states) == 0L) {
    return(invisible(FALSE))
  }
  rxode2::.rxInjectMatExpOdes(s)
  .ddt <- stats::setNames(rep("0", length(.states)), .states)
  for (.p in ls(envir = s, all.names = TRUE)) {
    .m <- regexec("^k[_.]([^_.]+)[_.]([^_.]+)$", .p)[[1L]]
    if (length(.m) == 1L) {
      next
    }
    .from <- substring(.p, .m[2L], .m[2L] + attr(.m, "match.length")[2L] - 1L)
    .to <- substring(.p, .m[3L], .m[3L] + attr(.m, "match.length")[3L] - 1L)
    if (.from %in% .states) {
      .ddt[[.from]] <- base::paste0(.ddt[[.from]], "-(", .p, ")*", .from)
    }
    if (.to %in% .states) {
      .ddt[[.to]] <- base::paste0(.ddt[[.to]], "+(", .p, ")*", .from)
    }
  }
  # Append any indLin() forcing functions (e.g. Michaelis-Menten elimination)
  # captured by rxode2::rxS() (stored as per-state rx__indLinForce_<state>__
  # symengine variables) so the emitted d/dt() includes the nonlinear term.
  for (.st in .states) {
    .forceName <- base::paste0("rx__indLinForce_", .st, "__")
    if (base::exists(.forceName, envir = s, inherits = FALSE)) {
      .force <- base::get(.forceName, envir = s, inherits = FALSE)
      .ddt[[.st]] <- base::paste0(.ddt[[.st]], "+(",
                                  rxode2::rxFromSE(.force), ")")
    }
  }
  s$..ddt <- base::paste0("d/dt(", .states, ")=", .ddt)
  invisible(TRUE)
}

#' Get the THETA/ETA lines from rxode2 UI
#'
#' @param rxui This is the rxode2 ui object
#' @return The theta/eta lines
#' @author Matthew L. Fidler
#' @noRd
.uiGetThetaEta <- function(rxui) {
  .iniDf <- rxui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  .etas <- NULL
  if (length(.w) > 0) {
    .thetas <- lapply(.w, function(i) {
      eval(parse(text=paste0("quote(", .iniDf$name[i], " <- THETA[", .iniDf$ntheta[i],"])")))
    })
    .i2 <- .iniDf[-.w, ]
  } else {
    .i2 <- .iniDf
    .thetas <- NULL
  }
  if (length(.i2$name) > 0) {
    .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
    .etas <- lapply(seq_along(.i2$name), function(i) {
      eval(parse(text=paste0("quote(", .i2$name[i], " <- ETA[", .i2$neta1[i], "])")))
    })
  }
  c(.thetas, .etas)
}

#' Get the THETA/ETA params from the rxode2 UI
#'
#' @param rxui This is the rxode2 ui object
#' @return The params eirxode2 UI
#' @author Matthew L. Fidler
#' @noRd
.uiGetThetaEtaParams <- function(rxui, str=FALSE) {
  .iniDf <- rxui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  .etas <- NULL
  if (length(.w) > 0) {
    .thetas <- vapply(.w, function(i) {
      paste0("THETA[", .iniDf$ntheta[i],"]")
    }, character(1), USE.NAMES=FALSE)
    .i2 <- .iniDf[-.w, ]
  } else {
    .etas <- NULL
    .i2 <- .iniDf
    .thetas <- character(0)
  }
  if (length(.i2$name) > 0) {
    .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
    .etas <- vapply(seq_along(.i2$name), function(i) {
      paste0("ETA[", .i2$neta1[i],"]")
    }, character(1), USE.NAMES=FALSE)
  }
  .str <- paste(c(.thetas, .etas, rxui$covariates), collapse=", ")
  if (str) {
    paste0("params(", .str, ")")
  } else {
    eval(parse(text=paste0("quote(params(", .str, "))")))
  }
}

#' @export
rxUiGet.foceiParams <- function(x, ...) {
  .ui <- x[[1]]
  .uiGetThetaEtaParams(.ui, str=TRUE)
}
attr(rxUiGet.foceiParams, "rstudio") <- "params(THETA[1], ETA[1])"

#' @export
rxUiGet.foceiCmtPreModel <- function(x, ...) {
  .ui <- x[[1]]
  .state <- .rxode2stateOdeNoOutput(.ui$mv0)
  if (length(.state) == 0) return("")
  paste(paste0("cmt(", .state, ")"), collapse="\n")
}
attr(rxUiGet.foceiCmtPreModel, "rstudio") <- ""

# This handles the errors for focei
.createFoceiLineObject <- function(x, line) {
  .predDf <- rxUiGet.predDfFocei(list(x, TRUE))
  if (line > nrow(.predDf)) {
    return(NULL)
  }
  .predLine <- .predDf[line, ]
  .ret <- list(x, .predLine, line)
  class(.ret) <- c(paste(.predLine$distribution), "rxGetDistributionFoceiLines")
  .ret
}

#' This is a S3 method for getting the distribution lines for a base rxode2 focei problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the focei. This is based
#'   on the idea that the focei parameters are defined
#' @author Matthew Fidler
#' @keywords internal
#' @export
rxGetDistributionFoceiLines <- function(line) {
  UseMethod("rxGetDistributionFoceiLines")
}

#' Get pred only options
#'
#' @param env  rxode2 environment option
#'
#' @return  If the current method is requesting loglik instead of pred/r
#'  (required for cwres)
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.getRxPredLlikOption <-function() {
  nlmixr2global$rxPredLlik
}

#' @export
rxGetDistributionFoceiLines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  .errNum <- line[[3]]
  if (rxode2hasLlik()) {
    rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1, .errNum,
                                                 rxPredLlik=.getRxPredLlikOption())
  } else {
    rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1)
  }
}

#' @export
rxGetDistributionFoceiLines.t <- function(line) {
  if (rxode2hasLlik()) {
    env <- line[[1]]
    pred1 <- line[[2]]
    .errNum <- line[[3]]
    rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1, .errNum,
                                                 rxPredLlik=.getRxPredLlikOption())
  } else {
    stop("t is not supported", call.=FALSE)
  }
}

#' @export
rxGetDistributionFoceiLines.cauchy <- function(line) {
  if (rxode2hasLlik()) {
    env <- line[[1]]
    pred1 <- line[[2]]
    .errNum <- line[[3]]
    rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1, .errNum,
                                                 rxPredLlik=.getRxPredLlikOption())
  } else {
    stop("t is not supported", call.=FALSE)
  }
}

#' @export
rxGetDistributionFoceiLines.default  <- function(line) {
  if (rxode2hasLlik()) {
    env <- line[[1]]
    pred1 <- line[[2]]
    .errNum <- line[[3]]
    rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1, .errNum,
                                                 rxPredLlik=.getRxPredLlikOption())
  } else {
    stop("unknown distribution", call.=FALSE)
  }
}

#' @export
rxGetDistributionFoceiLines.rxUi <- function(line) {
  .predDf <- rxUiGet.predDfFocei(list(line, TRUE))
  lapply(seq_along(.predDf$cond), function(c) {
    .mod <- .createFoceiLineObject(line, c)
    rxGetDistributionFoceiLines(.mod)
  })
}

#' @export
rxUiGet.foceiModel0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=rxGetDistributionFoceiLines(.f),
                              prefixLines=.uiGetThetaEta(.f),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE)
}
#attr(rxUiGet.foceiModel0, "desc") <- "FOCEi model base"
attr(rxUiGet.foceiModel0, "rstudio") <- quote(rxModelVars({}))

#' @export
rxUiGet.foceiModel0ll <- function(x, ...) {
  nlmixr2global$rxPredLlik <- TRUE
  on.exit(nlmixr2global$rxPredLlik <- FALSE)
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=rxGetDistributionFoceiLines(.f),
                              prefixLines=.uiGetThetaEta(.f),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE)
}

attr(rxUiGet.foceiModel0ll, "rstudio") <- quote(rxModelVars({}))


.foceiPrune <- function(x, fullModel=TRUE) {
  .x <- x[[1]]
  .x <- .x$foceiModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  if (.getRxPredLlikOption()) {
    if (fullModel) {
      .malert(("pruning branches ({.code if}/{.code else}) of llik full model..."))
    } else {
      .malert("pruning branches ({.code if}/{.code else}) of llik model...")
    }
  } else {
    if (fullModel) {
      .malert(("pruning branches ({.code if}/{.code else}) of full model..."))
    } else {
      .malert("pruning branches ({.code if}/{.code else}) of model...")
    }
  }
  .ret <- rxode2::.rxPrune(.x, envir = .env,
                           strAssign=rxode2::rxModelVars(x[[1]])$strAssign)
  .mv <- rxode2::rxModelVars(.ret)
  ## Need to convert to a function
  if (rxode2::.rxIsLinCmt() == 1L) {
    .vars <- c(.mv$params, .mv$lhs, .mv$slhs)
    .mv <- rxode2::.rxLinCmtGen(length(.mv$state), .vars)
  }
  .msuccess("done")
  rxode2::rxNorm(.mv)
}

.loadSymengine <- function(newmod, promoteLinSens = TRUE, fullModel = FALSE) {
  if (.getRxPredLlikOption()) {
    if (fullModel) {
      .malert("loading full llik model into {.pkg symengine} environment...")
    } else {
      .malert("loading llik model into {.pkg symengine} environment...")
    }
  } else {
    if (fullModel) {
      .malert("loading full model into {.pkg symengine} environment...")
    } else {
      .malert("loading into {.pkg symengine} environment...")
    }
  }
  .ret <- rxode2::rxS(newmod, TRUE, promoteLinSens = promoteLinSens)
  if (inherits(.ret$rx_r_, "numeric")) {
    assign("rx_r_", symengine::S(as.character(.ret$rx_r_)), envir=.ret)
  }
  .ret
}

#' @export
rxUiGet.loadPruneSens <- function(x, ...) {
  .loadSymengine(.foceiPrune(x), promoteLinSens = TRUE)
}
#attr(rxUiGet.loadPruneSens, "desc") <- "load sensitivity with linCmt() promoted"
attr(rxUiGet.loadPruneSens, "rstudio") <- emptyenv()

#' @export
rxUiGet.loadPrune <- function(x, ...) {
  .loadSymengine(.foceiPrune(x), promoteLinSens = FALSE)
}
#attr(rxUiGet.loadPrune, "desc") <- "load sensitivity without linCmt() promoted"
attr(rxUiGet.loadPrune, "rstudio") <- emptyenv()

.sensEtaOrTheta <- function(s, theta=FALSE) {
  .etaVars <- NULL
  if (theta && exists("..maxTheta", s)) {
    .etaVars <- paste0("THETA_", seq(1, s$..maxTheta), "_")
  } else if (exists("..maxEta", s)) {
    .etaVars <- paste0("ETA_", seq(1, s$..maxEta), "_")
  }
  if (length(.etaVars) == 0L) {
    stop("cannot identify parameters for sensitivity analysis\n   with nlmixr2 an 'eta' initial estimate must use '~'", call. = FALSE)
  }
  .stateVars <- .rxode2stateOdeNoOutput(s)
  # matExp() models are handled transparently here: rxode2::.rxJacobian calls
  # .rxInjectMatExpOdes(), which materializes the implied d/dt() from the
  # k_from_to rate constants so the standard ODE Jacobian/sensitivity machinery
  # applies.  The original-state d/dt() lines are emitted later by
  # .rxInjectMatExpDdt() in the .rxFinalize* functions.
  rxode2::.rxJacobian(s, c(.stateVars, .etaVars))
  rxode2::.rxSens(s, .etaVars)
  s
}

#' @export
rxUiGet.foceiEtaS <- function(x, ..., theta=FALSE) {
  .s <- rxUiGet.loadPruneSens(x, ...)
  .sensEtaOrTheta(.s)
}
#attr(rxUiGet.foceiEtaS, "desc") <- "Get symengine environment with eta sensitivities"
attr(rxUiGet.foceiEtaS, "rstudio") <- emptyenv()


#' @export
rxUiGet.foceiThetaS <- function(x, ..., theta=FALSE) {
  .s <- rxUiGet.loadPruneSens(x, ...)
  .sensEtaOrTheta(.s, theta=TRUE)
}
#attr(rxUiGet.foceiEtaS, "desc") <- "Get symengine environment with eta sensitivities"
attr(rxUiGet.foceiThetaS, "rstudio") <- emptyenv()

#' @export
rxUiGet.foceiHdEta <- function(x, ...) {
  .s <- rxUiGet.foceiEtaS(x)
  .stateVars <- .rxode2stateOdeNoOutput(.s)
  # FIXME: take out pred.minus.dv
  .predMinusDv <- rxode2::rxGetControl(x[[1]], "predMinusDv", TRUE)
  .grd <- rxode2::rxExpandFEta_(
    .stateVars, .s$..maxEta,
    ifelse(.predMinusDv, 1L, 2L)
  )
  if (rxode2::.useUtf()) {
    .malert("calculate \u2202(f)/\u2202(\u03B7)")
  } else {
    .malert("calculate d(f)/d(eta)")
  }
  rxode2::rxProgress(dim(.grd)[1])
  on.exit({
    rxode2::rxProgressAbort()
  })
  .any.zero <- FALSE
  .all.zero <- TRUE
  .ret <- apply(.grd, 1, function(x) {
    .l <- x["calc"]
    .l <- eval(parse(text = .l))
    .ret <- paste0(x["dfe"], "=", rxode2::rxFromSE(.l))
    .zErr <- suppressWarnings(try(as.numeric(get(x["dfe"], .s)), silent = TRUE))
    if (identical(.zErr, 0)) {
      .any.zero <<- TRUE
    } else if (.all.zero) {
      .all.zero <<- FALSE
    }
    rxode2::rxTick()
    .ret
  })
  if (.all.zero) {
    stop("none of the predictions depend on 'ETA'", call. = FALSE)
  }
  if (.any.zero) {
    warning("some of the predictions do not depend on 'ETA'", call. = FALSE)
  }
  .s$..HdEta <- .ret
  .s$..pred.minus.dv <- .predMinusDv
  rxode2::rxProgressStop()
  .s
}
attr(rxUiGet.foceiHdEta, "desc") <- "Generate the d(err)/d(eta) values for FO related methods"
attr(rxUiGet.foceiHdEta, "rstudio") <- emptyenv()


#' Finalize inner rxode2 based on symengine saved info
#'
#' @param .s Symengine/rxode2 object
#' @return Nothing
#' @author Matthew L Fidler
#' @noRd
.rxFinalizeInner <- function(.s, sum.prod = FALSE,
                             optExpression = TRUE) {
  .isMatExp <- isTRUE(.rxInjectMatExpDdt(.s))
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .r <- get("rx_r_", envir = .s)
  .r <- paste0("rx_r_=", rxode2::rxFromSE(.r))
  .yj <- paste(get("rx_yj_", envir = .s))
  .yj <- paste0("rx_yj_~", rxode2::rxFromSE(.yj))
  .lambda <- paste(get("rx_lambda_", envir = .s))
  .lambda <- paste0("rx_lambda_~", rxode2::rxFromSE(.lambda))
  .hi <- paste(get("rx_hi_", envir = .s))
  .hi <- paste0("rx_hi_~", rxode2::rxFromSE(.hi))
  .low <- paste(get("rx_low_", envir = .s))
  .low <- paste0("rx_low_~", rxode2::rxFromSE(.low))
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- character(0)
  .lhs <- .s$..lhs
  if (is.null(.lhs)) .lhs <- character(0)
  .sens <- .s$..sens
  if (is.null(.sens)) .sens <- character(0)
  # Only matExp() models need the model LHS here: it defines the k_from_to rate
  # constants that the materialized d/dt() lines reference.  For ordinary models
  # the d/dt()/sensitivity equations are self-contained, so the LHS is omitted.
  # The LHS is emitted as suppressed assignments ('~' not '=') so it does not add
  # output columns -- extra output columns shift the column layout the FOCEi C++
  # reads and corrupt the inner objective.
  .preLhs <- if (.isMatExp) sub("^([^=]+)=", "\\1~", .lhs) else character(0)
  .s$..inner <- paste(c(
    .preLhs,
    .ddt,
    .sens,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .s$..HdEta,
    .r,
    .s$..REta,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")
  .s$..innerOeta <- paste(c(
    .preLhs,
    .ddt,
    .sens,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .s$..HdEta,
    .r,
    .s$..REta,
    paste0("rx__ETA", seq_len(.s$..maxEta), "=ETA[",seq_len(.s$..maxEta), "]"),
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""))
  if (sum.prod) {
    .malert("stabilizing round off errors in inner problem...")
    .s$..inner <- rxode2::rxSumProdModel(.s$..inner)
    .s$..innerOeta <- rxode2::rxSumProdModel(.s$..innerOeta)
    .msuccess("done")
  }
  if (optExpression) {
    .s$..inner <- rxode2::rxOptExpr(.s$..inner,
                                    ifelse(.getRxPredLlikOption(),
                                           "inner llik model",
                                           "inner model"))
    suppressMessages(.s$..innerOeta <- rxode2::rxOptExpr(.s$..innerOeta,
                                                         ifelse(.getRxPredLlikOption(),
                                                                "inner llik model",
                                                                "inner model")))
  }
}

#' @export
rxUiGet.foceiEnv <- function(x, ...) {
  .s <- rxUiGet.foceiHdEta(x, ...)
  .stateVars <- .rxode2stateOdeNoOutput(.s)
  .grd <- rxode2::rxExpandFEta_(.stateVars, .s$..maxEta, FALSE)
  if (rxode2::.useUtf()) {
    .malert("calculate \u2202(R\u00B2)/\u2202(\u03B7)")
  } else {
    .malert("calculate d(R^2)/d(eta)")
  }
  rxode2::rxProgress(dim(.grd)[1])
  on.exit({
    rxode2::rxProgressAbort()
  })
  .ret <- apply(.grd, 1, function(x) {
    .l <- x["calc"]
    .l <- eval(parse(text = .l))
    .ret <- paste0(x["dfe"], "=", rxode2::rxFromSE(.l))
    rxode2::rxTick()
    .ret
  })

  .s$..REta <- .ret
  rxode2::rxProgressStop()
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizeInner(.s, .sumProd, .optExpression)
  .rxFinalizePred(.s, .sumProd, .optExpression)
  .s$..outer <- NULL
  .s
}
#attr(rxUiGet.foceiEnv, "desc") <- "Get the focei environment"
attr(rxUiGet.foceiEnv, "rstudio") <- emptyenv()

#' @export
rxUiGet.foceEnv <- function(x, ...) {
  .s <- rxUiGet.foceiHdEta(x, ...)
  .s$..REta <- NULL
  ## Take etas from rx_r
  eval(parse(text = rxode2::rxRepR0_(.s$..maxEta)))
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizeInner(.s, .sumProd, .optExpression)
  .rxFinalizePred(.s, .sumProd, .optExpression)
  .s$..outer <- NULL
  .s
}
#attr(rxUiGet.foceEnv, "desc") <- "Get the foce environment"
attr(rxUiGet.foceEnv, "rstudio") <- emptyenv()


#' @export
rxUiGet.getEBEEnv <- function(x, ...) {
  .s <- rxUiGet.loadPrune(x, ...)
  .s$..inner <- NULL
  .s$..innerOeta <- NULL
  .s$..outer <- NULL
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizePred(.s, .sumProd, .optExpression)
  .s
}
#attr(rxUiGet.getEBEEnv, "desc") <- "Get the EBE environment"
attr(rxUiGet.getEBEEnv, "rstudio") <- emptyenv()

.toRx <- function(x, msg, eventSens = "fd") {
  if (is.null(x)) {
    return(NULL)
  }
  .malert(msg)
  ## eventSens="jump" attaches rxode2's analytic event ("jump") sensitivity
  ## information to the model so the dosing-parameter (alag/F/rate/dur/...)
  ## sensitivities are computed analytically rather than by finite differences.
  ## Passed only for models that carry the sensitivity equations (the inner
  ## model); "fd" everywhere else preserves the legacy behavior.
  .ret <- rxode2::rxode2(paste(nlmixr2global$toRxParam, x,
                               nlmixr2global$toRxDvidCmt),
                         eventSens = eventSens)
  .msuccess("done")
  .ret
}

.nullInt <- function(x) {
  if (rxode2::rxIs(x, "integer") || rxode2::rxIs(x, "numeric")) {
    as.integer(x)
  } else {
    integer(0)
  }
}

#' @export
rxUiGet.predDfFocei <- function(x, ...) {
  .ui <- x[[1]]
  if (exists(".predDfFocei", envir=.ui)) {
    get(".predDfFocei", envir=.ui)
  } else {
    .predDf <- .ui$predDf
    if (all(.predDf$distribution == "norm")) {
      assign(".predDfFocei", .predDf, envir=.ui)
      .predDf
    } else {
      .w <- which(.predDf$distribution == "norm")
      if (length(.w) > 0) {
        .predDf$distribution[.w] <- "dnorm"
      }
      assign(".predDfFocei", .predDf, envir=.ui)
      .predDf
    }
  }
}
attr(rxUiGet.predDfFocei, "rstudio") <- NA



.rxFinalizePred <- function(.s, sum.prod = FALSE,
                            optExpression = TRUE) {
  .isMatExp <- isTRUE(.rxInjectMatExpDdt(.s))
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .r <- get("rx_r_", envir = .s)
  .r <- paste0("rx_r_=", rxode2::rxFromSE(.r))
  .yj <- paste(get("rx_yj_", envir = .s))
  .yj <- paste0("rx_yj_~", rxode2::rxFromSE(.yj))
  .lambda <- paste(get("rx_lambda_", envir = .s))
  .lambda <- paste0("rx_lambda_~", rxode2::rxFromSE(.lambda))
  .hi <- paste(get("rx_hi_", envir = .s))
  .hi <- paste0("rx_hi_~", rxode2::rxFromSE(.hi))
  .low <- paste(get("rx_low_", envir = .s))
  .low <- paste0("rx_low_~", rxode2::rxFromSE(.low))
  .lhs0 <- .s$..lhs0
  if (is.null(.lhs0)) .lhs0 <- ""
  .lhs <- .s$..lhs
  if (is.null(.lhs)) .lhs <- ""
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- ""
  # For matExp() models the model LHS defines the k_from_to rate constants that
  # the materialized d/dt() lines reference, so the LHS must precede the d/dt().
  # It is emitted suppressed ('~' not '=') so it does not add output columns.
  # Other models keep the LHS after the prediction (some error-model LHS depend
  # on rx_pred_).
  .preLhs <- if (.isMatExp) sub("^([^=]+)=", "\\1~", .lhs) else character(0)
  .postLhs <- if (.isMatExp) character(0) else .lhs
  .s$..pred <- paste(c(
    .s$..stateInfo["state"],
    .lhs0,
    .preLhs,
    .ddt,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .r,
    .postLhs,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    "tad=tad()",
    "dosenum=dosenum()",
    ""
  ), collapse = "\n")
  .s$..pred.nolhs <- paste(c(
    .s$..stateInfo["state"],
    .lhs0,
    .preLhs,
    .ddt,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .r,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")
  if (sum.prod) {
    .malert("stabilizing round off errors in predictions or EBE model...")
    .s$..pred <- rxode2::rxSumProdModel(.s$..pred)
    .msuccess("done")
  }
  if (optExpression) {
    .s$..pred <- rxode2::rxOptExpr(.s$..pred,
                                   ifelse(.getRxPredLlikOption(),"Llik EBE model","EBE model"))
  }
}

.innerInternal <- function(ui, s) {
  .cmt <-  ui$foceiCmtPreModel
  .interp <- ui$interpLinesStr
  if (.interp != "") {
    .cmt <-paste0(.cmt, "\n", .interp)
  }
  nlmixr2global$toRxParam <-
    paste0(.uiGetThetaEtaParams(ui, TRUE), "\n",
           .cmt, "\n")
  nlmixr2global$toRxDvidCmt <- .foceiToCmtLinesAndDvid(ui)
  if (exists("..maxTheta", s)) {
    .eventTheta <- rep(0L, s$..maxTheta)
  } else {
    .eventTheta <- integer()
  }
  if (exists("..maxEta", s)) {
    .eventEta <- rep(0L, s$..maxEta)
  } else {
    .eventEta <- integer()
  }
  ## Event-sensitivity method.  "jump" enables rxode2's analytic dosing-parameter
  ## (alag/F/rate/dur) sensitivities.
  .eventSens <- rxode2::rxGetControl(ui, "eventSens", "jump")
  ## `eventEta`/`eventTheta` flag the parameters that enter a dosing expression
  ## (alag/F/rate/dur).  In the legacy "fd" path inner.cpp computes their
  ## sensitivity by finite differences (predOde) because the analytic `rx__sens`
  ## states miss the event jump.  Under "jump" rxode2 injects the analytic jump
  ## into those `rx__sens` states, so the analytic gradient is now correct and
  ## the finite-difference fallback must be turned OFF -- otherwise the
  ## jump-corrected sensitivity is computed but never used.  Leaving the flags at
  ## zero routes every parameter through the analytic innerOde sensitivity.
  if (!identical(.eventSens, "jump")) {
    for (.v in s$..eventVars) {
      .vars <- as.character(get(.v, envir = s))
      .vars <- rxode2::rxGetModel(paste0("rx_lhs=", rxode2::rxFromSE(.vars)))$params
      for (.v2 in .vars) {
        .reg <- rex::rex(start, "ETA[", capture(any_numbers), "]", end)
        if (regexpr(.reg, .v2) != -1) {
          .num <- as.numeric(sub(.reg, "\\1", .v2))
          .eventEta[.num] <- 1L
        }
        .reg <- rex::rex(start, "THETA[", capture(any_numbers), "]", end)
        if (regexpr(.reg, .v2) != -1) {
          .num <- as.numeric(sub(.reg, "\\1", .v2))
          .eventTheta[.num] <- 1L
        }
      }
    }
  }
  pred.opt <- NULL
  ## Build the inner (sensitivity) model with the requested event-sensitivity
  ## method.  "jump" enables rxode2's analytic dosing-parameter sensitivities.
  inner <- .toRx(s$..inner, "compiling inner model...", eventSens = .eventSens)
  innerOeta <- s$..innerOeta
  .sumProd <- rxode2::rxGetControl(ui, "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(ui, "optExpression", TRUE)
  .predMinusDv <- rxode2::rxGetControl(ui, "predMinusDv", TRUE)
  if (!is.null(inner)) {
    if (.sumProd) {
      .malert("stabilizing round off errors in FD model...")
      s$..pred.nolhs <- rxode2::rxSumProdModel(s$..pred.nolhs)
      .msuccess("done")
    }
    if (.optExpression) {
      s$..pred.nolhs <- rxode2::rxOptExpr(s$..pred.nolhs,
                                          ifelse(.getRxPredLlikOption(),"Llik FD model","FD model"))
    }
    s$..pred.nolhs <- paste(c(
      paste0("params(", paste(inner$params, collapse = ","), ")"),
      s$..pred.nolhs
    ), collapse = "\n")
    pred.opt <- s$..pred.nolhs
  }
  # For mixture models build predOnly from the pruned model (which preserves the
  # mix() call and therefore gives nMix > 0 in the compiled model).  This lets
  # rxode2 accept per-individual mixest from iCov so that me/mn/mu are correct
  # and IPRED uses the right mixture branch for each subject.
  # NOTE: the *inner* model intentionally keeps the mixest==k symengine form
  # (nMix == 0) because inner.cpp manages mixture selection itself; using mix()
  # there would trigger a double-optimisation conflict.
  .mixProbs <- try(ui$mixProbs, silent=TRUE)
  .hasMix <- !inherits(.mixProbs, "try-error") && length(.mixProbs) > 0L
  .predOnly <- if (.hasMix) {
    .prunedStr <- paste(c(.foceiPrune(list(ui)), "tad=tad()", "dosenum=dosenum()", ""),
                        collapse="\n")
    .toRx(.prunedStr, ifelse(.getRxPredLlikOption(),
                             "compiling Llik EBE model (mixture)...",
                             "compiling EBE model (mixture)..."))
  } else {
    .toRx(s$..pred, ifelse(.getRxPredLlikOption(),
                           "compiling Llik EBE model...",
                           "compiling EBE model..."))
  }
  .ret <- list(
    inner = inner,
    innerOeta = innerOeta,
    predOnly = .predOnly,
    extra.pars = s$..extraPars,
    outer = .toRx(s$..outer),
    predNoLhs = .toRx(pred.opt, ifelse(.getRxPredLlikOption(),
                                       "compiling events Llik FD model...",
                                       "compiling events FD model...")),
    theta = NULL,
    ## warn=.zeroSens,
    pred.minus.dv = .predMinusDv,
    log.thetas = .nullInt(s$..extraTheta[["exp"]]),
    log.etas = .nullInt(s$..extraEta[["exp"]]),
    extraProps = s$..extraTheta,
    eventTheta = .eventTheta,
    eventEta = .eventEta
    ## ,
    ## cache.file=cache.file
  )
  class(.ret) <- "foceiModelList"
  .ret
}

#' @export
rxUiGet.focei <- function(x, ...) {
  .ui <- x[[1]]
  # For t/cauchy/dnorm, predOnly model
  nlmixr2global$rxPredLlik <- FALSE
  on.exit(nlmixr2global$rxPredLlik <- FALSE)
  .s <- rxUiGet.foceiEnv(x, ...)
  .ret <-  .innerInternal(.ui, .s)
  .predDf <- .ui$predDfFocei
  if (any(.predDf$distribution %in% c("t", "cauchy", "dnorm"))) {
    nlmixr2global$rxPredLlik <- TRUE
    .s <- rxUiGet.foceiEnv(x, ...)
    .s2 <- .innerInternal(.ui, .s)
    .w <- vapply(seq_along(.s2),
                 function(i) {
                   inherits(.s2[[i]], "rxode2")
                 }, logical(1), USE.NAMES=FALSE)
    .s2 <- .s2[.w]
    names(.s2) <- paste0(names(.s2), "Llik")
    .cls <- class(.ret)
    .ret <- c(.ret, .s2)
    class(.ret) <-.cls
  }
  .ret
}
#attr(rxUiGet.focei, "desc") <- "Get the FOCEi foceiModelList object"

#' @export
rxUiGet.foce <- function(x, ...) {
  .ui <- x[[1]]
  nlmixr2global$rxPredLlik <- FALSE
  on.exit(nlmixr2global$rxPredLlik <- FALSE)
  .s <- rxUiGet.foceEnv(x, ...)
  .ret <- .innerInternal(.ui, .s)
  .predDf <- .ui$predDfFocei
  if (any(.predDf$distribution %in% c("t", "cauchy", "dnorm"))) {
    nlmixr2global$rxPredLlik <- TRUE
    .s <- rxUiGet.foceEnv(x, ...)
    .s2 <- .innerInternal(.ui, .s)
    .w <- vapply(seq_along(.s2),
                 function(i) {
                   inherits(.s2[[i]], "rxode2")
                 }, logical(1), USE.NAMES=FALSE)
    .s2 <- .s2[.w]
    names(.s2) <- paste0(names(.s2), "Llik")
    .cls <- class(.ret)
    .ret <- c(.ret, .s2)
    class(.ret) <-.cls
  }
  .ret
}
#attr(rxUiGet.foce, "desc") <- "Get the FOCE foceiModelList object"


#' @export
rxUiGet.ebe <- function(x, ...) {
  .ui <-x[[1]]
  nlmixr2global$rxPredLlik <- FALSE
  on.exit(  nlmixr2global$rxPredLlik <- FALSE)
  .s <- rxUiGet.getEBEEnv(x, ...)
  .ret <- .innerInternal(.ui, .s)
  .predDf <- .ui$predDfFocei
  if (any(.predDf$distribution %in% c("t", "cauchy", "dnorm"))) {
    nlmixr2global$rxPredLlik <- TRUE
    .s <- rxUiGet.getEBEEnv(x, ...)
    .s2 <- .innerInternal(.ui, .s)
    .w <- vapply(seq_along(.s2),
                 function(i) {
                   inherits(.s2[[i]], "rxode2")
                 }, logical(1), USE.NAMES=FALSE)
    .s2 <- .s2[.w]
    names(.s2) <- paste0(names(.s2), "Llik")
    .cls <- class(.ret)
    .ret <- c(.ret, .s2)
    class(.ret) <-.cls
  }
  .ret
}
#attr(rxUiGet.ebe, "desc") <- "Get the EBE foceiModelList object"

#' @export
rxUiGet.foceiModelDigest <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- get("iniDf", .ui)
  .sumProd <- rxode2::rxGetControl(.ui, "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(.ui, "optExpression", TRUE)
  .predMinusDv   <- rxode2::rxGetControl(.ui, "predMinusDv", TRUE)
  ## eventSens changes the inner model codegen (analytic jump sensitivities) and
  ## the eventEta/eventTheta finite-difference flags, so it must be part of the
  ## cache key -- otherwise a "jump" build would reuse a cached "fd" model.
  .eventSens <- rxode2::rxGetControl(.ui, "eventSens", "jump")
  digest::digest(c(all(is.na(.iniDf$neta1)),
                   rxode2::rxGetControl(.ui, "interaction", 1L),
                   .iniDf$name,
                   .sumProd, .optExpression, .predMinusDv,
                   .eventSens,
                   rxode2::rxGetControl(.ui, "addProp", getOption("rxode2.addProp", "combined2")),
                   .ui$lstExpr))
}
#attr(rxUiGet.foceiModelDigest, "desc") <- "Get the md5 digest for the focei model"
attr(rxUiGet.foceiModelDigest, "rstudio") <- "hash"

#' @export
rxUiGet.foceiModelCache <- function(x, ...) {
  file.path(rxode2::rxTempDir(),
            paste0("focei-", rxUiGet.foceiModelDigest(x, ...), ".qs2"))
}
#attr(rxUiGet.foceiModelCache, "desc") <- "Get the focei cache file for a model"
attr(rxUiGet.foceiModelCache, "rstudio") <- "file"

#' @export
rxUiGet.foceiModel <- function(x, ...) {
  .cacheFile <- rxUiGet.foceiModelCache(x, ...)
  if (file.exists(.cacheFile)) {
    .ret <- qs2::qs_read(.cacheFile)
    lapply(seq_along(.ret), function(i) {
      if (inherits(.ret[[i]], "rxode2")) {
        rxode2::rxLoad(.ret[[i]])
      }
    })
    return(.ret)
  }
  .ui <- x[[1]]
  .iniDf <- get("iniDf", .ui)
  if (all(is.na(.iniDf$neta1))) {
    .ret <- rxUiGet.ebe(x, ...)
  } else {
    if (rxode2::rxGetControl(.ui, "interaction", 1L)) {
      .ret <- rxUiGet.focei(x, ...)
    } else {
      .ret <- rxUiGet.foce(x, ...)
    }
  }
  qs2::qs_save(.ret, .cacheFile)
  .ret
}
# attr(rxUiGet.foceiModel, "desc") <- "Get focei model object"

#' @export
rxUiGet.foceiFixed <- function(x, ...) {
  .x <- x[[1]]
  .df <- get("iniDf", .x)
  .dft <- .df[!is.na(.df$ntheta), ]
  .fix <- .dft$fix
  .dft <- .df[is.na(.df$ntheta), ]
  c(.fix, .dft$fix)
}
#attr(rxUiGet.foFixed, "desc") <- "focei theta fixed vector"
attr(rxUiGet.foceiFixed, "rstudio") <- c(FALSE, TRUE)

#' @export
rxUiGet.foceiEtaNames <- function(x, ...) {
  .x <- x[[1]]
  .df <- get("iniDf", .x)
  .dft <- .df[is.na(.df$ntheta), ]
  .dft[.dft$neta1 == .dft$neta2, "name"]
}
#attr(rxUiGet.foceiEtaNames, "desc") <- "focei eta names"
attr(rxUiGet.foceiEtaNames, "rstudio") <- c("eta.ka", "eta.cl", "eta.vc")

#' This assigns the tolerances based on a different tolerance for the
#' sensitivity equations
#'
#' It will update and modify the control inside of the UI.
#'
#' It also updates the predNeq that is needed for numeric derivatives
#'
#' @param ui rxode2 UI object
#' @param env focei environment for solving
#' @return Called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiOptEnvAssignTol <- function(ui, env) {
  .len <- length(env$model$predNoLhs$state)
  rxode2::rxAssignControlValue(ui, "predNeq", .len)
  if (!is.null(env$model$inner)) {
    .len0 <- length(env$model$inner$state)
    .len2 <- .len0 - .len
    if (.len2 > 0) {
      .env <- nlmixr2global$nlmixrEvalEnv$envir
      if (!is.environment(.env)) {
        .env <- parent.frame(1)
      }
      .rxControl <- rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl())
      rxode2::rxAssignControlValue(ui, "rxControl", rxode2::rxControlUpdateSens(.rxControl, .len2, .len0))
    }
  }
}
#' Assign the number of log likelihood items that need to be allocated
#'
#' @param ui rxode2 ui
#' @param env optimization environment
#' @return Nothing called for side effects.  Will update env$rxControl
#'   to have the maximum number of llik items in the model set.
#' @author Matthew L. Fidler
#' @noRd
.foceiOptEnvAssignNllik <- function(ui, env) {
  if (rxode2hasLlik()) {
    .maxLl <- max(vapply(seq_along(env$model), function(i) {
      .model <- env$model[[i]]
      if (inherits(.model, "rxode2")) {
        rxode2::rxModelVars(.model)$flags["nLlik"]
      } else {
        0L
      }
    }, integer(1), USE.NAMES=FALSE))
    if (.maxLl > 0) {
      .env <- nlmixr2global$nlmixrEvalEnv$envir
      if (!is.environment(.env)) {
        .env <- parent.frame(1)
      }
      .rxControl <- rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl())
      .rxControl$nLlikAlloc <- .maxLl
      rxode2::rxAssignControlValue(ui, "rxControl", .rxControl)
    }
  }
}

#'  This sets up the initial omega/eta estimates and the boundaries for the whole system
#'
#' @param ui rxode2 UI object
#' @param env focei solving environment
#' @return NoHing, called for side effecs
#' @author Matthew L. Fidler
#' @noRd
.foceiOptEnvSetupBounds <- function(ui, env) {
  .iniDf <- ui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  if (length(.w) > 0) {
    .lower <- vapply(.w,
                     function(i) {
                       .low <- .iniDf$lower[i]
                       .zeroRep <- rxode2::rxGetControl(ui, "sdLowerFact", 0.001)
                       if (.zeroRep <= 0) return(.low)
                       if (.low <= 0 &&
                             .iniDf$err[i] %in% c("add",
                                                  "lnorm", "logitNorm", "probitNorm",
                                                  "prop", "propT", "propF",
                                                  "pow", "powF", "powT")) {
                         .low <- .iniDf$est[i] * 0.001
                       }
                       .low
                     }, numeric(1), USE.NAMES=FALSE)
    .upper <- .iniDf$upper[.w]
    env$thetaIni <- ui$thetaIniMix
    env$mixIdx <- ui$thetaMixIndex
    env$thetaIni <- setNames(env$thetaIni, paste0("THETA[", seq_along(env$thetaIni), "]"))
  } else {
    .lower <- numeric(0)
    .upper <- numeric(0)
    env$mixIdx <- integer(0)
    env$thetaIni <- setNames(numeric(0), character(0))
  }
  rxode2::rxAssignControlValue(ui, "nfixed", sum(ui$iniDf$fix))
  .mixed <- !is.null(env$etaNames)
  if (.mixed && length(env$etaNames) == 0L) .mixed <- FALSE
  if (!.mixed) {
    rxode2::rxAssignControlValue(ui, "nomega", 0)
    rxode2::rxAssignControlValue(ui, "neta", 0)
    env$xType <- -1
    rxode2::rxAssignControlValue(ui, "ntheta", length(ui$iniDf$lower))
  } else {
    .om0 <- ui$omega
    .diagXform <- rxode2::rxGetControl(ui, "diagXform", "sqrt")
    env$rxInv <- rxode2::rxSymInvCholCreate(mat = .om0, diag.xform = .diagXform)
    env$xType <- env$rxInv$xType
    .om0a <- .om0
    .om0a <- .om0a / rxode2::rxGetControl(ui, "diagOmegaBoundLower", 100)
    .om0b <- .om0
    .om0b <- .om0b * rxode2::rxGetControl(ui, "diagOmegaBoundUpper", 5)
    .om0a <- rxode2::rxSymInvCholCreate(mat = .om0a, diag.xform = .diagXform)
    .om0b <- rxode2::rxSymInvCholCreate(mat = .om0b, diag.xform = .diagXform)
    .omdf <- data.frame(a = .om0a$theta, m = env$rxInv$theta, b = .om0b$theta, diag = .om0a$theta.diag)
    .omdf$lower <- with(.omdf, ifelse(a > b, b, a))
    .omdf$lower <- with(.omdf, ifelse(lower == m, -Inf, lower))
    .omdf$lower <- with(.omdf, ifelse(!diag, -Inf, lower))
    .omdf$upper <- with(.omdf, ifelse(a < b, b, a))
    .omdf$upper <- with(.omdf, ifelse(upper == m, Inf, upper))
    .omdf$upper <- with(.omdf, ifelse(!diag, Inf, upper))
    rxode2::rxAssignControlValue(ui, "nomega", length(.omdf$lower))
    rxode2::rxAssignControlValue(ui, "neta", sum(.omdf$diag))
    rxode2::rxAssignControlValue(ui, "ntheta", length(.lower))
    .lower <- c(.lower, .omdf$lower)
    .upper <- c(.upper, .omdf$upper)
  }
  env$lower <- .lower
  env$upper <- .upper
  .etaMat <- rxode2::rxGetControl(ui, "etaMat", NULL)
  if (length(.etaMat) == 1L && is.na(.etaMat)) .etaMat <- NULL
  env$etaMat <- .etaMat
  env
}

#' Setup the scaleC
#'
#' @param ui rxode2 UI
#' @param env Focei setup environment
#' @return NoHing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiOptEnvSetupScaleC <- function(ui, env) {
  .controlScaleC <- rxode2::rxGetControl(ui, "scaleC", NULL)
  .len <- length(env$lower)
  if (is.null(.controlScaleC)) {
    .scaleC <- rep(NA_real_, .len)
  } else {
    .scaleC <- as.double(.controlScaleC)
  }
  .lenC <- length(.scaleC)
  if (.len > .lenC) {
    .scaleC <- c(.scaleC, rep(NA_real_, .len - .lenC))
  } else if (.len < .lenC) {
    .scaleC <- .scaleC[seq(1, .lenC)]
    warning("'scaleC' control option has more options than estimated population parameters, please check",
            call.=FALSE)
  }

  .ini <- ui$iniDf
  .ini <- .ini[!is.na(.ini$err), c("est", "err", "ntheta")]
  for (.i in seq_along(.ini$err)) {
    if (is.na(.scaleC[.ini$ntheta[.i]])) {
      if (any(.ini$err[.i] == c("boxCox", "yeoJohnson", "pow2", "tbs", "tbsYj"))) {
        .scaleC[.ini$ntheta[.i]] <- 1
      } else if (any(.ini$err[.i] == c("prop", "add", "norm", "dnorm", "logn", "dlogn", "lnorm", "dlnorm"))) {
        .scaleC[.ini$ntheta[.i]] <- 0.5 * abs(.ini$est[.i])
      }
    }
  }
  .muRefCurEval <- ui$muRefCurEval
  .ini <- ui$iniDf
  for (.i in seq_along(.muRefCurEval$parameter)) {
    .curEval <- .muRefCurEval$curEval[.i]
    .par <- .muRefCurEval$parameter[.i]
    .w <- which(.ini$name == .par)
    if (length(.w) == 1) {
      if (!is.na(.ini$ntheta[.w])) {
        .j <- .ini$ntheta[.w]
        if (is.na(.scaleC[.j])) {
          # These have similar deriavtes on a log scale.
          if (.curEval == "exp") {
            # Hence D(S("log(exp(x))"}, "x")
            .scaleC[.j] <- 1 # log scaled
          } else if (.curEval == "factorial") {
            # Hence 1/D(S("log(factorial(x))"}, "x"):
            .scaleC[.j] <- abs(1 / digamma(.ini$est[.j] + 1))
          } else if (.curEval == "gamma") {
            #1/D(log(gamma(x)), x)
            .scaleC[.j] <- abs(1 / digamma(.ini$est[.j]))
          } else if (.curEval == "log") {
            #1/D(log(log(x)), x)
            .scaleC[.j] <- log(abs(.ini$est[.j])) * abs(.ini$est[.j])
          } else if (.curEval == "logit") {
            # 1/D(log(logit(x, a, b)))
            .a <- .muRefCurEval$low[.i]
            .b <- .muRefCurEval$hi[.i]
            .x <- .ini$est[.j]
            .scaleC[.j] <- -1.0*(-.a + .x)^2*(-1.0 + 1.0*(-.a + .b)/(-.a + .x))*log(abs(-1.0 + 1.0*(-.a + .b)/(-.a + .x)))/(-.a + .b)
          } else if (.curEval == "expit") {
            # 1/D(log(expit(x, a, b)))
            .a <- .muRefCurEval$low[.i]
            .b <- .muRefCurEval$hi[.i]
            .x <- .ini$est[.j]
            .scaleC[.j] <- 1.0*exp(.x)*(1.0 + exp(-.x))^2*(.a + 1.0*(-.a + .b)/(1.0 + exp(-.x)))/(-.a + .b)
          } else if (.curEval == "probitInv") {
            .a <- .muRefCurEval$low[.i]
            .b <- .muRefCurEval$hi[.i]
            .x <- .ini$est[.j]
            .scaleC[.j] <- 1.4142135623731*exp(0.5*.x^2)*sqrt(pi)*(.a + 0.5*(-.a + .b)*(1.0 + rxode2::erf(0.707106781186547*.x)))/(-.a + .b)
          } else if (.curEval == "probit") {
            .a <- .muRefCurEval$low[.i]
            .b <- .muRefCurEval$hi[.i]
            .x <- .ini$est[.j]
            erfinvF  <- function(y) {
              if(abs(y) > 1) return(NA_real_)
              sqrt(qchisq(abs(y),1)/2) * sign(y)
            }
            .scaleC[.j] <- sqrt(2)*(-.a+.b)*erfinvF(-1+2*(-.a+.x)/(-.a+.b))/sqrt(pi)/2*exp(((erfinvF(-1+2*(-.a+.x)/(-.a+.b))) ^ 2))
          } else if (is.na(.curEval) || .curEval == "") {
            # Additive mu-reference (theta + eta).  Scale by the
            # magnitude of the initial estimate so a unit step in the
            # internal (scaled) parameter is proportional to the
            # parameter's natural scale.  Falling through to the C++
            # default of 1/|init| produced steps so small the optimizer
            # could not move parameters with large-magnitude initials
            # (issue 641).
            .val <- abs(.ini$est[.j])
            if (.val > 1) .scaleC[.j] <- .val
          }
        }
      }
    }
  }
  env$scaleC <- .scaleC
}

#' @export
rxUiGet.scaleCtheta <- function(x, ...) {
  .ui <- x[[1]]
  .env <- new.env(parent=emptyenv())
  .env$lower <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), "lower"]
  .foceiOptEnvSetupScaleC(.ui, .env)
  .env$scaleC[!.ui$iniDf$fix]
}
attr(rxUiGet.scaleCtheta, "rstudio") <- c(1.0, NA_real_)

#' @export
rxUiGet.scaleCnls <- function(x, ...) {
  .ui <- x[[1]]
  .env <- new.env(parent=emptyenv())
  .env$lower <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), "lower"]
  .foceiOptEnvSetupScaleC(.ui, .env)
  .env$scaleC[!.ui$iniDf$fix & !(.ui$iniDf$err %in% c("add", "prop", "pow"))]
}
attr(rxUiGet.scaleCnls, "rstudio") <- c(1.0, NA_real_)

# focei.mu.ref
# eta# and the corresponding theta number

#' @export
rxUiGet.foceiMuRefVector <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .muRefDataFrame <- .ui$muRefDataFrame
  .w <- which(!is.na(.iniDf$ntheta))
  .i2 <- .iniDf[-.w, ]
  if (length(.i2$name) > 0) {
    .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
    .i2 <- .i2[order(.i2$neta1), ]
    vapply(seq_along(.i2$neta1), function(i) {
      if (.i2$fix[i]) return(-1L)
      .name <- .i2$name[i]
      .w <- which(.muRefDataFrame$eta == .name)
      if (length(.w) != 1) return(-1L)
      .name <- .muRefDataFrame$theta[.w]
      .w <- which(.iniDf$name == .name)
      if (length(.w) != 1) return(-1L)
      if (.iniDf$fix[.w]) return(-1L)
      .iniDf$ntheta[.w] - 1L
    }, integer(1))
  } else {
    integer(0)
  }
}
#attr(rxUiGet.foceiMuRefVector, "desc") <- "focei mu ref vector"
attr(rxUiGet.foceiMuRefVector, "rstudio") <- c(0L, -1L)

# focei.mu.cov.eta
# For the mu-referenced FOCEI family (mufocei/irlsfocei/...): a 0/1 flag per
# eta, same length/ordering as foceiMuRefVector, marking which etas are
# mu-ref-covariate-eligible (see .muRefClassify()) and therefore must be
# protected from FOCEI's internal eta-drift reset mechanisms in src/inner.cpp
# (they are only ever updated by the restart-loop's linear-model step).

#' @export
rxUiGet.foceiMuCovEtaVector <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  .i2 <- .iniDf[-.w, ]
  # Only mu-referenced-FOCEI-family methods (muModel != "none") protect
  # mu-ref-covariate etas from the drift-reset mechanisms; every other
  # method (focei/foce/fo/foi/agq/laplace/etc, muModel="none" default) must
  # see the same all-zero vector it always has, so behavior is unchanged.
  .muModel <- rxode2::rxGetControl(.ui, "muModel", "none")
  if (length(.i2$name) > 0 && !identical(.muModel, "none")) {
    .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
    .i2 <- .i2[order(.i2$neta1), ]
    .muCovEtas <- .muRefClassify(.ui)$muCovEtas
    vapply(seq_along(.i2$neta1), function(i) {
      if (.i2$name[i] %in% .muCovEtas) 1L else 0L
    }, integer(1))
  } else {
    integer(0)
  }
}
attr(rxUiGet.foceiMuCovEtaVector, "rstudio") <- c(0L, 1L)

#' @export
rxUiGet.foceiSkipCov <- function(x, ...) {
  .ui <- x[[1]]
  .maxTheta <- max(.ui$iniDf$ntheta, na.rm=TRUE)
  if (!is.finite(.maxTheta)) {
    logical(0)
  } else {
    .theta <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), ]
    .skipCov <- rep(FALSE, .maxTheta)
    .skipCov[which(!is.na(.theta$err))] <- TRUE
    .skipCov[.theta$fix] <- TRUE
    if (length(.uiIovEnv$iovVars) > 0) {
      .skipCov[which(.theta$name %in% .uiIovEnv$iovVars)] <- TRUE
    }
    # Mixture probability parameters are estimated on the mlogit scale; their
    # covariance cannot be meaningfully interpreted, so skip them.
    if (length(.ui$mixProbs) > 0) {
      .skipCov[which(.theta$name %in% .ui$mixProbs)] <- TRUE
    }
    .skipCov
  }
}
#attr(rxUiGet.foceiSkipCov, "desc") <- "what covariance elements to skip"
attr(rxUiGet.foceiSkipCov, "rstudio") <- c(FALSE, TRUE)

#'  Setup the skip covariate function
#'
#'
#' @param ui rxode2 parsed function
#' @param env environment
#' @return Nothing called for side effects.
#' @author Matthew L. Fidler
#' @noRd
.foceiSetupSkipCov <- function(ui, env) {
  env$skipCov <- rxode2::rxGetControl(ui, "skipCov", NULL)
  if (is.null(env$skipCov)) {
    env$skipCov <- ui$foceiSkipCov
  }
  .maxTheta <- max(ui$iniDf$ntheta, na.rm=TRUE)
  if (!is.finite(.maxTheta)) {
    .maxTheta <- 0
  }
  if (length(env$skipCov) > .maxTheta) {
    # The covariance step writes a full theta+Omega skip mask back to the fit
    # env (#694: the Omega tail is no longer all-TRUE); the R-side skipCov is
    # theta-only, and foceiCalcCov rebuilds the full mask from thetaFixed, so
    # drop the trailing Omega entries here unconditionally.
    assign("skipCov", env$skipCov[seq_len(.maxTheta)], env)
  }
  assign("nEstOmega", length(which(!is.na(ui$iniDf$neta1) & !ui$iniDf$fix)),
         env)
  if (length(env$skipCov) != .maxTheta) {
    .iniTheta <- ui$iniDf[!is.na(ui$iniDf$ntheta), ]
    env$skipCov <- is.na(.iniTheta$err)
    warning("'skipCov' improperly specified, reset", call.=FALSE)
  }
}

.foceiOptEnvLik <- function(ui, env) {
  #if (!exists("noLik", envir = env)){
  if (!exists("model", envir=env)) {
    env$model <- rxUiGet.foceiModel(list(ui))
  }
  #} else {
  #env$model <- rxUiGet.ebe(list(ui))
  #}
  .foceiOptEnvAssignTol(ui, env)
  .foceiOptEnvAssignNllik(ui, env)
  .foceiOptEnvSetupBounds(ui, env)
  .foceiOptEnvSetupScaleC(ui, env)
  # Theta-side transform codes (log/logit/probit) come from the
  # shared pure inspector .iterPrintXParFromUi.  Default call gives
  # vectors of length ntheta_total in ntheta order — exactly what
  # focei's C-side consumes (both the iteration printer and the
  # final-fit-summary back-transform read xform$xPar / xform$probitIdx
  # by ntheta index, and xform$logitThetaLow/Hi / probitThetaLow/Hi
  # by the codes those arrays encode).  Single transport sub-list,
  # consistent with how saem (.cfg$xform) and nlm (.ctl$iterPrintXform)
  # ship the same data.
  env$xform <- .iterPrintXParFromUi(ui)
  .foceiSetupSkipCov(ui, env)
  env$control <- get("control", envir=ui)
  env$control$nF <- 0
  env$control$printTop <- TRUE
  env
}

#' @export
rxUiGet.foceiOptEnv <- function(x, ...) {
  .x <- x[[1]]
  if (exists("foceiEnv", envir=.x)) {
    .env <- get("foceiEnv", envir=.x)
    rm("foceiEnv", envir=.x)
  } else {
    .env <- new.env(parent=emptyenv())
  }
  .env$etaNames <- rxUiGet.foceiEtaNames(x, ...)
  .env$thetaFixed <- rxUiGet.foceiFixed(x, ...)
  rxode2::rxAssignControlValue(.x, "foceiMuRef", .x$foceiMuRefVector)
  rxode2::rxAssignControlValue(.x, "foceiMuCovEta", .x$foceiMuCovEtaVector)
  # Mu-referenced-FOCEI-family (mufocei/irlsfocei/...): the theta/eta index
  # arrays are purely UI-derived (no dataset needed) and wired here exactly
  # like foceiMuRef/foceiMuCovEta above; the covariate *values* matrix
  # needs the dataset and is wired separately in .foceiFamilyReturn() once
  # env$dataSav exists.
  .muModelStr <- rxode2::rxGetControl(.x, "muModel", "none")
  rxode2::rxAssignControlValue(.x, "foceiMuModel",
                               c(none = 0L, lin = 1L, irls = 2L)[[.muModelStr]])
  if (!identical(.muModelStr, "none")) {
    .muGroupSetup <- .muRefCppGroupSetup(.x)
  } else {
    .muGroupSetup <- list(muGroupTheta = integer(0), muGroupEta = integer(0),
                          muGroupCovStart = integer(0), muGroupCovCount = integer(0),
                          muGroupCovTheta = integer(0), muGroupCovUserFixed = integer(0),
                          muGroupCovBounded = integer(0), muGroupCovNames = character(0))
  }
  rxode2::rxAssignControlValue(.x, "foceiMuGroupTheta", .muGroupSetup$muGroupTheta)
  rxode2::rxAssignControlValue(.x, "foceiMuGroupEta", .muGroupSetup$muGroupEta)
  rxode2::rxAssignControlValue(.x, "foceiMuGroupCovStart", .muGroupSetup$muGroupCovStart)
  rxode2::rxAssignControlValue(.x, "foceiMuGroupCovCount", .muGroupSetup$muGroupCovCount)
  rxode2::rxAssignControlValue(.x, "foceiMuGroupCovTheta", .muGroupSetup$muGroupCovTheta)
  rxode2::rxAssignControlValue(.x, "foceiMuGroupCovUserFixed", .muGroupSetup$muGroupCovUserFixed)
  # Bounded covariate coefficients (Phase 8): excluded from the design
  # matrix like a user-fixed one, but NOT excluded from the outer
  # optimizer's free-parameter set (foceiSetupTheta_()'s
  # isMuGroupSkip skips this array specifically) -- see
  # .muRefGroups()'s docs (R/muRefClassify.R) for the full rationale.
  rxode2::rxAssignControlValue(.x, "foceiMuGroupCovBounded", .muGroupSetup$muGroupCovBounded)
  # Reuse the existing, documented muModelTol/muModelMaxCycles foceiControl()
  # fields (originally written for the superseded R-level restart loop) to
  # bound the in-C++ inner regress/re-optimize cycle (updateMuGroups(),
  # src/inner.cpp) that now runs once per real outer iteration.
  rxode2::rxAssignControlValue(.x, "foceiMuGroupTol",
                               rxode2::rxGetControl(.x, "muModelTol", 1e-3))
  rxode2::rxAssignControlValue(.x, "foceiMuGroupMaxCycles",
                               rxode2::rxGetControl(.x, "muModelMaxCycles", 10L))
  # Stash the covariate names on the ui so .foceiFamilyReturn() can build
  # the values matrix once the dataset is available, without recomputing
  # .muRefCppGroupSetup() a second time.
  assign(".muGroupCovNames", .muGroupSetup$muGroupCovNames, envir = .x)
  .env$adjLik <- rxode2::rxGetControl(.x, "adjLik", TRUE)
  .env$diagXformInv <- c("sqrt" = ".square", "log" = "exp", "identity" = "identity")[rxode2::rxGetControl(.x, "diagXform", "sqrt")]
  .env$thetaNames <- .x$iniDf[!is.na(.x$iniDf$ntheta), "name"]
  # FIXME is ODEmodel needed?
  .env$ODEmodel <- TRUE
  .foceiOptEnvLik(.x, .env)
  .env
}
attr(rxUiGet.foceiOptEnv, "desc") <- "Get focei optimization environment"
attr(rxUiGet.foceiOptEnv, "rstudio") <- emptyenv()

#' This function process the data for use in focei
#'
#' The $origData is the data that is fed into the focei before modification
#' The $dataSav is the data saved for focei
#'
#' @param data Input dataset
#' @param env focei environment where focei family is run
#' @param ui rxode2 ui
#' @param rxControl is the rxode2 control that is used to translate to the modeling dataset
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.foceiPreProcessData <- function(data, env, ui, rxControl=NULL) {
  if (is.null(rxControl)) {
    .env <- nlmixr2global$nlmixrEvalEnv$envir
    if (!is.environment(.env)) {
      .env <- parent.frame(1)
    }
    rxControl <- rxControl()
  }
  if (inherits(data, "data.frame")) {
    env$origData <- as.data.frame(data[, names(data), drop = FALSE])
  } else {
    env$origData <- as.data.frame(data)
  }
  data <- env$origData
  .covNames <- ui$covariates
  colnames(data) <- vapply(names(data), function(x) {
    if (any(x == .covNames)) {
      x
    } else {
      toupper(x)
    }
  }, character(1))
  if (is.null(data$ID)) data$ID <- 1L
  if (is.null(data$EVID) && is.null(data$AMT)) data$EVID <- 0
  if (is.null(data$AMT)) data$AMT <- 0
  checkmate::assert_names(names(data), must.include = c("DV", "TIME"))
  ## Make sure they are all double amounts.
  for (.v in c("DV", "TIME")) {
    data[[.v]] <- as.double(data[[.v]])
  }
  .mod <- rxode2::rxModelVars(paste0(ui$mv0$model["normModel"], "\n", .foceiToCmtLinesAndDvid(ui)))
  .strCmpP <- .mod$strCmpParams
  .strCmpPNames <- tolower(names(.strCmpP))
  .lvls <- NULL
  for (.v in .covNames) {
    .d <- data[[.v]]
    .strCmpIdx <- match(tolower(.v), .strCmpPNames)
    if (!is.na(.strCmpIdx)) {
      .modelLvls <- levels(.strCmpP[[.strCmpIdx]])
      .extraLvls <- sort(setdiff(unique(as.character(.d)), .modelLvls))
      .fullLvls <- c(.modelLvls, .extraLvls)
      if (inherits(.d, "character") || inherits(.d, "factor")) {
        .l <- factor(as.character(.d), levels = .fullLvls)
        data[[.v]] <- .l
        .lvls <- c(.lvls, setNames(list(.fullLvls), .v))
      }
    } else if (inherits(.d, "character")) {
      .l <- factor(.d)
      data[[.v]] <- .l
      .lvls <- c(.lvls, setNames(list(levels(.l)), .v))
    } else if (inherits(.d, "factor")) {
      .lvls <- c(.lvls, setNames(list(levels(.d)), .v))
    }
  }
  data$nlmixrRowNums <- seq_len(nrow(data))
  .keep <- unique(c("nlmixrRowNums", env$table$keep))
  .et <- rxode2::etTrans(inData=data, obj=.mod,
                         addCmt=TRUE, dropUnits=TRUE,
                         keep=unique(c("nlmixrRowNums", env$table$keep)),
                         allTimeVar=TRUE, keepDosingOnly=FALSE,
                         addlKeepsCov = rxControl$addlKeepsCov,
                         addlDropSs = rxControl$addlDropSs,
                         ssAtDoseTime = rxControl$ssAtDoseTime)
  .lst <- attr(class(.et), ".rxode2.lst")
  .keepL <- .lst$keepL[[1]]
  .idLvl <- .lst$idLvl
  .dat <- cbind(as.data.frame(.et), .keepL)
  env$dataSav <- .dat
  env$idLvl <- .idLvl
  env$covLvl <- .lvls
}

.thetaReset <- new.env(parent = emptyenv())
#' Name the Omega/covariance block of the focei `$cov`
#'
#' `foceiCalcCov` delta-transforms the Omega block of `$cov` to the natural
#' variance-covariance scale and stashes the free eta pair indices (lower
#' triangle, a >= b) as `.covOmegaPairs`; the C++ names only the theta + sigma
#' block, so name the trailing Omega rows/cols here as `om.<theta>` (variance) or
#' `cov.<theta>.<theta>` (covariance) via the mu-reference mapping, falling back to
#' the eta name.  Leaves the printed parameter table untouched (it matches
#' population thetas by name).
#' @param .ret focei fit environment
#' @return nothing; mutates the dimnames of `.ret$cov`
#' @author Hidde van de Beek
#' @noRd
.foceiNameOmegaCov <- function(.ret) {
  if (!exists(".covOmegaPairs", envir = .ret, inherits = FALSE)) {
    return(invisible())
  }
  if (!exists("cov", envir = .ret, inherits = FALSE)) {
    return(invisible())
  }
  .cov <- .ret$cov
  if (!is.matrix(.cov)) {
    return(invisible())
  }
  .prs <- get(".covOmegaPairs", envir = .ret)
  .nom <- nrow(.prs)
  if (is.null(.nom) || .nom == 0L || .nom > nrow(.cov)) {
    return(invisible())
  }
  .ui <- rxode2::rxUiDecompress(get("ui", envir = .ret))
  .map <- .foceiEtaThetaMap(.ui)
  .nm <- ifelse(is.na(.map$thetaForEta), .map$etaNames, .map$thetaForEta)   # eta name if not mu-ref'd
  .omNm <- .foceiOmegaCovNames(.prs, .nm)
  .cn <- colnames(.cov)
  if (is.null(.cn)) {
    .cn <- rep("", nrow(.cov))
  }
  .cn[seq.int(nrow(.cov) - .nom + 1L, nrow(.cov))] <- .omNm
  dimnames(.cov) <- list(.cn, .cn)
  .ret$cov <- .cov
  invisible()
}

#' Is a model structurally within the analytic-covariance scope
#'
#' Checks only the UI-derivable (pre-fit) conditions: rxode2 sensitivity support,
#' a supported residual error model, at least one eta, and no fixed structural
#' theta that would break eta indexing.  A non-mu-referenced eta is not
#' disqualifying: the uniform direction engine keeps its ETA_i_ sensitivity
#' direction and names its Omega by the eta.  The remaining condition (the
#' augmented ODE actually integrating) is only knowable post-fit.  Used to decide
#' up front whether to run the analytic engine and skip the finite-difference cov
#' step, or leave FD as the fallback.
#' @param ui an rxode2 UI object
#' @return TRUE if the model is structurally in analytic scope
#' @author Hidde van de Beek
#' @noRd
.foceiAnalyticInScope <- function(ui) {
  tryCatch({
    if (!.hasRxExpandSens2()) {
      return(FALSE)
    }
    if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) {
      return(FALSE)                                    # FO/FOI: first-order marginal, not conditional
    }
    if (identical(as.integer(rxode2::rxGetControl(ui, "interaction", 1L)), 0L)) {
      return(FALSE)                                    # FOCE -> finite-difference fallback (analytic is FOCEI only)
    }
    if (is.null(.foceiAnalyticErrFull(ui))) {
      return(FALSE)
    }
    .map <- .foceiEtaThetaMap(ui)
    if (nrow(.map$etaRows) == 0L) {
      return(FALSE)
    }
    if (any(.iniIsFixed(ui$iniDf, .map$thetaForEta))) {
      return(FALSE)                                     # fixed (mu-ref) structural theta
    }
    TRUE
  }, error = function(e) FALSE)
}

#' Finite-difference fallback when the analytic covariance is unavailable
#'
#' `covType="analytic"` in scope skips the FD cov step pre-fit (analytic first), so if
#' the analytic engine then cannot produce a cov there is none.  This computes the
#' requested finite-difference cov as a zero-iteration re-fit -- the same mechanism as
#' [setCov], with `covType="fd"` so it does not re-enter the analytic path -- at the
#' \emph{fitted} estimates (the caller passes the UI with the fixed effects injected).
#' @param .ret focei fit environment
#' @param ui the fitted rxode2 UI (fixed effects already injected)
#' @return TRUE if a finite-difference covariance was installed, FALSE otherwise
#' @author Hidde van de Beek
#' @noRd
.foceiAnalyticFdFallback <- function(.ret, ui) {
  .control <- tryCatch(.ret$control, error = function(e) NULL)
  .dat <- tryCatch(.ret$origData, error = function(e) NULL)
  if (is.null(.control) || is.null(.dat)) {
    return(FALSE)
  }
  .control$covType <- "fd"                               # do not re-enter the analytic path
  .control$.covAnalyticMode <- NULL
  .control$.covAnalyticFdMethod <- NULL
  .control$covMethod <- rxode2::rxGetControl(ui, ".covAnalyticFdMethod", 1L)
  .control$maxInnerIterations <- 0L
  .control$maxOuterIterations <- 0L
  .control$calcTables <- FALSE
  .control$skipCov <- .ret$skipCov
  .control$etaMat <- .ret$etaMat
  .fit2 <- tryCatch(
    nlmixr2CreateOutputFromUi(ui, data = .dat, control = .control,
                              table = .ret$table, env = new.env(parent = emptyenv()),
                              est = "none"),
    error = function(e) NULL)
  .cov <- tryCatch(.fit2$cov, error = function(e) NULL)
  if (is.null(.cov) || !is.matrix(.cov)) {
    return(FALSE)
  }
  .ret$cov <- .cov
  .ret$covMethod <- .fit2$covMethod
  if (!is.null(.fit2$parFixedDf)) .ret$parFixedDf <- .fit2$parFixedDf
  if (!is.null(.fit2$parFixed)) .ret$parFixed <- .fit2$parFixed
  TRUE
}

#' Post-fit half of the `covType="analytic"` engine
#'
#' `.covAnalyticMode` from [.foceiFamilyControl] drives it: `"compute"` (in analytic
#' scope, the FD cov step was skipped) computes the exact analytic cov and installs it;
#' if the analytic engine cannot produce one (out of scope on the data, or an augmented
#' ODE that will not solve) it falls back to the finite-difference cov
#' ([.foceiAnalyticFdFallback]).  `"fd"` (out of scope up front) means FD already ran,
#' so this just announces it.  The analytic engine [.foceiCov] is a post-fit oracle
#' over a fit object, so a light shim exposes the pieces it needs.
#' @param .ret focei fit environment
#' @return nothing; may set `.ret$cov`
#' @author Hidde van de Beek
#' @noRd
.foceiAnalyticCovOverride <- function(.ret) {
  if (!exists("ui", envir = .ret, inherits = FALSE)) {
    return(invisible())
  }
  .ui <- rxode2::rxUiDecompress(get("ui", envir = .ret))
  .mode <- rxode2::rxGetControl(.ui, ".covAnalyticMode", "")
  if (identical(.mode, "fd")) {                         # out of scope: FD ran as the fallback
    .minfo("analytic covariance out of scope for this model; using finite differences")
    return(invisible())
  }
  if (!identical(.mode, "compute")) {
    return(invisible())                                 # covType="fd", or no cov requested
  }
  .covFull <- isTRUE(rxode2::rxGetControl(.ui, "covFull", FALSE))
  # the analytic engine (and the FD fallback) read the population thetas from
  # ui$iniDf$est, so the UI must carry the fitted estimates (at this point .ret$ui may
  # still hold the initial values); inject the fitted fixed effects.
  .fx <- tryCatch(.ret$fixef, error = function(e) NULL)
  if (!is.null(.fx) && length(.fx)) {
    .idf <- .ui$iniDf
    .mi <- match(names(.fx), .idf$name)
    .ok <- !is.na(.mi)
    .idf$est[.mi[.ok]] <- as.numeric(.fx)[.ok]
    .ui <- tryCatch({
      .ui$iniDf <- .idf
      .ui
    }, error = function(e) .ui)
  }
  # Bounded-parameter transforms (preProcessBoundedTransform) reparametrize thetas onto
  # an unconstrained scale; the FD cov applies the natural-scale Jacobian but the
  # analytic path does not, so fall back to the finite-difference cov (#697 finding 9).
  if (length(.ui$boundedTransforms)) {
    .minfo("analytic covariance does not support bounded-parameter transforms; using finite differences")
    .foceiAnalyticFdFallback(.ret, .ui)
    return(invisible())
  }
  .shim <- list(finalUi = .ui, omega = .ret$omega, eta = .ret$ranef,
                dataSav = .ret$dataSav)
  .ac <- tryCatch(.foceiCov(.shim, covFull = .covFull), error = function(e) NULL)
  if (is.null(.ac)) {
    # analytic could not produce a cov (out of scope on the data, or an augmented ODE
    # that would not solve) -> compute the finite-difference cov now as the fallback
    # (#697 findings 5,12; the FD re-run applies its own boundary guard).
    if (.foceiAnalyticFdFallback(.ret, .ui)) {
      .minfo("analytic covariance not available for this fit; used the finite-difference covariance")
    } else {
      warning(paste0("nlmixr2: covType=\"analytic\" produced no covariance and the ",
                     "finite-difference fallback did not either (e.g. a boundary ",
                     "estimate) for this fit."), call. = FALSE)
    }
    return(invisible())
  }
  .ret$cov <- .ac$cov
  .ret$covMethod <- .ac$method                          # "analytic-fd2" / "analytic"
  assign(".covMethodEngine", paste0("analytic (", .ac$method, ")"), envir = .ret)
  # the FD cov step was skipped, so popDf came back with NA SEs; fill the SE/%RSE/CI
  # columns from the analytic cov so the printed table matches an FD fit.
  .foceiAnalyticFillParFixed(.ret, .ac$se)
  invisible()
}

#' Back-transform a theta from the estimation scale to the natural scale
#'
#' Mirrors the C++ `scaleBackTransform` (src/scale.h): `exp` for a log theta,
#' `expit` for a logit theta, probit-inverse for a probit theta, identity
#' otherwise, using rxode2's `expit` / `probitInv`.
#' @param i the 1-based theta index
#' @param x the value on the estimation scale
#' @param xf the fit env's `xform` list (log/logit/probit theta indices + bounds)
#' @return the back-transformed value on the natural scale
#' @author Hidde van de Beek
#' @noRd
.scaleBackTransformR <- function(i, x, xf) {
  if (is.null(xf)) {
    return(x)
  }
  if (!is.null(xf$logNthetas) && i %in% xf$logNthetas) {
    return(exp(x))
  }
  .w <- match(i, xf$logitNthetas)
  if (!is.na(.w)) {
    return(rxode2::expit(x, xf$logitNthetasLow[.w], xf$logitNthetasHi[.w]))
  }
  .w <- match(i, xf$probitNthetas)
  if (!is.na(.w)) {
    return(rxode2::probitInv(x, xf$probitNthetasLow[.w], xf$probitNthetasHi[.w]))
  }
  x
}

#' Fill the parameter table from an externally-computed covariance
#'
#' For the `covType="analytic"` path (where the analytic cov overwrote the FD one),
#' rebuilds the numeric SE / %RSE / back-transformed CI columns of `popDf` from the
#' analytic covariance and the character `popDfSig` in the same form
#' `foceiFinalizeCov` (src/inner.cpp) builds them for a finite-difference cov, so
#' the printed table is identical.  Runs before [.updateParFixed], which then
#' appends the BSV / shrinkage columns as usual.
#' @param .ret focei fit environment
#' @param seNamed the analytic SE vector named by parameter; only names matching a
#'   theta row are used (the residual sigma SE lands only when `covFull`; Omega SEs
#'   are ignored, as they are not table rows)
#' @return nothing; mutates `.ret$popDf`, `.ret$popDfSig` and `.ret$se`
#' @author Hidde van de Beek
#' @noRd
.foceiAnalyticFillParFixed <- function(.ret, seNamed) {
  if (!exists("popDf", envir = .ret, inherits = FALSE)) {
    return(invisible())
  }
  .pd <- .ret$popDf
  .tn <- rownames(.pd)
  if (is.null(.tn) || is.null(.pd$Estimate)) {
    return(invisible())
  }
  .est <- .pd$Estimate
  .bt <- .pd[["Back-transformed"]]
  .ci <- tryCatch(.ret$control$ci, error = function(e) 0.95)
  if (is.null(.ci)) {
    .ci <- 0.95
  }
  .qn <- stats::qnorm(1 - (1 - .ci) / 2)
  .xf <- if (exists("xform", envir = .ret, inherits = FALSE)) .ret$xform else NULL
  .sdig <- tryCatch(as.integer(.ret$control$sigdig), error = function(e) 3L)
  if (length(.sdig) != 1L || is.na(.sdig)) {
    .sdig <- 3L
  }
  .se <- setNames(rep(NA_real_, length(.tn)), .tn)
  .mi <- match(.tn, names(seNamed))
  .ok <- !is.na(.mi)
  .se[.ok] <- as.numeric(seNamed)[.mi[.ok]]
  .rse <- abs(.se / .est) * 100
  .lo <- .hi <- rep(NA_real_, length(.tn))
  for (.i in seq_along(.tn)) {
    if (is.finite(.se[.i])) {
      .lo[.i] <- .scaleBackTransformR(.i, .est[.i] - .se[.i] * .qn, .xf)
      .hi[.i] <- .scaleBackTransformR(.i, .est[.i] + .se[.i] * .qn, .xf)
    }
  }
  .pd$SE <- .se
  .pd[["%RSE"]] <- .rse
  .pd[["CI Lower"]] <- .lo
  .pd[["CI Upper"]] <- .hi
  .ret$popDf <- .pd
  .ret$se <- .se
  # character popDfSig: Est. | SE | %RSE | Back-transformed(ci%CI), matching the FD form
  .g <- function(v) ifelse(is.finite(v), formatC(v, digits = .sdig, format = "g"), "")
  .btName <- paste0("Back-transformed(", as.integer(.ci * 100), "%CI)")
  .btCi <- vapply(seq_along(.tn), function(.i) {
    .b <- formatC(.bt[.i], digits = .sdig, format = "g")
    if (is.finite(.se[.i])) {
      paste0(.b, " (", formatC(.lo[.i], digits = .sdig, format = "g"), ", ",
             formatC(.hi[.i], digits = .sdig, format = "g"), ")")
    } else {
      .b
    }
  }, character(1))
  .sig <- data.frame(`Est.` = formatC(.est, digits = .sdig, format = "g"),
                     `SE` = .g(.se),
                     `%RSE` = .g(.rse),
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
  .sig[[.btName]] <- .btCi
  rownames(.sig) <- .tn
  .ret$popDfSig <- .sig
  invisible()
}

#' Internal focei fit function in R
#'
#' @param .ret Internal focei environment
#' @return Modified focei environment with fit information (from C++)
#' @author Matthew L. Fidler
#' @noRd
.foceiFitInternal <- function(.ret) {
  if (exists("objective", .ret)) {
    checkmate::assertNumeric(.ret$objective, len=1, .var.name="fitEnv$objective")
  }
  if (exists("etaObf", .ret)) {
    checkmate::assertDataFrame(.ret$etaObf, .var.name="fitEnv$etaObf")
    if (!(names(.ret$etaObf)[1] == "ID")) {
      stop("the first column of fitEnv$etaObj needs to be an integer and named ID",
           call.=FALSE)
    }
    if (is.factor(.ret$etaObf$ID)) {
      .ret$etaObf$ID <- as.integer(.ret$etaObf$ID)
    }
    checkmate::assertInteger(.ret$etaObf$ID, any.missing=FALSE, min=1, .var.name="fitEnv$etaObj$ID")
  }
  this.env <- new.env(parent=emptyenv())
  assign("err", "theta reset", this.env)
  ## Event ("jump") sensitivities: when requested, point rxode2's event-
  ## sensitivity globals at the inner (sensitivity) model right before the C++
  ## fit, which solves the inner model through a direct ind_solve() loop (so it
  ## never goes through rxSolve()/.rxSetEventSensDims()).  The handle_evid jump
  ## injection is compartment-count guarded, so the smaller pred model solved in
  ## the same loop skips it safely.  Reset on exit.
  .eventSens <- tryCatch(.ret$control$eventSens, error=function(e) "jump")
  if (identical(.eventSens, "jump") &&
        exists("model", .ret) && !is.null(.ret$model$inner)) {
    .esLoaded <- tryCatch(
      rxode2::rxEventSensLoadModel(.ret$model$inner),
      error=function(e) FALSE)
    if (isTRUE(.esLoaded)) {
      on.exit(rxode2::rxEventSensDeactivate(), add=TRUE)
    }
  }
  .thetaReset$thetaNames <- .ret$thetaNames
  nResets <- 0L
  if (getOption("nlmixr2.retryFocei", TRUE)) {
    while (this.env$err == "theta reset") {
      nResets <- nResets + 1L
      if (nResets > 10L) {
        stop("Maximum number of theta resets (10) exceeded; fit is unstable.", call. = FALSE)
      }
      assign("err", "", this.env)
      .ret0 <- tryCatch(
      {
        foceiFitCpp_(.ret)
      },
      error = function(e) {
        if (regexpr("theta reset", e$message) != -1) {
          assign("zeroOuter", FALSE, this.env)
          assign("zeroGrad", FALSE, this.env)
          if (regexpr("theta reset0", e$message) != -1) {
            assign("zeroGrad", TRUE, this.env)
          }  else if (regexpr("theta resetZ", e$message) != -1) {
            assign("zeroOuter", TRUE, this.env)
          }
          assign("err", "theta reset", this.env)
        } else {
          assign("err", e$message, this.env)
        }
      })
      if (this.env$err == "theta reset") {
        .nm <- names(.ret$thetaIni)
        .ret$thetaIni <- setNames(.thetaReset$thetaIni + 0.0, .nm)
        .ret$rxInv$theta <- .thetaReset$omegaTheta
        .ret$control$printTop <- FALSE
        .ret$etaMat <- .thetaReset$etaMat
        .ret$control$etaMat <- .thetaReset$etaMat
        .ret$control$maxInnerIterations <- .thetaReset$maxInnerIterations
        .ret$control$nF <- .thetaReset$nF
        #.ret$control$gillRetC <- .thetaReset$gillRetC
        #.ret$control$gillRet <- .thetaReset$gillRet
        #.ret$control$gillRet <- .thetaReset$gillRet
        #.ret$control$gillDf <- .thetaReset$gillDf
        #.ret$control$gillDf2 <- .thetaReset$gillDf2
        #.ret$control$gillErr <- .thetaReset$gillErr
        #.ret$control$rEps <- .thetaReset$rEps
        #.ret$control$aEps <- .thetaReset$aEps
        #.ret$control$rEpsC <- .thetaReset$rEpsC
        #.ret$control$aEpsC <- .thetaReset$aEpsC
        .ret$control$c1 <- .thetaReset$c1
        .ret$control$c2 <- .thetaReset$c2
        if (this.env$zeroOuter) {
          message("Posthoc reset")
          warning("Posthoc reset")
          .ret$control$maxOuterIterations <- 0L
        } else if (this.env$zeroGrad && isTRUE(.ret$control$zeroGradBobyqa)) {
          message("Theta reset (zero/bad gradient values); Switch to bobyqa")
          warning("Theta reset (zero/bad gradient values); Switch to bobyqa")
          rxode2::rxReq("minqa")
          .ret$control$outerOptFun <- .bobyqa
          .ret$control$outerOpt <- -1L
          .ret$control$outerOptTxt <- "bobyqa"
        } else {
          message("Theta reset (ETA drift)")
          warning("Theta reset (ETA drift)")
        }
      } else if (this.env$err != "") {
        stop(this.env$err)
      } else {
        return(.ret0)
      }
    }
  } else {
    foceiFitCpp_(.ret)
  }
}

.nlmixrCheckFoceiEnvironment <- function(ret) {
  checkmate::assertDataFrame(ret$dataSav, .var.name="focei$dataSav")
  checkmate::assertNumeric(ret$thetaIni, any.missing=FALSE,
                           null.ok=TRUE, .var.name="focei$thetaIni")
  checkmate::assertLogical(ret$skipCov, null.ok=TRUE,
                           any.missing=FALSE, .var.name="focei$skipCov")
  if (!inherits(ret$rxInv, "rxSymInvCholEnv")) {
    stop("focei$rxInv needs to be of class'rxSymInvCholEnv'",
         call.=FALSE)
  }
  checkmate::assertNumeric(ret$lower, null.ok=TRUE,
                           any.missing=FALSE, .var.name="focei$lower")
  checkmate::assertNumeric(ret$upper, null.ok=TRUE,
                           any.missing=FALSE, .var.name="focei$upper")
  if (length(ret$etaMat) == 1L && is.na(ret$etaMat)) {
    ret$etaMat <- NULL
  }
  checkmate::assertMatrix(ret$etaMat, mode="double", null.ok=TRUE,
                          any.missing=FALSE, .var.name="focei$etaMat")
  if (!inherits(ret$control, "foceiControl")) {
    stop("focei$control must be a focei control object",
         call.=FALSE)
  }
}

#'  Restart the estimation if it wasn't successful by moving the parameters (randomly)
#'
#' @param .ret0 Fit
#' @param .ret Input focei environment
#' @param control Control represents the foceiControl to restart the fit
#' @return final focei fit, may still not work
#' @author Matthew L. Fidler
#' @noRd
.nlmixrFoceiRestartIfNeeded <- function(.ret0, .ret, control) {
  .n <- 1
  .est0 <- .ret$thetaIni
  lower <- .ret$lower
  upper <- .ret$upper
  while (inherits(.ret0, "try-error") && control$maxOuterIterations != 0 && .n <= control$nRetries) {
    .draw <- TRUE
    if (isFALSE(control$zeroGradFirstReset) && grepl("bobyqa", attr(.ret0, "condition")$message)) {
      message("Changing to \"bobyqa\"")
      rxode2::rxReq("minqa")
      .ret$control$outerOpt <- -1L
      .ret$control$outerOptFun <- .bobyqa
      .ret$control$outerOptTxt <- "bobyqa"
      .draw <- FALSE
    }
    ## Maybe change scale?
    message(sprintf("Restart %s/%s", .n, control$nRetries))
    if (is.na(control$zeroGradFirstReset) && .n ==control$nRetries) {
      .ret$control$zeroGradFirstReset <- TRUE
    }
    .ret$control$nF <- 0
    .estNew <- .est0 + 0.2 * .n * abs(.est0) * stats::runif(length(.est0)) - 0.1 * .n
    .estNew <- vapply(
      seq_along(.est0),
      function(.i) {
        if (!.draw || .ret$thetaFixed[.i]) {
          .est0[.i]
        } else if (.estNew[.i] < lower[.i]) {
          lower[.i] + (.Machine$double.eps)^(1 / 7)
        } else if (.estNew[.i] > upper[.i]) {
          upper[.i] - (.Machine$double.eps)^(1 / 7)
        } else {
          .estNew[.i]
        }
      }, numeric(1), USE.NAMES=FALSE)
    .ret$thetaIni <- setNames(.estNew, names(.est0))
    .nlmixrCheckFoceiEnvironment(.ret)
    if (getOption("nlmixr2.retryFocei", TRUE)) {
      .ret0 <- try(.foceiFitInternal(.ret))
    } else {
      .ret0 <- .foceiFitInternal(.ret)
    }

    .n <- .n + 1
  }
  .ret0
}
#'  Assign the control to the ui
#'
#' @param env Estimation/output environment
#' @param ... Other arguments
#' @return nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiFamilyControl <- function(env, ..., type="foceiControl") {
  .ui <- get("ui", envir=env)
  .control <- env$control
  if (is.null(.control)) {
    .control <- do.call(type)
  }
  if (!inherits(.control, type)) {
    .control <- do.call(type, .control)
  }
  if (exists("est", envir = env)) {
    .control$est <- env$est
  }
  if (inherits(nlmixr2global$etaMat, "nlmixr2FitCore") &&
        is.null(.control[["etaMat"]])) {
    warning("Passed the initial etas from the last fit",
            call.=FALSE)
    .control[["etaMat"]] <- nlmixr2global$etaMat$etaMat
  }
  # Change control when there is only 1 item being optimized
  .iniDf <- get("iniDf", envir=.ui)
  .est <- .iniDf[!.iniDf$fix,,drop=FALSE]
  if (length(.est$name) == 0L) {
    .etas <- .iniDf[!is.na(.iniDf$neta1),, drop = FALSE]
    if (length(.etas$name) == 0L) {
      stop("no parameters to estimate", call.=FALSE)
    } else {
      .minfo("no population parameters to estimate; changing to a EBE estimation")
      .control$maxOuterIterations <- 0L # no outer optimization
      .control$normType <- 6L #"constant"
      .control$interaction <- 0L # focei
      .control$covMethod <- 0L # ""
      warning("no population parameters to estimate; changing to a EBE estimation",
              call.=FALSE)
    }
  } else if (length(.est$name) == 1L) {
    .minfo("only one parameter to estimate, using stats::optimize")
    .control$outerOpt <- -1L
    .control$outerOptFun <- .optimize
    .control$normType <- 6L #"constant"
    .control$outerOptTxt <- "stats::optimize"
  }
  .optimHess <- any(.ui$predDfFocei$distribution != "norm")
  if (length(.optimHess) != 1) {
    .optimHess <- FALSE
  }
  .control$needOptimHess <- .optimHess
  if (.control$needOptimHess) {
    .control$interaction <- 0L
  }
  # Analytic covariance engine (covType="analytic").  covType="fd" runs the finite-
  # difference cov as requested and is untouched.  covType="analytic": when the model
  # is in analytic scope, zero covMethod so the (expensive) FD cov step is skipped and
  # the exact analytic cov is computed post-fit ([.foceiAnalyticCovOverride]); the
  # requested FD method is remembered so that IF the analytic engine cannot produce a
  # cov (out of scope on the data, or an augmented ODE that will not solve) the FD cov
  # is computed then, as a fallback.  Out of scope up front -> leave covMethod so FD
  # just runs.  A no-cov request (covMethod==0) is left alone.
  if (identical(.control$covType, "analytic") &&
        !is.null(.control$covMethod) && .control$covMethod != 0L) {
    assign("control", .control, envir = .ui)   # so the scope check sees addProp etc.
    if (.foceiAnalyticInScope(.ui)) {
      .control$.covAnalyticMode <- "compute"
      .control$.covAnalyticFdMethod <- .control$covMethod  # remembered for the FD fallback
      .control$covMethod <- 0L                             # skip the FD cov step (analytic first)
    } else {
      .control$.covAnalyticMode <- "fd"
    }
  }
  assign("control", .control, envir=.ui)
}

#' Get the cmt() and dvid() lines
#'
#' @param ui rxode UI
#' @return cmt() and dvid() string
#' @author Matthew L. Fidler
#' @noRd
.foceiToCmtLinesAndDvid <- function(ui) {
  .cmtLines <- ui$cmtLines
  paste(c("", vapply(seq_along(.cmtLines),
                     function(i){deparse1(.cmtLines[[i]])},
                     character(1), USE.NAMES=FALSE),
          deparse1(ui$dvidLine)),
        collapse="\n")
}
#' Calculate the parameter history
#'
#' @param .ret return data
#' @return parameter history data frame
#' @noRd
#' @author Matthew L. Fidler
.parHistCalc <- function(.ret) {
  .tmp <- .ret$parHistData
  .tmp <- .tmp[.tmp$type == "Unscaled", names(.tmp) != "type"]
  .iter <- .tmp$iter
  .tmp <- .tmp[, names(.tmp) != "iter"]
  data.frame(iter = .iter, .tmp, check.names=FALSE)
}

#' Setup the par history information
#'
#' @param .ret Return data
#' @return Nothing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiSetupParHistData <- function(.ret) {
  if (exists("parHistData", envir=.ret)) {
    .ret$parHistData$type <- factor(.ret$parHistData$type,
                                    levels=c("Gill83 Gradient", "Mixed Gradient", "Forward Difference",
                                             "Central Difference", "Scaled", "Unscaled",
                                             "Back-Transformed", "Forward Sensitivity"))
    .ret$parHistData$iter <- as.integer(.ret$parHistData$iter)
    .ret$parHist <- .parHistCalc(.ret)
  }
}
#' Strip fastmatch properties out of matrix dimensions
#'
#' @param mat matrix, data.frame list or other object to process
#' @return matrix with fastmatch attributes removed from dimnames, if
#'   the object is a list of matrices, it also strips the fastmatch
#'   attributes from each matrix
#' @noRd
#' @author Matthew L. Fidler
.stripFastmatchItem <- function(mat) {
  if (inherits(mat, "data.frame")) {
    for (.n in names(mat)) {
      attr(mat[[.n]], ".match.hash") <- NULL
    }
    return(mat)
  }
  if (is.list(mat)) {
    .n <- names(mat)
    return(stats::setNames(lapply(seq_along(.n), function(i) {
      .stripFastmatchItem(mat[[i]])
    }), .n))
  }
  if (is.character(mat)) {
    .ret <- mat
    attr(.ret, ".match.hash") <- NULL
    return(.ret)
  }
  if (!is.matrix(mat)) {
    return(mat)
  }
  .dn <- dimnames(mat)
  attr(.dn[[1]], ".match.hash") <- NULL
  attr(.dn[[2]], ".match.hash") <- NULL
  dimnames(mat) <- .dn
  mat
}

#' Strips fastmatch hash from dimnames
#'
#'
#' @param ret fit environment to modify
#' @return modified fit environment (though since it is in an environment, it is modified in place)
#' @noRd
#' @author Matthew L. Fidler
.stripFastmatchHash <- function(ret) {
  for (v in c("omega", "phiC", "phiH")) {
    if (exists(v, ret)) {
      ret[[v]] <- .stripFastmatchItem(ret[[v]])
    }
  }
  .ui <- ret$ui
  for (v in c("predDf", "muRefDataFrame", "level")) {
    .ui[[v]] <- .stripFastmatchItem(.ui[[v]])
  }
  .ui$control <- NULL
  ret$ui <- .ui
  ret
}

.foceiFamilyReturn <- function(env, ui, ..., method=NULL, est="none") {
  .control <- ui$control
  .control$est <- est
  ui$control <- .control
  .env <- ui$foceiOptEnv
  .env$table <- env$table
  .data <- env$data
  .env$ui <- ui
  .env$est <- est
  if (!is.null(.env$control)) {
    .env$control$est <- est
  }
  nlmixrWithTiming("setup", {
    .foceiPreProcessData(.data, .env, ui, .control$rxControl)
  })
  # Mu-referenced-FOCEI-family (mufocei/irlsfocei/...): the covariate
  # *values* matrix needs the dataset, which only exists after
  # .foceiPreProcessData() populates .env$dataSav -- the index arrays
  # (foceiMuGroupTheta/Eta/CovTheta/...) were already wired in
  # rxUiGet.foceiOptEnv() (UI-only, no dataset needed).
  if (exists(".muGroupCovNames", envir = ui)) {
    .muGroupCovNames <- get(".muGroupCovNames", envir = ui)
    if (length(.muGroupCovNames) > 0L) {
      .ctl <- .env[["control"]]
      .ctl$foceiMuGroupCovData <- .muRefCppCovData(.muGroupCovNames, .env[["dataSav"]])
      .env[["control"]] <- .ctl
    }
  }
  if (!is.null(.env$cov)) {
    if (!checkmate::testMatrix(.env$cov, any.missing=FALSE, min.rows=1, #.var.name="env$cov",
                               row.names="strict", col.names="strict")) {
      .env$covDebug <- .env$cov
      .minfo(paste0("covariance not in proper form, can access value in ", crayon::bold$blue("$covDebug")))
      warning(paste0("covariance not in proper form, can access value in $covDebug"))
      .env$cov <- NULL
    }
  }
  if (.control$nAGQ > 0) {
    .ag <- .agq(length(ui$eta), .control$nAGQ)
    .env$aqn <- as.integer(.ag$n)
    .env$qx <- .ag$x
    .env$qw <- .ag$w
    .env$qfirst <- .ag$first
    .env$nAGQ <- .control$nAGQ
    .env$aqLow <- .control$agqLow
    .env$aqHi <- .control$agqHi
  } else {
    .env$aqn <- 0L
    .env$qx <- double(0)
    .env$qw <- double(0)
    .env$qfirst <- FALSE
    .env$nAGQ <- 0L
    .env$aqLow <- -Inf
    .env$aqHi <- Inf
  }
  # Mu-referenced-FOCEI-family (mufocei/irlsfocei/...): the regression
  # update now runs natively in C++ (updateMuGroups(), src/inner.cpp),
  # driven entirely by the muModel/foceiMuGroup* control values wired in
  # rxUiGet.foceiOptEnv above -- .foceiFitInternal() is called exactly the
  # same way as every other FOCEI-family method, no separate engine.
  if (getOption("nlmixr2.retryFocei", TRUE)) {
    .ret0 <- try(.foceiFitInternal(.env))
  } else {
    .ret0 <- .foceiFitInternal(.env)
  }
  .ret0 <- .nlmixrFoceiRestartIfNeeded(.ret0, .env, .control)
  if (inherits(.ret0, "try-error")) {
    stop("Could not fit data\n  ", attr(.ret0, "condition")$message, call.=FALSE)
  }
  .ret <- nlmixrWithTiming("postprocess", {
    .ret <- .ret0
    if (!is.null(method))
      .ret$method <- method
    .priorEnvTolFactor <- NULL
    if (is.environment(ui) && exists("foceiEnv", envir=ui, inherits=FALSE)) {
      .priorEnv <- ui$foceiEnv
      if (is.environment(.priorEnv) && exists("tolFactor", envir=.priorEnv, inherits=FALSE)) {
        .priorEnvTolFactor <- .priorEnv$tolFactor
      }
    }
    if (exists("ui", envir=.ret)) {
      ui <- rxode2::rxUiDecompress(get("ui", envir=.ret))
    } else {
      ui <- rxode2::rxUiDecompress(ui)
    }
    if (exists(".predDfFocei", envir=ui)) {
      rm(".predDfFocei", envir=ui)
    }
    ui <- rxode2::rxUiCompress(ui)
    .ret$ui <- ui
    .foceiSetupParHistData(.ret)
    # For mixture models: fix ranef (remove MIXEST), build mixList and mixNum
    .mixFix(.ret, ui)
    if (!all(is.na(ui$iniDf$neta1))) {
      if (exists("etaExpected", envir=.ret)) {
        .etas <- .ret$etaExpected
      } else {
        .etas <- .ret$ranef
      }
      .w <- which(names(.etas) %in% c("mixnum", "MIXEST"))
      if (length(.w) > 0L) {
        .etas <- .etas[, -.w, drop=FALSE]
      }
      .thetas <- .ret$fixef
      .pars <- .Call(`_nlmixr2est_nlmixr2Parameters`, .thetas, .etas)
      .ret$shrink <- .Call(`_nlmixr2est_calcShrinkOnly`, .ret$omega, .pars$eta.lst, length(.etas$ID))
    }
    assign("est", est, envir=.ret)
    .foceiNameOmegaCov(.ret)
    .foceiAnalyticCovOverride(.ret)   # covType="analytic": replace FD $cov with the exact analytic cov
    .updateParFixed(.ret)
    if (!exists("table", .ret)) {
      .ret$table <- tableControl()
    }
    .nlmixr2FitUpdateParams(.ret)
    .ret$IDlabel <- rxode2::.getLastIdLvl()
    .idLvl <- if (exists("idLvl", envir=.ret)) .ret$idLvl else character(0)
    if (exists("tolFactor", envir=.ret)) {
      .tf <- .ret$tolFactor
      if (length(.tf) == length(.idLvl)) {
        .tf <- setNames(.tf, .idLvl)
      }
      .ret$tolFactor <- .tf
    }
    if (!is.null(.priorEnvTolFactor) && length(.priorEnvTolFactor) == length(.idLvl)) {
      .foceiTf <- if (exists("tolFactor", envir=.ret)) unname(.ret$tolFactor) else rep(1.0, length(.priorEnvTolFactor))
      .ret$tolFactor <- setNames(pmax(.foceiTf, .priorEnvTolFactor), .idLvl)
    }
    if (exists("skipTable", envir=.ret)) {
      if (is.na(.ret$skipTable)) {
      } else if (.ret$skipTable) {
        .control$calcTables <- FALSE
      }
    }
    assign("skipCov", .env$skipCov, envir=.ret)
    nmObjHandleModelObject(.ret$model, .ret)
    nmObjHandleControlObject(get("control", envir=.ret), .ret)
    .ret
  })
  nlmixr2global$currentTimingEnvironment <- .ret # add environment for updating timing info
  if (.control$calcTables) {
    .tmp <- try(addTable(.ret,
                         updateObject="no",
                         keep=.ret$table$keep,
                         drop=.ret$table$drop,
                         table=.ret$table), silent=TRUE)
    if (inherits(.tmp, "try-error")) {
      warning("error calculating tables, returning without table step", call.=FALSE)
    } else {
      .ret <- .mixFixTable(.tmp, .env, ui)
    }
  }
  assign("sessioninfo", .sessionInfo(), envir=.env)
  nlmixrWithTiming("compress", {
    if (exists("saem", .env)) {
      .saem <- get("saem", envir=.env)
      .saemCfg <- attr(.saem, "saem.cfg")
      # Delete unneeded variables
      .saemCfg2 <- list()
      for (.v in c("i1", "nphi1", "nphi0", "N", "ntotal", "ix_endpnt", "y", "nmc", "niter", "opt", "inits", "Mcovariables")) {
        .saemCfg2[[.v]] <- .saemCfg[[.v]]
      }
      attr(.saem, "saem.cfg") <- .saemCfg2
      rm(list="saem", envir=.env)
      .env$saem0 <- .saem
    }
    if (.control$compress) {
      for (.item in c("origData", "parHistData", "phiM")) {
        if (exists(.item, .env)) {
          .obj <- get(.item, envir=.env)
          .size <- utils::object.size(.obj)
          .type <- rxode2::rxGetDefaultSerialize()
          .objC <- switch(.type,
                 qs2 = {
                   qs2::qs_serialize(.obj)
                 },
                 qdata = {
                   qs2::qd_serialize(.obj)
                 },
                 bzip2 = {
                   memCompress(serialize(.obj, NULL), type="bzip2")
                 },
                 xz = {
                   memCompress(serialize(.obj, NULL), type="xz")
                 },
                 base = {
                   serialize(.obj, NULL)
                 },
                 stop("unknown serialization type") # nocov
                 )
          .size2 <- utils::object.size(.objC)
          if (.size2 < .size) {
            .size0 <- (.size - .size2)
            .malert("compress {  .item } in nlmixr2 object, save { .size0 }" )
            assign(.item, .objC, envir=.env)
          }
        }
      }
    }
    for (.item in c("adj", "adjLik", "diagXformInv", "etaMat", "etaNames",
                    "fullTheta", "scaleC", "gillRet", "gillRetC",
                    "xform",
                    "lower", "noLik", "objf", "OBJF",
                    "rxInv", "scaleC", "se", "skipCov", "thetaFixed", "thetaIni", "thetaNames", "upper",
                    "xType", "IDlabel", "ODEmodel", "model",
                    # times
                    "optimTime", "setupTime", "covTime",
                    "parHist", "dataSav", "idLvl", "theta",
                    "missingTable", "missingControl", "missingEst")) {
      if (exists(.item, .env)) {
        rm(list=.item, envir=.env)
      }
    }
    assign("ui", rxode2::rxUiCompress(.env$ui), envir=.env)
  })
  if (any(names(.ret) == "CWRES") && regexpr("^fo", est) == -1) {
    # focei is available; add objective function
    .setOfvFo(.ret, "focei")
  }
  .postFinalObjectHooksRun(.ret)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.focei <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'focei'",
                             .var.name=.ui$modelName)
  if (!rxode2hasLlik()) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'focei'",
                                      .var.name=.ui$modelName)
  }
  .foceiFamilyControl(env, ...)
  on.exit({
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm("control", envir=.ui)
    }
  })
  .ui <- env$ui
  .ret <- .foceiFamilyReturn(env, .ui, ..., est="focei")
  .ret
}
attr(nlmixr2Est.focei, "covPresent") <- TRUE
attr(nlmixr2Est.focei, "unbounded") <- .foUnbounded
attr(nlmixr2Est.focei, "iov") <- TRUE

#' Add objective function line to the return object
#'
#' @param ret Return object
#' @param objDf Objective function data frame to add
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.addObjDfToReturn <- function(ret, objDf) {
  if (inherits(ret, "nlmixr2FitData")) {
    ret <- attr(class(ret), ".foceiEnv")
  }
  .objDf1 <- get("objDf", ret)
  if (any(names(.objDf1) == "Condition#(Cov)")) {
    if (!any(names(objDf) == "Condition#(Cov)")) {
      objDf[["Condition#(Cov)"]] <- NA_real_
    }
  } else if (any(names(objDf) == "Condition#(Cov)")) {
    if (!any(names(.objDf1) == "Condition#(Cov)")) {
      .objDf1[["Condition#(Cov)"]] <- NA_real_
    }
  }
  if (any(names(.objDf1) == "Condition#(Cor)")) {
    if (!any(names(objDf) == "Condition#(Cor)")) {
      objDf[["Condition#(Cor)"]] <- NA_real_
    }
  } else if (any(names(objDf) == "Condition#(Cor)")) {
    if (!any(names(.objDf1) == "Condition#(Cor)")) {
      .objDf1[["Condition#(Cor)"]] <- NA_real_
    }
  }
  assign("objDf", rbind(.objDf1, objDf), envir=ret)
}


#'@rdname nlmixr2Est
#'@export
nlmixr2Est.output <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'output'", .var.name=.ui$modelName)
  if (!rxode2hasLlik()) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'output'", .var.name=.ui$modelName)
  }

  .foceiFamilyControl(env, ...)
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
  rxode2::rxAssignControlValue(.ui, "maxOuterIterations", 0L)
  rxode2::rxAssignControlValue(.ui, "maxInnerIterations", 0L)
  on.exit({
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm("control", envir=.ui)
    }
  })
  if (!exists("est", envir=env)) env$est <- "posthoc"
  .foceiFamilyReturn(env, .ui, ..., est=env$est)
}

#' Create nlmixr output from the UI
#'
#'
#' @param ui This is the UI that will be used for the translation
#' @param data This has the data
#' @param control focei control for data creation
#' @param table Table options
#' @param env Environment setup which needs the following:
#' - `$table` for table options
#' - `$origData` -- Original Data
#' - `$dataSav` -- Processed data from .foceiPreProcessData
#' - `$idLvl` -- Level information for ID factor added
#' - `$covLvl` -- Level information for items to convert to factor
#' - `$ui` for ui object
#' - `$fullTheta` Full theta information
#' - `$etaObf` data frame with ID, etas and OBJI
#' - `$cov` For covariance
#' - `$covMethod` for the method of calculating the covariance
#' - `$adjObf` Should the objective function value be adjusted
#' - `$objective` objective function value
#' - `$extra` Extra print information
#' - `$method` Estimation method (for printing)
#' - `$omega` Omega matrix
#' - `$theta` Is a theta data frame
#' - `$model` a list of model information for table generation.  Needs a `predOnly` model
#' - `$message` Message for display
#' - `$est` estimation method
#' - `$ofvType` (optional) tells the type of ofv is currently being use
#'
#' There are some more details that need to be described here
#'
#' @param est Estimation method
#' @return nlmixr fit object
#' @author Matthew L. Fidler
#' @export
nlmixr2CreateOutputFromUi <- function(ui, data=NULL, control=NULL, table=NULL, env=NULL, est="none") {
  nlmixr2global$finalUiCompressed <- FALSE
  on.exit(nlmixr2global$finalUiCompressed <- TRUE)
  if (inherits(ui, "function")) {
    ui <- rxode2::rxode2(ui)
  }
  if (!inherits(ui, "rxUi")) {
    stop("the first argument needs to be from rxode2 ui", call.=FALSE)
  }
  ui <- rxode2::rxUiDecompress(ui)
  if (inherits(env, "environment")) {
    assign("foceiEnv", env, envir=ui)
  }
  if (!inherits(data, "data.frame")) {
    stop("the 'data' argument must be a data.frame", call.=FALSE)
  }
  .env <- new.env(parent=emptyenv())
  assign("ui", ui, envir=.env)
  .env$data <- data
  .env$control <- control
  .env$table <- table
  .env$est <- est
  class(.env) <- c("output", "nlmixr2Est")
  nlmixr2Est(.env)
}
