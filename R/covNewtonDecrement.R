#' Newton decrement of an objective at an estimate
#'
#' `0.5 * g' H^-1 g` with `H^-1 = cov / 2` (the OFV is `-2 log L`): the OFV a
#' quadratic model says the estimate is above its minimum, which is also half
#' the squared distance to that minimum in the metric of `cov`.
#' @param g named OFV gradient
#' @param cov named covariance of the estimate
#' @return the decrement over the parameters in both (the thetas: the gradient
#'   has no Omega element on the covariance's scale); `NA` when there are none
#' @noRd
.covNewtonDecrement <- function(g, cov) {
  if (!is.numeric(g) || is.null(names(g)) || !is.matrix(cov)) {
    return(NA_real_)
  }
  .nm <- intersect(names(g), rownames(cov))
  .nm <- .nm[is.finite(g[.nm]) & is.finite(diag(cov)[.nm])]
  if (length(.nm) == 0L) {
    return(NA_real_)
  }
  .g <- g[.nm]
  .c <- cov[.nm, .nm, drop = FALSE]
  if (!all(is.finite(.c))) {
    return(NA_real_)
  }
  0.25 * drop(crossprod(.g, .c %*% .g))
}

#' Estimation methods whose estimate is not the optimum of a FOCEi-type
#' objective, so the decrement has something to check
#' @param est the estimation method
#' @return logical
#' @noRd
.covNewtonEst <- function(est) {
  .skip <- c(
    "focei",
    "foce",
    "focep",
    "foi",
    "fo",
    "laplace",
    "agq",
    "posthoc",
    "ifocei",
    "ifoce",
    "ifocep",
    "ilaplace",
    "iagq",
    "mfocei",
    "mfoce",
    "mfocep",
    "mlaplace",
    "magq",
    "none",
    "output",
    "predict",
    "simulate",
    "simulation",
    "rxSolve"
  )
  length(est) == 1L && is.character(est) && nzchar(est) && !(est %in% .skip)
}

#' The control of a zero-outer-iteration FOCEi pass at a fit's estimates
#'
#' The etas are optimized starting from the fit's.
#' @param fit the finalized fit
#' @param fast turn the analytic outer gradient on
#' @return a `foceiControl()`, or `NULL`
#' @noRd
.covNewtonCtl <- function(fit, fast) {
  .ctl <- fit$foceiControl
  if (!is.list(.ctl)) {
    return(NULL)
  }
  .ctl$etaMat <- fit$etaMat
  # a sampler's foceiControl skips the inner problem; the gradient needs the modes
  .ctl$maxInnerIterations <- .covInnerIterations(.ctl)
  .ctl$maxOuterIterations <- 0L
  .ctl$calcTables <- FALSE
  .ctl$covMethod <- 0L
  .ctl$compress <- FALSE
  .ctl$print <- 0L
  .ctl$fo <- FALSE
  .ctl$interaction <- TRUE
  .ctl$nAGQ <- 0L
  .ctl$literalFix <- TRUE
  # mu profiling (the imp family's "lin") would drop thetas from the gradient
  .ctl$muModel <- "none"
  # a derivative-free outerOpt would clear `fast`; no outer step is taken anyway
  .ctl$fast <- fast
  .ctl$outerOpt <- "nlminb"
  .ctl$outerOptTxt <- NULL
  .ctl$outerOptFun <- NULL
  suppressWarnings(do.call(foceiControl, .ctl))
}

#' Run a nested FOCEi pass quietly
#' @param ui the model (or fit) to evaluate
#' @param fit the finalized fit, for its data
#' @param ctl a `.covNewtonCtl()` control
#' @return the nested fit
#' @noRd
.covNewtonFocei <- function(ui, fit, ctl) {
  suppressMessages(suppressWarnings(.nlmixr2PriorGateBypass(
    nlmixr2(ui, nlme::getData(fit), "focei", control = ctl)
  )))
}

#' The FOCEi objective's gradient at a fit's estimate
#'
#' The analytic outer gradient when the model is in its scope, otherwise
#' central differences in the free thetas of `cov`.
#' @param fit the finalized fit
#' @param cov the fit's covariance
#' @return the OFV gradient named by parameter, or `NULL`
#' @noRd
.covNewtonGrad <- function(fit, cov) {
  # the nested fits reset the mu-referencing global state
  .savedMuRef <- .muRefTrans$cur
  on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
  .ctl <- .covNewtonCtl(fit, TRUE)
  if (is.null(.ctl)) {
    return(NULL)
  }
  .g <- tryCatch(.foceiGradDirect(.covNewtonFocei(fit, fit, .ctl)), error = function(e) NULL)
  if (!is.null(names(.g))) {
    return(.g)
  }
  .covNewtonGradFd(fit, cov)
}

#' Central-difference FOCEi gradient in the free thetas of a covariance
#' @param fit the finalized fit
#' @param cov the fit's covariance
#' @return the named OFV gradient, or `NULL`
#' @noRd
.covNewtonGradFd <- function(fit, cov) {
  .ui <- fit$finalUi
  .ini <- .ui$iniDf
  .th <- .ini[!is.na(.ini$ntheta) & !.ini$fix, , drop = FALSE]
  .nm <- intersect(.th$name, rownames(cov))
  if (length(.nm) == 0L) {
    return(NULL)
  }
  .ctl <- .covNewtonCtl(fit, FALSE)
  .est <- stats::setNames(.th$est, .th$name)[.nm]
  .ofv <- function(nm, val) {
    tryCatch(
      {
        .u <- suppressMessages(do.call(rxode2::ini, c(list(.ui), stats::setNames(list(val), nm))))
        .covNewtonFocei(.u, fit, .ctl)$objf
      },
      error = function(e) NA_real_
    )
  }
  .f0 <- NULL
  .g <- vapply(
    .nm,
    function(nm) {
      .h <- 1e-3 * max(1, abs(.est[[nm]]))
      .up <- .ofv(nm, .est[[nm]] + .h)
      .dn <- .ofv(nm, .est[[nm]] - .h)
      if (!is.finite(.up) && !is.finite(.dn)) {
        # both rejected (a proportion near both of its bounds): a smaller step
        .h <- .h / 100
        .up <- .ofv(nm, .est[[nm]] + .h)
        .dn <- .ofv(nm, .est[[nm]] - .h)
      }
      if (is.finite(.up) && is.finite(.dn)) {
        return((.up - .dn) / (2 * .h))
      }
      # a probe ini() rejects (a mixture proportion past 1): one-sided difference
      if (is.null(.f0)) {
        .f0 <<- .ofv(nm, .est[[nm]])
      }
      if (is.finite(.up)) (.up - .f0) / .h else (.f0 - .dn) / .h
    },
    double(1)
  )
  stats::setNames(.g, .nm)
}

#' The covariance whose inverse stands in for the OFV Hessian
#'
#' A sandwich (`"r,s"`) or S covariance is not `2 H^-1`, so the R-only
#' covariance is preferred, and an S covariance without one gives `NULL`.
#' @param fit the finalized fit
#' @param env the fit environment
#' @return a named covariance matrix, or `NULL`
#' @noRd
.covNewtonCov <- function(fit, env) {
  # the R-only covariance: the fit's own, or the "r" one setCov() keeps
  .cl <- .covCacheGet(env)
  for (.r in list(get0("covR", envir = env, inherits = FALSE), .cl[["r"]], .cl[["r (full)"]])) {
    if (is.matrix(.r) && !is.null(rownames(.r))) {
      return(.r)
    }
  }
  # an S (score) covariance alone carries no curvature
  .lab <- tryCatch(as.character(fit$covMethod), error = function(e) "")
  if (length(.lab) == 1L && grepl("(^|[^a-z])s([^a-z]|$)", .lab)) {
    return(NULL)
  }
  tryCatch(fit$cov, error = function(e) NULL)
}

#' `.covNewtonGrad()`, keeping the outer fit's globals
#'
#' A nested `nlmixr2()` resets `nlmixr2global`, including the outer fit's
#' timing; it is restored on any exit, an interrupt included.
#' @param fit the finalized fit
#' @param cov the fit's covariance
#' @return the gradient, or `NULL`
#' @noRd
.covNewtonGradKeep <- function(fit, cov) {
  .savedGlobal <- as.list(nlmixr2global, all.names = TRUE)
  on.exit({
    rm(list = ls(nlmixr2global, all.names = TRUE), envir = nlmixr2global)
    for (.gn in names(.savedGlobal)) {
      assign(.gn, .savedGlobal[[.gn]], envir = nlmixr2global)
    }
  })
  tryCatch(.covNewtonGrad(fit, cov), error = function(e) NULL)
}

#' The `covNewtonWarn` threshold of a fit
#' @param fit the finalized fit
#' @return the threshold; `NA` when the method has none or it is disabled
#' @noRd
.covNewtonWarnLimit <- function(fit) {
  .ctl <- tryCatch(fit$control, error = function(e) NULL)
  .v <- if (is.list(.ctl)) .ctl[["covNewtonWarn"]] else NULL
  .v <- suppressWarnings(as.double(.v))
  if (length(.v) != 1L || !is.finite(.v)) {
    return(NA_real_)
  }
  .v
}

#' Warn when a covariance belongs to an estimate that is not stationary
#'
#' For a method whose estimate is not a FOCEi optimum, stores the Newton
#' decrement of the FOCEi objective at the estimate (`fit$covNewtonDecrement`)
#' and warns (collected into `$runInfo`) when it exceeds `covNewtonWarn`.
#' @param fit the finalized fit
#' @param est the estimation method
#' @return the decrement, invisibly (`NA` when not computed)
#' @noRd
.covNewtonWarn <- function(fit, est) {
  .env <- tryCatch(if (is.environment(fit)) fit else fit$env, error = function(e) NULL)
  if (!is.environment(.env) || !.covNewtonEst(est)) {
    return(invisible(NA_real_))
  }
  if (exists("covNewtonDecrement", envir = .env, inherits = FALSE)) {
    return(invisible(.env$covNewtonDecrement))
  }
  .limit <- .covNewtonWarnLimit(fit)
  .cov <- .covNewtonCov(fit, .env)
  .om <- tryCatch(fit$omega, error = function(e) NULL)
  # an IOV fit's omega is a list of matrices
  if (is.na(.limit) || !is.matrix(.cov) || length(.om) == 0L) {
    return(invisible(NA_real_))
  }
  .g <- .covNewtonGradKeep(fit, .cov)
  .d <- .covNewtonDecrement(.g, .cov)
  if (is.na(.d)) {
    return(invisible(NA_real_))
  }
  assign("covNewtonDecrement", .d, envir = .env)
  if (.d > .limit) {
    warning(.covNewtonWarnMsg(.d), call. = FALSE)
  }
  invisible(.d)
}

#' The non-stationary covariance warning, kept to one `$runInfo` line
#' @param d the Newton decrement
#' @return the message
#' @noRd
.covNewtonWarnMsg <- function(d) {
  sprintf("SEs at a non-stationary estimate (decrement %s); polish with focei", format(signif(d, 3)))
}
