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

#' The FOCEi objective's analytic gradient at a fit's estimate
#'
#' A zero-outer-iteration FOCEi pass at the fit's estimates, optimizing the
#' etas from the fit's, with the analytic outer gradient on.
#' @param fit the finalized fit
#' @return the OFV gradient named by parameter, or `NULL` when the analytic
#'   gradient is out of scope for the model
#' @noRd
.covNewtonGrad <- function(fit) {
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
  .ctl$fast <- TRUE
  .ctl$outerOpt <- "nlminb"
  .ctl$outerOptTxt <- NULL
  .ctl$outerOptFun <- NULL
  # the nested fits reset the mu-referencing global state
  .savedMuRef <- .muRefTrans$cur
  on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
  .ctl <- suppressWarnings(do.call(foceiControl, .ctl))
  .f2 <- suppressMessages(suppressWarnings(.nlmixr2PriorGateBypass(
    nlmixr2(fit, nlme::getData(fit), "focei", control = .ctl)
  )))
  .g <- .foceiGradDirect(.f2)
  if (is.null(names(.g))) {
    return(NULL)
  }
  .g
}

#' The covariance whose inverse stands in for the OFV Hessian
#'
#' A sandwich (`"r,s"`) is not `2 H^-1`, so its R-only covariance is preferred.
#' @param fit the finalized fit
#' @param env the fit environment
#' @return a named covariance matrix, or `NULL`
#' @noRd
.covNewtonCov <- function(fit, env) {
  .r <- get0("covR", envir = env, inherits = FALSE)
  if (is.matrix(.r) && !is.null(rownames(.r))) {
    return(.r)
  }
  tryCatch(fit$cov, error = function(e) NULL)
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
  # a nested nlmixr2() resets nlmixr2global, including the outer fit's timing
  .savedGlobal <- as.list(nlmixr2global, all.names = TRUE)
  .g <- tryCatch(.covNewtonGrad(fit), error = function(e) NULL)
  rm(list = ls(nlmixr2global, all.names = TRUE), envir = nlmixr2global)
  for (.gn in names(.savedGlobal)) {
    assign(.gn, .savedGlobal[[.gn]], envir = nlmixr2global)
  }
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
