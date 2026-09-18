#' Options for the finite-difference covariance in setCov()
#'
#' Used by \code{setCov(fit, "r,s")}, \code{"r"} and \code{"s"}.  Each option
#' left \code{NULL} keeps the value the fit was estimated with.
#'
#' @inheritParams foceiControl
#' @return \code{rsControl} object
#' @author Matt Fidler
#' @seealso \code{\link{setCov}()}
#' @examples
#' rsControl(hessEps = 1e-4)
#' @export
rsControl <- function(
  hessEps = NULL,
  gillKcov = NULL,
  gillStepCov = NULL,
  gillFtolCov = NULL,
  covGillF = NULL,
  covSmall = NULL,
  rmatNorm = NULL,
  smatNorm = NULL
) {
  checkmate::assertNumeric(hessEps, lower = 0, len = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertIntegerish(gillKcov, lower = 0, len = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertNumeric(gillStepCov, lower = 1, len = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertNumeric(gillFtolCov, lower = 0, len = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertLogical(covGillF, len = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertNumeric(covSmall, lower = 0, len = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertLogical(rmatNorm, len = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertLogical(smatNorm, len = 1, any.missing = FALSE, null.ok = TRUE)
  .ret <- list(
    hessEps = hessEps,
    gillKcov = gillKcov,
    gillStepCov = gillStepCov,
    gillFtolCov = gillFtolCov,
    covGillF = covGillF,
    covSmall = covSmall,
    rmatNorm = rmatNorm,
    smatNorm = smatNorm
  )
  if (!is.null(.ret$gillKcov)) {
    .ret$gillKcov <- as.integer(.ret$gillKcov)
  }
  .ret <- .ret[!vapply(.ret, is.null, logical(1))]
  class(.ret) <- "rsControl"
  .ret
}

#' Options for the SAEM stochastic-approximation covariance in setCov()
#'
#' Used by \code{setCov(fit, "sa")}, which runs a short SAEM at the fit's
#' estimates before the covariance phase.
#'
#' @param nBurn,nEm warm-up iterations that equilibrate the MCMC chains before
#'   the covariance phase
#' @param nSaCov iterations in the covariance phase; more gives a less noisy
#'   covariance
#' @param seed random seed
#' @return \code{saControl} object
#' @author Matt Fidler
#' @seealso \code{\link{setCov}()}, \code{\link{saemControl}()}
#' @examples
#' saControl(nSaCov = 1000)
#' @export
saControl <- function(nBurn = 100L, nEm = 100L, nSaCov = 500L, seed = 99L) {
  checkmate::assertIntegerish(nBurn, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nEm, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nSaCov, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(seed, len = 1, any.missing = FALSE)
  .ret <- list(nBurn = as.integer(nBurn), nEm = as.integer(nEm), nSaCov = as.integer(nSaCov), seed = as.integer(seed))
  class(.ret) <- "saControl"
  .ret
}

#' Options for the importance-sampling covariance in setCov()
#'
#' Used by \code{setCov(fit, "imp")}, which runs frozen importance-sampling EM
#' iterations at the fit's estimates.
#'
#' @param nIter frozen EM iterations (\code{0} is an E-step-only evaluation)
#' @inheritParams impmapControl
#' @return \code{impCovControl} object
#' @author Matt Fidler
#' @seealso \code{\link{setCov}()}, \code{\link{impmapControl}()}
#' @examples
#' impCovControl(isample = 1000)
#' @export
impCovControl <- function(nIter = 1L, isample = 300L, impSeed = 42L) {
  checkmate::assertIntegerish(nIter, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(isample, lower = 1, min.len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(impSeed, len = 1, any.missing = FALSE)
  .ret <- list(nIter = as.integer(nIter), isample = as.integer(isample), impSeed = as.integer(impSeed))
  class(.ret) <- "impCovControl"
  .ret
}

# storage mode of each rsControl() option, so values from the fit's
# foceiControl (which may round-trip logicals as integers) compare exactly
.rsControlMode <- c(
  hessEps = "double",
  gillKcov = "integer",
  gillStepCov = "double",
  gillFtolCov = "double",
  covGillF = "logical",
  covSmall = "double",
  rmatNorm = "logical",
  smatNorm = "logical"
)

#' The cache key of a covariance method's options
#'
#' \code{setCov()} records the result of this generic for every covariance it
#' computes, and reuses a cached covariance only when the key it asks for is
#' \code{identical()} to the recorded one.  The default key is the control
#' object itself, as a plain list.
#'
#' A package whose covariance depends on more than its control -- a covariance
#' seeded from another covariance on the fit, say -- adds a method for its
#' control class that puts that state in the key, so a change in it recomputes
#' the covariance instead of reinstalling a stale one.  Options that do not
#' change the result (parallel workers, say) can be left out.
#'
#' @param control covariance control object (for example \code{rsControl()})
#' @param fit nlmixr2 fit, or its environment when the key of a built-in
#'   method's defaults is taken
#' @param ... ignored
#' @return named list, compared with \code{identical()}
#' @author Matt Fidler
#' @seealso \code{\link{setCov}()}
#' @examples
#' setCovOptions(saControl(), NULL)
#' @export
setCovOptions <- function(control, fit, ...) {
  UseMethod("setCovOptions")
}

#' @rdname setCovOptions
#' @export
setCovOptions.default <- function(control, fit, ...) {
  if (is.null(control)) {
    return(list())
  }
  unclass(control)
}

#' @rdname setCovOptions
#' @export
setCovOptions.rsControl <- function(control, fit, ...) {
  # entries left NULL are taken from the fit's foceiControl
  .env <- if (is.environment(fit)) fit else .setCovEnv(fit)
  .fc <- if (is.environment(.env)) .env$foceiControl else NULL
  .ret <- lapply(names(.rsControlMode), function(.n) {
    .v <- control[[.n]]
    if (is.null(.v)) {
      .v <- .fc[[.n]]
    }
    if (is.null(.v)) {
      return(NULL)
    }
    switch(.rsControlMode[[.n]], double = as.double(.v), integer = as.integer(.v), logical = as.logical(.v))
  })
  names(.ret) <- names(.rsControlMode)
  .ret
}

#' Options a covariance method is computed with, as a comparable plain list
#'
#' @param env fit environment
#' @param control covariance control object
#' @param fit nlmixr2 fit (the environment when no fit is at hand)
#' @return named list
#' @noRd
.covOptionsResolve <- function(env, control, fit = env) {
  if (is.null(control)) {
    return(list())
  }
  .ret <- setCovOptions(control, fit)
  if (is.null(.ret)) list() else .ret
}

#' Options `setCov()` is asked to compute `method` with
#'
#' The supplied `control`, else the default of the method's own `control`
#' argument.  A method without a `control` argument has no options, so a
#' supplied `control` is ignored for it (as the method itself ignores it).
#' @param fit nlmixr2 fit
#' @param env fit environment
#' @param method covariance-method name
#' @param args list of the arguments passed through `setCov()`'s `...`
#' @return list(options = named list, explicit = whether a used `control` was
#'   supplied)
#' @noRd
.covOptionsRequested <- function(fit, env, method, args) {
  .m <- utils::getS3method("setCov", .covBaseName(method), optional = TRUE)
  if (is.null(.m) || !("control" %in% names(formals(.m)))) {
    return(list(options = list(), explicit = FALSE))
  }
  # match the arguments the way the method will, so a positional control counts
  .cl <- as.call(c(list(as.name("setCov"), fit = quote(fit), method = quote(method)), args))
  .ctl <- as.list(match.call(.m, .cl))$control
  .explicit <- !is.null(.ctl)
  if (!.explicit) {
    .ctl <- eval(formals(.m)$control, list(fit = fit, method = method), environment(.m))
  }
  list(options = .covOptionsResolve(env, .ctl, fit), explicit = .explicit)
}

#' Default options of a built-in covariance method on a fit
#'
#' The finite-difference methods use the fit's own `foceiControl`, the analytic
#' one has no options, and "sa"/"imp" use their control defaults.
#' @param env fit environment
#' @param name covariance-method name
#' @return named list, or `NULL` when not a built-in method
#' @noRd
.covOptionsDefault <- function(env, name) {
  if (nzchar(.covFdType(name))) {
    return(.covOptionsResolve(env, rsControl()))
  }
  switch(
    .covBaseName(name),
    analytic = list(),
    sa = .covOptionsResolve(env, saControl()),
    imp = .covOptionsResolve(env, impCovControl()),
    NULL
  )
}

#' Record the options of every covariance a fresh fit holds
#'
#' Only the estimation-time finite-difference and analytic covariances are
#' recorded; any other unrecorded covariance stays unknown.
#' @param env fit environment
#' @return invisibly `TRUE`
#' @noRd
.covOptionsRecordEstimation <- function(env) {
  .names <- names(.covCacheGet(env))
  if (
    exists("cov", envir = env, inherits = FALSE) &&
      .covIsName(env$covMethod)
  ) {
    .names <- c(env$covMethod, .names)
  }
  for (.n in .names) {
    if (!is.null(env$covOptions[[.n]])) {
      next
    }
    if (nzchar(.covFdType(.n)) || identical(.covBaseName(.n), "analytic")) {
      .covOptionsSet(env, .n, .covOptionsDefault(env, .n))
    }
  }
  invisible(TRUE)
}

#' Options a covariance already on the fit was computed with
#'
#' Unrecorded finite-difference and analytic covariances (a fit from before
#' options were recorded) are taken to use the fit's own settings; any other
#' unrecorded covariance is unknown (`NULL`).
#' @param env fit environment
#' @param name covariance-method name
#' @return named list, or `NULL` when unknown
#' @noRd
.covOptionsRecorded <- function(env, name) {
  .rec <- env$covOptions
  if (!is.null(.rec[[name]])) {
    return(.rec[[name]])
  }
  if (nzchar(.covFdType(name)) || identical(.covBaseName(name), "analytic")) {
    return(.covOptionsDefault(env, name))
  }
  NULL
}

#' Does a covariance on the fit match the requested options?
#' @param env fit environment
#' @param name covariance-method name on the fit
#' @param requested requested options
#' @param explicit whether a `control` was supplied
#' @return single logical
#' @noRd
.covOptionsMatch <- function(env, name, requested, explicit) {
  .rec <- .covOptionsRecorded(env, name)
  # an unrecorded estimation-time covariance used the fit's own settings, which
  # is what a call without a control asks for
  if (is.null(.rec)) {
    return(!explicit)
  }
  identical(.rec, requested)
}

#' Record the options a covariance was computed with
#' @param env fit environment
#' @param name covariance-method name
#' @param options options list
#' @return invisibly `TRUE`
#' @noRd
.covOptionsSet <- function(env, name, options) {
  .rec <- env$covOptions
  if (is.null(.rec)) {
    .rec <- list()
  }
  .rec[[name]] <- if (is.null(options)) list() else options
  assign("covOptions", .rec, envir = env)
  invisible(TRUE)
}
