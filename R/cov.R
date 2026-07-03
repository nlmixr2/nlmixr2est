#' Re-run the covariance step post-fit and install the result
#'
#' Shared core of [setCov] and the analytic-covariance FD fallback
#' ([.foceiAnalyticFdFallback]): given a fitted UI, the original data and a control
#' whose `covType`/`covMethod` select the covariance, re-runs
#' [nlmixr2CreateOutputFromUi] with `est="none"` (a zero-iteration re-fit that only
#' recomputes the cov) and installs `$cov` / `$covMethod` / the parameter table into
#' the fit environment.  Does not overwrite `$cov` when the re-run yields no matrix.
#' @param env fit environment to install into
#' @param ui fitted rxode2 UI
#' @param data original data (`getData(fit)` / `fit$origData`)
#' @param control foceiControl with the covariance selection already set
#' @param table table control (may be NULL)
#' @param outEnv environment passed to [nlmixr2CreateOutputFromUi] as its output env
#'   (setCov threads its own so a supplied covariance matrix is carried in)
#' @return TRUE if a covariance matrix was installed, FALSE otherwise
#' @author Hidde van de Beek
#' @noRd
.foceiRecomputeCov <- function(env, ui, data, control, table, outEnv = new.env(parent = emptyenv())) {
  .fit2 <- tryCatch(
    nlmixr2CreateOutputFromUi(ui, data = data, control = control,
                              table = table, env = outEnv, est = "none"),
    error = function(e) NULL)
  .cov <- tryCatch(.fit2$cov, error = function(e) NULL)
  if (is.null(.cov) || !is.matrix(.cov)) {
    return(FALSE)
  }
  env$cov <- .cov
  env$covMethod <- .fit2$covMethod
  if (!is.null(.fit2$parFixedDf)) env$parFixedDf <- .fit2$parFixedDf
  if (!is.null(.fit2$parFixed)) env$parFixed <- .fit2$parFixed
  TRUE
}

.setCov <- function(obj, ...) {
  .pt <- proc.time()
  .env <- obj
  if (rxode2::rxIs(obj, "nlmixr2FitData")) {
    .env <- obj$env
  }
  if (exists("cov", .env)) {
    .cur <- list(.env$cov)
    names(.cur) <- .env$covMethod
    if (exists("covList", .env)) {
      if (is.null(.env$covList[[.env$covMethod]])) {
        .covList <- c(.env$covList, .cur)
      } else {
        .covList <- .env$covList
      }
    } else {
      .covList <- .cur
    }
    assign("covList", .covList, .env)
  }
  .control <- .env$foceiControl
  .control$maxInnerIterations <- 0L
  .control$maxOuterIterations <- 0L
  .control$boundTol <- 0 # turn off boundary
  .control$calcTables <- FALSE
  .lst <- list(...)
  .env2 <- new.env(parent = emptyenv())
  for (.n in names(.lst)) {
    .control[[.n]] <- .lst[[.n]]
  }
  if (!is.null(.lst$covMethod)) {
    if (rxode2::rxIs(.lst$covMethod, "character")) {
      .lst$covMethod <- match.arg(.lst$covMethod, c("r,s", "r", "s"))
      .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
      .control$covMethod <- .covMethodIdx[.lst$covMethod]
    } else if (inherits(.lst$covMethod, "matrix")) {
      .env2$cov <- as.matrix(.lst$covMethod)
      .env2$ui <- obj$ui
      .control$covMethod <- 0L
    } else if (length(.lst$covMethod) == 1) {
      if (.lst$covMethod == "") {
        .control$covMethod <- 0L
      }
    }
  } else if (.control$covMethod == 0L) {
    .control$covMethod <- 1L
  }
  ## covDerivMethod=c("central", "forward"),
  if (!is.null(.lst$hessEps)) {
    .control$hessEps <- .lst$hessEps
    .lst$hessEps <- NULL
  }
  if (!is.null(.lst$gillKcov)) {
    .control$gillKcov <- .lst$gillKcov
    .lst$gillKcov <- NULL
  }
  if (!is.null(.lst$gillStepCov)) {
    .control$gillStepCov <- .lst$gillStepCov
    .lst$gillStepCov <- NULL
  }
  if (!is.null(.lst$gillFtolCov)) {
    .control$gillFtolCov <- .lst$gillFtolCov
    .lst$gillFtolCov <- NULL
  }
  if (!is.null(.lst$rmatNorm)) {
    .control$rmatNorm <- .lst$rmatNorm
    .lst$rmatNorm <- NULL
  }
  if (!is.null(.lst$smatNorm)) {
    .control$smatNorm <- .lst$smatNorm
    .lst$smatNorm <- NULL
  }
  if (!is.null(.lst$covGillF)) {
    .control$covGillF <- .lst$covGillF
    .lst$covGillF <- NULL
  }
  if (!is.null(.lst$covSmall)) {
    .control$covSmall <- .lst$covSmall
    .lst$covSmall <- NULL
  }
  .dat <- getData(obj)
  .ui <- obj$ui
  .control$skipCov <- obj$skipCov
  .control$etaMat <- obj$etaMat # as.matrix(nlme::random.effects(obj)[, -1])
  .foceiRecomputeCov(.env, .ui, .dat, .control, .env$table, .env2)
  #.updateParFixed(.env)
  .parent <- parent.frame(2)
  .bound <- do.call("c", lapply(ls(.parent), function(.cur) {
    if (identical(.parent[[.cur]], obj)) {
      return(.cur)
    }
    return(NULL)
  }))
  message(paste0("Updated original fit object ", ifelse(is.null(.bound), "", crayon::yellow(.bound))))
  .env$time$covariance <- (proc.time() - .pt)["elapsed"]
  return(.env$cov)
}

#' Set the covariance type based on prior calculated covariances
#'
#' @param fit nlmixr2 fit
#' @param method covariance method (see the `covMethod` argument for the control
#'   options for the choices)
#' @return Fit object with covariance updated
#' @author Matt Fidler
#' @seealso \code{\link{foceiControl}()}, \code{\link{saemControl}()}
#' @export
setCov <- function(fit, method) {
  .pt <- proc.time()
  .env <- fit
  if (rxode2::rxIs(fit, "nlmixr2FitData")) {
    .env <- fit$env
  }
  if (method == .env$covMethod) {
    stop("no need to switch covariance methods, already set to '",
      method,
      "'",
      call. = FALSE
    )
  }
  if (exists("covList", .env)) {
    .covList <- .env$covList
    .cov <- .covList[[method]]
    if (!is.null(.cov)) {
      .setCov(fit, covMethod = .cov)
      .env$covMethod <- method
      return(invisible(fit))
    }
  }
  stop("different covariance types have not been calculated",
    call. = FALSE
  )
}

##' @export
getVarCov.nlmixr2FitCore <- function(obj, ...) {
  .env <- obj
  if (rxode2::rxIs(obj, "nlmixr2FitData")) {
    .env <- obj$env
  }
  .force <- FALSE
  .args <- list(...)
  if (!is.null(.args$force)) {
    .force <- .args$force
  }
  if (exists("cov", envir = .env) && !.force) {
    if (rxode2::rxIs(.env$cov, "matrix")) {
      return(.env$cov)
    }
  }
  .setCov(obj, ...)
}

##' @export
getVarCov.nlmixr2FitCoreSilent <- getVarCov.nlmixr2FitCore


.cov2cor <- function(cov) {
  .sd2 <- sqrt(diag(cov))
  .cor <- stats::cov2cor(cov)
  dimnames(.cor) <- dimnames(cov)
  diag(.cor) <- .sd2
  .cor
}
