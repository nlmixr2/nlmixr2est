#' Compute a covariance matrix from a stacked linearized design matrix,
#' degrading only ill-identified parameter(s) instead of the whole matrix
#'
#' Ordinary `crossprod()`/`backsolve()` covariance (`(X'X)^-1`) errors
#' entirely when any column of `X` is rank-deficient. This detects exactly
#' which column(s) are singular (via `qr()`'s pivoted rank detection) and
#' returns `NA` only for those. Reusable by any mixed-effects method's
#' R-level covariance step; `saem`'s `calc.COV()` is the first consumer.
#'
#' @param X n x p design/derivative matrix (rows = stacked per-subject
#'   observation blocks, columns = parameters being estimated)
#' @param tol numeric tolerance passed to `qr()` for rank detection (default
#'   matches `qr()`'s own default of `1e-07`)
#' @return p x p matrix. Entries are the usual linearized `(X'X)^-1`
#'   covariance for well-identified parameters; the row/column for any
#'   parameter whose contribution to `X` is (numerically) rank-deficient
#'   given the rest is `NA_real_` instead of causing the whole calculation
#'   to error.
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2RobustCov <- function(X, tol = 1e-07) {
  p <- ncol(X)
  qrX <- qr(X, tol = tol)
  rank <- qrX$rank
  if (rank == p) {
    Ri <- backsolve(qr.R(qrX), diag(p))
    return(crossprod(t(Ri)))
  }
  ret <- matrix(NA_real_, p, p)
  if (rank > 0) {
    keep <- sort(qrX$pivot[seq_len(rank)])
    Xkeep <- X[, keep, drop = FALSE]
    qrKeep <- qr(Xkeep, tol = tol)
    Rikeep <- backsolve(qr.R(qrKeep), diag(rank))
    ret[keep, keep] <- crossprod(t(Rikeep))
  }
  ret
}

#' Validate a covariance matrix that may have `NA` rows/columns for
#' ill-identified parameters (see `.nlmixr2RobustCov()`)
#'
#' `chol()` errors on any matrix containing `NA`; this runs it only on the
#' finite submatrix so partial identifiability doesn't discard the rest.
#'
#' @param covm p x p covariance matrix, possibly with `NA_real_` rows/columns
#' @return `try()`-wrapped result of `chol()` on the finite submatrix (or
#'   the whole matrix, if none of it is `NA`); inherits `"try-error"` the
#'   same way a direct `chol(covm)` call would if the finite part itself
#'   isn't valid
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2CholPartial <- function(covm) {
  .bad <- which(is.na(diag(covm)))
  if (length(.bad) == 0) {
    return(try(chol(covm), silent = TRUE))
  }
  .good <- setdiff(seq_len(nrow(covm)), .bad)
  if (length(.good) == 0) {
    return(try(stop("no identifiable parameters"), silent = TRUE))
  }
  try(chol(covm[.good, .good, drop = FALSE]), silent = TRUE)
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
  .mat <- obj$etaMat # as.matrix(nlme::random.effects(obj)[, -1])
  .control$skipCov <- obj$skipCov
  .control$etaMat <- .mat
  .fit2 <- nlmixr2CreateOutputFromUi(.ui, data=.dat, control=.control,
                                     table=.env$table,env=.env2, est="none")
  .env$cov <- .fit2$cov
  .env$parFixedDf <- .fit2$parFixedDf
  .env$parFixed <- .fit2$parFixed
  .env$covMethod <- .fit2$covMethod
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
