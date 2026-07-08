#' Linearized `(X'X)^-1` covariance, degrading only rank-deficient
#' parameter(s) instead of erroring on the whole matrix
#'
#' @param X n x p design/derivative matrix (rows = stacked per-subject
#'   observation blocks, columns = parameters being estimated)
#' @param tol numeric tolerance passed to `qr()` for rank detection (default
#'   matches `qr()`'s own default of `1e-07`)
#' @return p x p matrix; well-identified parameters get the usual `(X'X)^-1`
#'   covariance, rank-deficient ones get `NA_real_` rows/columns instead of
#'   the whole calculation erroring.
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

#' `chol()` on the finite submatrix of a covariance matrix that may have
#' `NA` rows/columns (see `.nlmixr2RobustCov()`)
#'
#' @param covm p x p covariance matrix, possibly with `NA_real_` rows/columns
#' @return `try()`-wrapped `chol()` result on the finite submatrix (or the
#'   whole matrix if none of it is `NA`)
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

#' Recompute a mu-referenced (lin/irls) fit's covariance on the full base model
#'
#' The estimation-time covariance step bails for mu-referenced-FOCEI-family fits
#' (muModel = "lin"/"irls") -- see `foceiCalcCov()` in `src/inner.cpp` -- because
#' the mu->phi reduced parameterization used during estimation yields incorrect
#' standard errors on the mu-referenced/linear parameters.  This recomputes the
#' covariance at the converged estimates with `muModel = "none"`, i.e. on the full
#' corresponding focei/foce/focep model (all structural thetas as ordinary
#' parameters), and installs it.  No-op for non-mu fits or when no covariance was
#' requested.  Must run before the fit env is compressed (needs `etaMat`).
#' @param .ret assembled fit environment
#' @return invisibly TRUE if the covariance was recomputed and installed
#' @noRd
.foceiRecomputeMuCov <- function(fit, est) {
  # only the mu-referenced (lin/irls) families (mufocei/irlsfocei/mufoce/mufocep/...)
  if (!grepl("^(mu|irls)", est)) return(NULL)
  .control <- tryCatch(fit$foceiControl, error = function(e) NULL)
  if (is.null(.control)) return(NULL)
  .cm <- .control$covMethod
  if (is.null(.cm) || identical(as.integer(.cm), 0L)) return(NULL)
  # deep-copy the UI (an environment) so the nested re-fit cannot mutate THIS fit's UI
  .ui <- tryCatch(rxode2::rxUiDecompress(unserialize(serialize(fit$ui, NULL))),
                  error = function(e) NULL)
  if (is.null(.ui)) return(NULL)
  .baseEst <- sub("^(mu|irls)", "", est)     # mufocei->focei, mufoce->foce, mufocep->focep
  .control$muModel <- "none"                 # recompute on the full (non-reduced) model
  # drop the mu-group wiring so the recompute treats every structural theta as an
  # ordinary parameter (nothing excluded/re-profiled by the mu machinery)
  for (.mn in grep("^foceiMu", names(.control), value = TRUE)) .control[[.mn]] <- NULL
  # Re-fit through the FULL nlmixr2() path (preprocessing hooks matter -- the leaner
  # nlmixr2CreateOutputFromUi posthoc cov is not faithful) with maxOuterIterations=0 so
  # it stays at the converged estimates and only recomputes the covariance.  This must
  # run on the COMPLETED fit (post mu-finalization); mid-assembly it corrupts the outer
  # mu-covariate fit's rewriting state.
  .control$est <- .baseEst
  .control$maxOuterIterations <- 0L
  .control$boundTol <- 0
  .control$calcTables <- FALSE
  .control$skipCov <- NULL                   # recompute skipCov for the full model (keep mu thetas)
  .em <- tryCatch(fit$etaMat, error = function(e) NULL)
  if (!is.null(.em)) .control$etaMat <- .em
  # the nested re-fit resets mu-referencing global state (.muRefTrans$cur); save + restore.
  .savedMuRef <- .muRefTrans$cur
  on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
  .fit2 <- try(suppressMessages(suppressWarnings(
    nlmixr2(.ui, data = getData(fit), est = .baseEst, control = .control))), silent = TRUE)
  if (inherits(.fit2, "try-error") || is.null(.fit2$cov)) return(NULL)
  .env2 <- .fit2$env
  # The base-model re-fit rendered a correct parameter table (SEs computed from its
  # cov) for the SAME parameters at the SAME estimates, so carry its cov + already-
  # rendered tables (popDf/popDfSig/parFixed/se + skipCov and cov diagnostics) over --
  # the mu fit's own popDf$SE is empty (its cov step bailed) and .updateParFixed only
  # reformats it, it does not recompute SEs from the cov matrix.
  .extras <- list()
  for (.n in c("covR", "covS", "covRS", "Rinv", "Sinv", "R", "S", "covLvl", "skipCov",
               "eigenCov", "eigenVecCov", "conditionNumberCov", "covList",
               "popDf", "popDfSig", "parFixedDf", "parFixed", "se")) {
    if (exists(.n, envir = .env2, inherits = FALSE)) .extras[[.n]] <- get(.n, envir = .env2)
  }
  list(cov = .fit2$cov, covMethod = .fit2$covMethod, extras = .extras)
}

#' Install the full-model mu covariance onto a completed mu/irls fit (post-fit).
#' @param fit completed nlmixr2 fit (object or its env)
#' @param est estimation method string
#' @return invisibly TRUE if installed
#' @noRd
.foceiInstallMuCov <- function(fit, est) {
  .r <- tryCatch(.foceiRecomputeMuCov(fit, est), error = function(e) NULL)
  if (is.null(.r)) return(invisible(FALSE))
  .env <- if (is.environment(fit)) fit else tryCatch(fit$env, error = function(e) NULL)
  if (!is.environment(.env)) return(invisible(FALSE))
  assign("cov", .r$cov, envir = .env)
  assign("covMethod", .r$covMethod, envir = .env)
  for (.n in names(.r$extras)) assign(.n, .r$extras[[.n]], envir = .env)
  invisible(TRUE)
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
