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

#' Recompute the covariance eigen diagnostics + objDf condition numbers from the
#' fit's installed covariance (method-independent).  NA rows/columns from a
#' rank-deficient covariance (see `.nlmixr2RobustCov()`) are dropped before the
#' eigen decomposition.  Adds the `Condition#(Cov)`/`Condition#(Cor)` columns
#' when missing so every installed covariance reports its condition numbers.
#' @param env fit environment with `cov` installed
#' @return invisibly TRUE if the diagnostics were updated
#' @noRd
.nlmixr2CovConditionUpdate <- function(env) {
  if (!exists("cov", envir = env, inherits = FALSE)) return(invisible(FALSE))
  .cov <- get("cov", envir = env)
  if (!inherits(.cov, "matrix") || nrow(.cov) == 0L) return(invisible(FALSE))
  .cn <- .cnr <- NA_real_
  .good <- which(!is.na(diag(.cov)))
  if (length(.good) > 0L) {
    .c <- .cov[.good, .good, drop = FALSE]
    .e <- try(eigen(.c, symmetric = TRUE), silent = TRUE)
    if (!inherits(.e, "try-error")) {
      assign("eigenCov", .e$values, envir = env)
      assign("eigenVecCov", .e$vectors, envir = env)
      .a <- abs(.e$values)
      if (length(.a) > 0L) .cn <- max(.a) / min(.a)
    }
    assign("fullCor", .cov2cor(.c), envir = env)
    .er <- try(eigen(stats::cov2cor(.c), symmetric = TRUE), silent = TRUE)
    if (!inherits(.er, "try-error")) {
      assign("eigenCor", .er$values, envir = env)
      assign("eigenVecCor", .er$vectors, envir = env)
      .a <- abs(.er$values)
      if (length(.a) > 0L) .cnr <- max(.a) / min(.a)
    }
  }
  assign("conditionNumberCov", .cn, envir = env)
  assign("conditionNumberCor", .cnr, envir = env)
  if (exists("objDf", envir = env, inherits = FALSE)) {
    .od <- get("objDf", envir = env)
    .od[["Condition#(Cov)"]] <- .cn
    .od[["Condition#(Cor)"]] <- .cnr
    assign("objDf", .od, envir = env)
  }
  invisible(TRUE)
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
      .lst$covMethod <- match.arg(.lst$covMethod, c("analytic", "r,s", "r", "s"))
      if (identical(.lst$covMethod, "analytic")) {
        .control$covMethod <- 2L
        .control$covType <- "analytic"
      } else {
        .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
        .control$covMethod <- .covMethodIdx[.lst$covMethod]
        # an explicit FD request must not re-run the analytic hook
        .control$covType <- "fd"
      }
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
  .nlmixr2CovConditionUpdate(.env)
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

#' Base est whose full model the post-fit covariance recompute runs on
#'
#' `NULL` for ests without a post-fit recompute.
#' @param est estimation method string
#' @return base est string or NULL
#' @noRd
.foceiRecomputeBaseEst <- function(est) {
  # mu-referenced (lin/irls) families recompute on their own base method
  if (grepl("^(m|i)(focei|foce|focep|agq|laplace)$", est)) {
    return(sub("^(m|i)", "", est))
  }
  # methods whose estimation pass cannot produce a covariance (EM table pass,
  # nonparametric engines, or an external engine) recompute on a
  # zero-iteration focei model
  if (est %in% c("imp", "impmap", "qrpem", "nlme") ||
        grepl("^(m|i)?(npag|npb)$", est)) {
    return("focei")
  }
  NULL
}

#' Recompute a fit's covariance on the full base model (post-fit)
#'
#' The estimation-time covariance step bails for mu-referenced-FOCEI-family fits
#' (muModel = "lin"/"irls") -- see `foceiCalcCov()` in `src/inner.cpp` -- because
#' the mu->phi reduced parameterization used during estimation yields incorrect
#' standard errors on the mu-referenced/linear parameters.  The imp family
#' (imp/impmap/qrpem) and nlme never compute a focei covariance during
#' estimation at all.  This recomputes the covariance at the converged
#' estimates with `muModel = "none"`, i.e. on the full base model (all
#' structural thetas as ordinary parameters), and installs it; the C++
#' covariance step runs the usual analytic -> "r,s" -> "r"/"s" chain.  No-op
#' when no covariance was requested.  Must run before the fit env is
#' compressed (needs `etaMat`).
#' @param .ret assembled fit environment
#' @return invisibly TRUE if the covariance was recomputed and installed
#' @noRd
.foceiRecomputeMuCov <- function(fit, est) {
  .baseEst <- .foceiRecomputeBaseEst(est)
  if (is.null(.baseEst)) return(NULL)
  # covMethod="imp" installed the Monte-Carlo importance-sampling covariance;
  # that explicit request wins over the recompute
  if (identical(tryCatch(fit$covMethod, error = function(e) NULL), "imp")) {
    return(NULL)
  }
  .control <- tryCatch(fit$foceiControl, error = function(e) NULL)
  if (is.null(.control)) return(NULL)
  .cm <- .control$covMethod
  if (is.null(.cm) || identical(as.integer(.cm), 0L)) return(NULL)
  # deep-copy the UI (an environment) so the nested re-fit cannot mutate THIS fit's UI
  .ui <- tryCatch(rxode2::rxUiDecompress(unserialize(serialize(fit$ui, NULL))),
                  error = function(e) NULL)
  if (is.null(.ui)) return(NULL)
  .control$muModel <- "none"                 # recompute the foceiModel on the full model
  # drop the mu-group wiring so the recompute treats every structural theta as an
  # ordinary parameter (nothing excluded/re-profiled by the mu machinery); keep `fast`
  # (and every other setting) as specified so the full model is rebuilt the same way.
  for (.mn in grep("^foceiMu", names(.control), value = TRUE)) .control[[.mn]] <- NULL
  # The covariance must be evaluated AT the mu fit's converged point -- NOT re-optimized to
  # a (possibly better) nearby point.  So freeze BOTH problems: maxOuterIterations=0 (final
  # thetas) AND maxInnerIterations=0 (final etas held at etaMat).  Runs through the FULL
  # nlmixr2() path (the leaner nlmixr2CreateOutputFromUi posthoc cov is not faithful); must
  # run on the COMPLETED fit (post mu-finalization) or it corrupts the mu-covariate rewrite.
  .control$est <- .baseEst
  .control$maxOuterIterations <- 0L
  .control$maxInnerIterations <- 0L
  .control$boundTol <- 0
  .control$calcTables <- FALSE
  .control$skipCov <- NULL                   # recompute skipCov for the full model (keep mu thetas)
  # explicitly pin the final thetas (on the UI) and the final etas (etaMat)
  .th <- tryCatch(fit$theta, error = function(e) NULL)
  if (!is.null(.th)) {
    .w <- match(names(.th), .ui$iniDf$name)
    .ok <- !is.na(.w)
    .ui$iniDf$est[.w[.ok]] <- as.numeric(.th)[.ok]
  }
  .eta <- tryCatch(fit$eta, error = function(e) NULL)
  if (!is.null(.eta)) {
    .etaCols <- setdiff(names(.eta), "ID")
    .control$etaMat <- as.matrix(.eta[, .etaCols, drop = FALSE])
  }
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
               "fullCor", "eigenCor", "eigenVecCor", "conditionNumberCor",
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
  .env <- if (is.environment(fit)) fit else tryCatch(fit$env, error = function(e) NULL)
  if (!is.environment(.env)) return(invisible(FALSE))
  if (is.null(.r)) {
    # message only when a recompute was requested but failed and a legacy
    # covariance (e.g. nlme's) is being kept
    if (identical(est, "nlme") &&
          exists("cov", envir = .env, inherits = FALSE)) {
      .cm <- tryCatch(fit$foceiControl$covMethod, error = function(e) 0L)
      if (!is.null(.cm) && !identical(as.integer(.cm), 0L)) {
        message("the analytic/finite-difference covariance could not be computed; keeping the \"nlme\" covariance")
      }
    }
    return(invisible(FALSE))
  }
  # keep the pre-existing covariance (e.g. nlme's tTable cov) recoverable via
  # covList/setCov()
  .stash <- NULL
  if (exists("cov", envir = .env, inherits = FALSE) &&
        exists("covMethod", envir = .env, inherits = FALSE)) {
    .stash <- list(get("cov", envir = .env))
    names(.stash) <- as.character(get("covMethod", envir = .env))
  }
  assign("cov", .r$cov, envir = .env)
  assign("covMethod", .r$covMethod, envir = .env)
  for (.n in names(.r$extras)) assign(.n, .r$extras[[.n]], envir = .env)
  if (!is.null(.stash) && !identical(names(.stash), .r$covMethod)) {
    .covList <- NULL
    if (exists("covList", envir = .env, inherits = FALSE)) {
      .covList <- get("covList", envir = .env)
    }
    if (is.null(.covList[[names(.stash)]])) .covList <- c(.covList, .stash)
    assign("covList", .covList, envir = .env)
  }
  # the mu fit's objDf was rendered before any cov existed; recompute the eigen
  # diagnostics + Condition#(Cov)/Condition#(Cor) from the installed cov
  .nlmixr2CovConditionUpdate(.env)
  invisible(TRUE)
}

#' Set the covariance type based on prior calculated covariances
#'
#' Switches a completed fit's covariance to \code{method}.  A previously
#' computed covariance is re-installed from the cache; otherwise it is
#' recomputed at the converged estimates: \code{"r,s"}/\code{"r"}/\code{"s"} and
#' \code{"analytic"} on a zero-iteration FOCEI model, and \code{"sa"} (SAEM
#' Louis FIM) / \code{"imp"} (importance-sampling Monte-Carlo) via the decoupled
#' recompute engine (the latter two require a mixed-effects fit).  When
#' \code{"sa"}/\code{"imp"}/\code{"analytic"} cannot be computed the covariance
#' is left unchanged (it is never silently downgraded to \code{"r,s"}).
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
  # analytic: compute the EXACT analytic observed information; on any failure the
  # covariance is left UNCHANGED -- it is never silently downgraded to the "r,s"
  # finite-difference covariance (which the C++ cov chain would otherwise fall
  # back to and mislabel "analytic").
  if (identical(method, "analytic")) {
    .r <- tryCatch(.foceiCovAnalyticCalc(fit), error = function(e) NULL)
    if (is.null(.r) || !is.matrix(.r$cov) || !all(is.finite(.r$cov)) ||
          !isTRUE(.covInstallResult(.env, list(cov = .r$cov, covMethod = "analytic")))) {
      stop("covMethod=\"analytic\" could not be computed for this fit; the covariance is left unchanged",
           call. = FALSE)
    }
    assign(".covAnalytic", .r, envir = .env)
    return(invisible(fit))
  }
  # not cached: the finite-difference methods can be recomputed on the full model
  # at the converged estimates (.setCov refits with est="none")
  if (method %in% c("r,s", "r", "s")) {
    .setCov(fit, covMethod = method)
    .env$covMethod <- method
    return(invisible(fit))
  }
  # sa/imp: decoupled recompute at the converged estimates (mixed-effects fits)
  if (method %in% c("sa", "imp")) {
    .r <- tryCatch(.covRecompute(fit, method), error = function(e) NULL)
    if (!isTRUE(.covInstallResult(.env, .r))) {
      stop("covMethod=\"", method, "\" could not be computed for this fit; the covariance is left unchanged",
           call. = FALSE)
    }
    return(invisible(fit))
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
