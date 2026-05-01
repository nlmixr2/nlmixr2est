# Pre-processing hook: automatic bounded parameter transformation
#
# For estimation methods that operate in unconstrained space, this hook
# automatically transforms bounded parameters by rewriting the model to
# include the appropriate back-transform:
#
#   - Two finite bounds (a, b): param <- expit(param_internal, a, b)
#   - Lower bound only (a, Inf): param <- a + exp(param_internal)
#   - Upper bound only (-Inf, b): param <- b - exp(param_internal)
#
# The model is rewritten before estimation. The estimator sees an
# unconstrained model. Results are back-transformed afterward, including
# Jacobian correction for the covariance matrix.
#
# Author: Hajar Besbassi & Matt Fidler

#' Identify which parameters need a bounded transform
#'
#' Scans the ui's iniDf for theta parameters with finite lower and/or upper
#' bounds that are not fixed and not residual error terms. Skips parameters
#' with mu-referenced transforms (warns that mu-referencing will be lost).
#'
#' @param ui rxode2 ui object
#' @param est estimation routine name
#' @param control control object
#' @return list of transform specifications (one element per bounded parameter)
#' @noRd
#' @author Hajar Besbassi & Matt Fidler
.getBoundedParams <- function(ui, est, control) {
  .iniDf <- ui$iniDf
  .thetaDf <- .iniDf[!is.na(.iniDf$ntheta), ]

  .transforms <- list()
  for (i in seq_len(nrow(.thetaDf))) {
    .name <- .thetaDf$name[i]
    .lo <- .thetaDf$lower[i]
    .hi <- .thetaDf$upper[i]
    .est <- .thetaDf$est[i]

    # Skip fixed params
    if (.thetaDf$fix[i]) next
    # Skip residual error params (have non-NA err column)
    if (!is.na(.thetaDf$err[i])) next
    # Skip synthetic IOV helper thetas; their dedicated back-transform/finalize
    # path is handled in R/iov.R and should not be rewrapped here.
    if (!is.na(.thetaDf$backTransform[i]) &&
          grepl("^nlmixr2iov", .thetaDf$backTransform[i])) next

    .hasLo <- is.finite(.lo)
    .hasHi <- is.finite(.hi)

    if (.hasLo && .hasHi) {
      .eps <- (.hi - .lo) * 1e-6
      .estClamped <- max(.lo + .eps, min(.hi - .eps, .est))
      .transforms[[length(.transforms) + 1]] <- list(
        name = .name,
        internalName = paste0("rxBoundedTr.", .name),
        type = "logit",
        lower = .lo,
        upper = .hi,
        initTrans = log((.estClamped - .lo) / (.hi - .estClamped)),
        initOrig = .est
      )
    } else if (.hasLo && !.hasHi) {
      .val <- .est - .lo
      if (.val <= 0) .val <- 1e-6
      .transforms[[length(.transforms) + 1]] <- list(
        name = .name,
        internalName = paste0("rxBoundedTr.", .name),
        type = "lower_exp",
        lower = .lo,
        upper = .hi,
        initTrans = log(.val),
        initOrig = .est
      )
    } else if (!.hasLo && .hasHi) {
      .val <- .hi - .est
      if (.val <= 0) .val <- 1e-6
      .transforms[[length(.transforms) + 1]] <- list(
        name = .name,
        internalName = paste0("rxBoundedTr.", .name),
        type = "upper_exp",
        lower = .lo,
        upper = .hi,
        initTrans = log(.val),
        initOrig = .est
      )
    }
  }
  .transforms
}

.warnOnLostBoundedMuRef <- function(ui, newUi, transforms, est, control) {
  if (!.isMuMethod(est, control)) {
    nlmixr2global$transformMu <- FALSE
    return(invisible(NULL))
  }
  .oldCe <- ui$muRefCurEval
  .newCe <- newUi$muRefCurEval
  .lostMuRef <- FALSE
  for (.tr in transforms) {
    .wOld <- which(.oldCe$parameter == .tr$name)
    if (length(.wOld) != 1L) next
    .oldEval <- .oldCe$curEval[.wOld]
    if (!.isCurEvalEncodedFunction(.oldEval) || nchar(.oldEval) == 0L) next
    .wNew <- which(.newCe$parameter == .tr$internalName)
    if (length(.wNew) == 1L &&
          .isCurEvalEncodedFunction(.newCe$curEval[.wNew]) &&
          nchar(.newCe$curEval[.wNew]) > 0L) {
      next
    }
    .lostMuRef <- TRUE
    warning(" mu-reference transform (", .oldEval,
            ") for `", .tr$name, "` lost since bounded (and performance degraded)",
            call. = FALSE)
  }
  nlmixr2global$transformMu <- .lostMuRef
  invisible(NULL)
}

#' Build the transform expression for a parameter
#'
#' Returns a language object for the back-transform line to prepend to the
#' model block, e.g. \code{td1 <- expit(rxBoundedTr.td1, 0, 1)}.
#'
#' @param tr transform specification (named list from .getBoundedParams)
#' @return language object
#' @noRd
#' @author Hajar Besbassi
.buildTransformExpr <- function(tr) {
  switch(tr$type,
    "logit" = str2lang(paste0(tr$name, " <- expit(", tr$internalName,
                              ", ", tr$lower, ", ", tr$upper, ")")),
    "lower_exp" = str2lang(paste0(tr$name, " <- ", tr$lower,
                                  " + exp(", tr$internalName, ")")),
    "upper_exp" = str2lang(paste0(tr$name, " <- ", tr$upper,
                                  " - exp(", tr$internalName, ")"))
  )
}

#' Rewrite the model UI with bounded-parameter transforms injected
#'
#' For each bounded parameter, renames it to \code{rxBoundedTr.<name>}
#' in iniDf, drops its bounds, and prepends a back-transform line to the
#' model block. Attaches the transform list and original UI state to the
#' returned ui so the post-estimation hook can restore natural-scale
#' values and parameter names.
#'
#' @param ui rxode2 ui object
#' @param transforms list of transform specifications
#' @return rewritten rxode2 ui with \code{boundedTransforms}
#' @noRd
#' @author Hajar Besbassi
.rewriteModelWithTransforms <- function(ui, transforms) {
  # Save original UI for post-estimation hook
  # Store in a new environment attached to the new UI so the post-estimation
  # hook can restore original parameter names and apply Jacobian corrections.
  .origEnv <- new.env(parent = emptyenv())
  .origEnv$iniDf <- ui$iniDf
  .origEnv$lstExpr <- ui$lstExpr
  .origEnv$muRefCurEval <- ui$muRefCurEval
  .origEnv$transforms <- transforms

  # Build new ini block: rename params, change bounds, set backTransform
  .iniDf <- ui$iniDf
  .env <- nlmixr2global$nlmixrEvalEnv$envir
  if (is.null(.env)) .env <- globalenv() # fallback
  for (.tr in transforms) {
    .w <- which(.iniDf$name == .tr$name)
    if (length(.w) == 1L) {
      .iniDf$name[.w] <- .tr$internalName
      .iniDf$lower[.w] <- -Inf
      .iniDf$upper[.w] <- Inf
      .iniDf$est[.w] <- .tr$initTrans
    }
  }

  # Build new model block: inject transform lines at the top
  # The transform line creates the original param name as a derived variable:
  #   td1 <- expit(rxBoundedTr.td1, 0, 1)
  # So existing model lines (e.g., f(depot) = td1) still reference td1 correctly.
  # We do NOT rename td1 in existing expressions - only in iniDf.
  .transformExprs <- lapply(transforms, .buildTransformExpr)

  .lstExpr <- ui$lstExpr

  # Prepend transform lines
  .newLstExpr <- c(.transformExprs, .lstExpr)

  # Reconstruct the model
  .modelStr <- paste0("model({",
    paste(vapply(.newLstExpr, deparse1, character(1)), collapse = "\n"),
    "})")
  .model <- str2lang(.modelStr)

  .ini <- as.expression(lotri::as.lotri(.iniDf))
  .ini[[1]] <- quote(`ini`)

  .mod <- .getUiFunFromIniAndModel(ui, .ini, .model)
  .newUi <- .mod()

  # Store transform info and original model for post-estimation hook
  .newUi$boundedTransforms <- transforms
  .newUi
}

#' Check if an estimation method is unbounded
#'
#' Uses the \code{"unbounded"} attribute on the \code{nlmixr2Est.<method>}
#' S3 method. This allows external packages (babelmixr2, etc.) to register
#' their own methods as bounded/unbounded without modifying this file.
#'
#' The attribute can be:
#' \itemize{
#'   \item \code{TRUE} - always unbounded
#'   \item \code{FALSE} - handles bounds natively
#'   \item \code{function(control)} returning logical - conditional
#'     (e.g., optim method-dependent, FOCEI outer optimizer-dependent)
#' }
#'
#' If \code{control$boundedTransform} is \code{FALSE}, returns \code{FALSE}
#' so the preprocessing hook is skipped entirely.
#'
#' @param est estimation routine name
#' @param control control object (may contain \code{boundedTransform} flag)
#' @return boolean
#' @noRd
#' @author Hajar Besbassi
.isUnboundedMethod <- function(est, control = NULL) {
  # Allow user to disable bounded-parameter transforms via control option
  if (!is.null(control) && isFALSE(control$boundedTransform)) return(FALSE)
  .v <- as.character(utils::methods("nlmixr2Est"))
  .method <- paste0("nlmixr2Est.", est)
  if (.method %in% .v) {
    .unbounded <- attr(utils::getS3method("nlmixr2Est", est), "unbounded")
    if (is.null(.unbounded)) return(FALSE)
    if (is.function(.unbounded)) return(isTRUE(.unbounded(control)))
    return(isTRUE(.unbounded))
  }
  FALSE
}

#' Pre-processing hook: inject bounded-parameter back-transforms
#'
#' Registered via \code{preProcessHooksAdd()}. Runs before estimation and,
#' for unbounded methods with bounded parameters, returns a rewritten ui
#' with internal unbounded parameters and prepended back-transform lines.
#' Returns \code{NULL} (no-op) for bounded-native methods or models with
#' no bounded parameters.
#'
#' @param ui rxode2 ui object
#' @param est estimation routine name
#' @param data dataset (unused, required by hook signature)
#' @param control control object
#' @return named list with \code{ui} element, or \code{NULL}
#' @noRd
#' @author Hajar Besbassi
.preProcessBoundedTransform <- function(ui, est, data, control) {
  nlmixr2global$preProcessBoundedTransform <- FALSE
  nlmixr2global$postEstimationBoundedTransform  <- FALSE
  nlmixr2global$transformMu <- FALSE

  if (!.isUnboundedMethod(est, control)) return(NULL)


  .transforms <- .getBoundedParams(ui, est, control)
  if (length(.transforms) == 0L) return(NULL)

  nlmixr2global$preProcessBoundedTransform <- TRUE

  .newUi <- .rewriteModelWithTransforms(ui, .transforms)
  .warnOnLostBoundedMuRef(ui, .newUi, .transforms, est, control)

  list(ui = .newUi)
}

# -----------------------------------------------------------------------
# Back-transform helpers
# -----------------------------------------------------------------------

#' Back-transform parameter history (parHist) from internal to natural scale
#'
#' @param parHist data.frame with one row per iteration, columns for parameters.
#'   Transformed columns use the internal name (e.g., \code{rxBoundedTr.td1}).
#' @param transforms list of transform specifications
#' @return modified data.frame with back-transformed values and original names
#' @author Hajar Besbassi
#' @noRd
.backTransformParHist <- function(parHist, transforms) {
  if (is.null(parHist) || length(transforms) == 0) return(parHist)
  for (.tr in transforms) {
    .col <- .tr$internalName
    if (!.col %in% names(parHist)) next
    .vals <- parHist[[.col]]
    parHist[[.col]] <- switch(.tr$type,
      "logit" = rxode2::expit(.vals, .tr$lower, .tr$upper),
      "lower_exp" = .tr$lower + exp(.vals),
      "upper_exp" = .tr$upper - exp(.vals),
      .vals
    )
    # Rename the column to the original parameter name
    names(parHist)[names(parHist) == .col] <- .tr$name
  }
  parHist
}
#' Post estimation back transform the thetaDf
#'
#' @param env environment to back-transform
#' @return nothing, called for side-effects
#' @noRd
#' @author Hajar Besbassi & Matthew L. Fidler
.postEstimationBoundedTransformThetaDf <- function(env, transforms) {
  # --- Back-transform theta ---
  # env$theta is a data.frame with columns: lower, theta, upper, fixed
  .thetaDf <- env$theta
  if (!is.null(.thetaDf) && is.data.frame(.thetaDf)) {
    .thetaNames <- if (exists("thetaNames", envir = env)) env$thetaNames else rownames(.thetaDf)
    for (.tr in transforms) {
      .w <- which(.thetaNames == .tr$internalName)
      if (length(.w) != 1L) next
      .val <- .thetaDf$theta[.w]
      # Back-transform the estimate
      .thetaDf$theta[.w] <- switch(.tr$type,
                                   "logit" = rxode2::expit(.val, .tr$lower, .tr$upper),
                                   "lower_exp" = .tr$lower + exp(.val),
                                   "upper_exp" = .tr$upper - exp(.val),
                                   .val
                                   )
      # Restore original bounds
      .thetaDf$lower[.w] <- .tr$lower
      .thetaDf$upper[.w] <- .tr$upper
      # Restore original parameter name
      .thetaNames[.w] <- .tr$name
    }
    rownames(.thetaDf) <- .thetaNames
    env$theta <- .thetaDf
    env$thetaNames <- .thetaNames
  }
}
#' Back transform cov matrix with Jacobian
#'
#' @param env environment for back-transform
#' @return nothing, called for side-effects
#' @noRd
#' @author Matthew L. Fidler
.postEstimationBoundedTransformJacobian <- function(env, transforms) {
  # --- Jacobian correction on covariance matrix ---
  .thetaDf <- env$theta
  if (exists("cov", envir = env) && !is.null(env$cov) && is.matrix(env$cov)) {
    # Compute Jacobian diagonal from natural-scale values
    # For logit: d(expit)/d(logit) = (natVal - lo) * (hi - natVal) / (hi - lo)
    # For lower_exp: d(lo + exp(x))/dx = natVal - lo
    # For upper_exp: d(hi - exp(x))/dx = -(hi - natVal) [negative sign preserved]
    .jdiag <- rep(1.0, nrow(.thetaDf))
    for (.tr in transforms) {
      .w <- which(env$thetaNames == .tr$name)
      if (length(.w) != 1L) next
      .natVal <- .thetaDf$theta[.w]
      .jdiag[.w] <- switch(.tr$type,
                           "logit" = (.natVal - .tr$lower) * (.tr$upper - .natVal) / (.tr$upper - .tr$lower),
                           "lower_exp" = .natVal - .tr$lower,
                           "upper_exp" = -(.tr$upper - .natVal),
                           1.0
                           )
    }
    # Apply Jacobian: Cov_natural = J * Cov_internal * J'
    # cov may be smaller than theta (fixed params are excluded via skipCov)
    .nCov <- nrow(env$cov)
    .nTheta <- length(.jdiag)
    if (.nCov <= .nTheta && exists("skipCov", envir = env)) {
      .skipCov <- env$skipCov
      .jCov <- numeric(0)
      for (k in seq_along(.jdiag)) {
        if (k <= length(.skipCov) && !.skipCov[k]) {
          .jCov <- c(.jCov, .jdiag[k])
        }
      }
      if (length(.jCov) == .nCov) {
        .J <- diag(.jCov)
        env$cov <- .J %*% env$cov %*% t(.J)
      }
    }
  }
}
#' Post-Estimation transform the UI back to the original
#'
#' @param env environment to back-transform
#' @return nothing, called for side effects
#' @noRd
#' @author Matthew L. Fidler
.postEstimationBoundedTransformUi <- function(env, transforms, ui) {
  # Remove transformation lines in the ui (and also the transformed injection)
  .trLhs <- lapply(seq_along(transforms), function(i) {
    str2lang(transforms[[i]]$name)
  })

  .keep <- which(vapply(ui$lstExpr,
                        function(expr) {
                          if (is.call(expr) && length(expr) >= 3 &&
                                identical(expr[[1]], quote(`<-`))) {
                            .lhs <- expr[[2]]
                            if (any(vapply(seq_along(.trLhs),
                                           function(i) {
                                             identical(.lhs, .trLhs[[i]])
                                           }, logical(1), USE.NAMES = FALSE))) {
                              return(FALSE) # Remove this expression
                            }
                          }
                          TRUE
                        }, logical(1), USE.NAMES = FALSE))
  .lstExpr <- lapply(.keep,
                     function(i) {
                       ui$lstExpr[[i]]
                     })

  .model <- rxode2::as.model(.lstExpr)


  # now restore the iniDf
  .iniDf <- ui$iniDf
  for (.tr in transforms) {
    .w <- which(.iniDf$name == .tr$internalName)
    if (length(.w) != 1L) next
    .val <- .iniDf$est[.w]
    # Back-transform the estimate
    .iniDf$est[.w] <- switch(.tr$type,
                             "logit" = rxode2::expit(.val, .tr$lower, .tr$upper),
                             "lower_exp" = .tr$lower + exp(.val),
                             "upper_exp" = .tr$upper - exp(.val),
                             .val
                             )
    # Restore original bounds
    .iniDf$lower[.w] <- .tr$lower
    .iniDf$upper[.w] <- .tr$upper
    # Restore original parameter name
    .iniDf$name[.w] <- .tr$name
  }

  .ini <- as.expression(lotri::as.lotri(.iniDf))
  .ini[[1]] <- quote(`ini`)
  rm("boundedTransforms", envir=ui$meta)
  .newUi <- .getUiFunFromIniAndModel(ui, .ini, .model)
  .newUi <- .newUi()
  assign("ui", .newUi, envir = env)
  if (nlmixr2global$transformMu) {
    warning("to keep mu-referencing remove bounds or use control=list(boundedTransform=FALSE)",
            call. = FALSE)
  }
}

#' Post-estimation hook: back-transform bounded parameters
#'
#' Registered via \code{preFinalParTableHooksAdd()}. Runs after estimation
#' but before the final parameter table is built (called from C++
#' \code{foceiFinalizeTables} via \code{.preFinalParTableHooksRun}).
#'
#' Modifies \code{env$theta}, \code{env$thetaNames}, and \code{env$cov} to
#' restore natural-scale parameter values, original parameter names, and
#' Jacobian-corrected covariance matrix. Also restores the original ui on
#' the fit object so downstream consumers see the user's model.
#'
#' @param env fit environment containing \code{ui}, \code{theta},
#'   \code{thetaNames}, \code{cov}
#' @return invisible \code{NULL}; called for its side effects on \code{env}
#' @noRd
#' @author Hajar Besbassi & Matt Fidler
.postEstimationBoundedTransform <- function(env) {
  on.exit({
    nlmixr2global$transformMu <- FALSE
  }, add = TRUE)
  .ui <- env$ui
  if (is.null(.ui)) return(invisible(NULL))

  .transforms <- .ui$boundedTransforms
  if (is.null(.transforms) || length(.transforms) == 0) return(invisible(NULL))

  nlmixr2global$postEstimationBoundedTransform <- TRUE

  .postEstimationBoundedTransformThetaDf(env, .transforms)

  .postEstimationBoundedTransformJacobian(env, .transforms)

  if (exists("parHistData", envir = env) && !is.null(env$parHistData)) {
    env$parHistData <- .backTransformParHist(env$parHistData, .transforms)
  }

  .postEstimationBoundedTransformUi(env, .transforms, .ui)

  invisible(NULL)
}
#' Add internal ability to see if the bounded transform hooks ran (for testing purposes, so it isn't exported)
#'
#' @return a named list with two logical elements: \code{pre} is \code{TRUE} if the pre-processing hook ran and injected transforms, and \code{post} is \code{TRUE} if the post-estimation hook ran and applied back-transforms. Both are \code{FALSE} if the hooks were not triggered (e.g., because the method is bounded-native or there are no bounded parameters).
#' @keywords internal
#' @noRd
#' @author Matthew L. Fidler
.testBoundedTransform <- function() {
  c(pre=nlmixr2global$preProcessBoundedTransform,
    post=nlmixr2global$postEstimationBoundedTransform)
}

# -----------------------------------------------------------------------
# Register hooks
# -----------------------------------------------------------------------
preProcessHooksAdd(".preProcessBoundedTransform", .preProcessBoundedTransform)
preFinalParTableHooksAdd(".postEstimationBoundedTransform", .postEstimationBoundedTransform)
