#' Pre-processing hook: automatic bounded parameter transformation
#'
#' For estimation methods that operate in unconstrained space,
#' this hook automatically transforms bounded parameters by rewriting
#' the model to include the appropriate back-transform:
#'
#' - Two finite bounds (a, b): param <- expit(param_internal, a, b)
#' - Lower bound only (a, Inf): param <- a + exp(param_internal)
#' - Upper bound only (-Inf, b): param <- b - exp(param_internal)
#'
#' The model is rewritten before estimation. The estimator sees an
#' unconstrained model. Results are back-transformed afterward,
#' including Jacobian correction for the covariance matrix.
#'
#' @author Hajar Besbassi
#' @noRd

# -----------------------------------------------------------------------
# Identify which parameters need transformation
# -----------------------------------------------------------------------
.getBoundedParams <- function(ui) {
  .iniDf <- ui$iniDf
  .thetaDf <- .iniDf[!is.na(.iniDf$ntheta), ]

  # Check which params already have a mu-ref transform
  .ce <- ui$muRefCurEval

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
    # Skip if already mu-referenced with any transform
    .w <- which(.ce$parameter == .name)
    if (length(.w) == 1L && nchar(.ce$curEval[.w]) > 0) next

    .hasLo <- is.finite(.lo)
    .hasHi <- is.finite(.hi)

    if (.hasLo && .hasHi) {
      .eps <- (.hi - .lo) * 1e-6
      .estClamped <- max(.lo + .eps, min(.hi - .eps, .est))
      .transforms[[length(.transforms) + 1]] <- list(
        name = .name,
        internalName = paste0(.name, "_untransformed"),
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
        internalName = paste0(.name, "_untransformed"),
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
        internalName = paste0(.name, "_untransformed"),
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

# -----------------------------------------------------------------------
# Build the transform expression for a parameter
# -----------------------------------------------------------------------
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

# -----------------------------------------------------------------------
# Rewrite the model UI with transforms injected
# -----------------------------------------------------------------------
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
      # Register back-transform function for results display
      .btFunName <- paste0(".nlmixr2BT_", .tr$internalName)
      .btFun <- switch(.tr$type,
        "logit" = local({
          lo <- .tr$lower; hi <- .tr$upper
          function(x) lo + (hi - lo) / (1 + exp(-x))
        }),
        "lower_exp" = local({
          lo <- .tr$lower
          function(x) lo + exp(x)
        }),
        "upper_exp" = local({
          hi <- .tr$upper
          function(x) hi - exp(x)
        })
      )
      assign(.btFunName, .btFun, envir = .env)
      .iniDf$backTransform[.w] <- .btFunName
    }
  }

  # Build new model block: inject transform lines at the top
  # The transform line creates the original param name as a derived variable:
  #   td1 <- expit(td1_untransformed, 0, 1)
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
  .newUi$boundedOriginal <- .origEnv
  .newUi
}

# -----------------------------------------------------------------------
# Check if estimation method is unbounded
#
# Uses the "unbounded" attribute on the nlmixr2Est.<method> S3 method.
# This allows external packages (babelmixr2, etc.) to register their
# own methods as bounded/unbounded without modifying this file.
#
# The attribute can be:
#   TRUE  - always unbounded
#   FALSE - handles bounds natively
#   function(control) returning logical - conditional (e.g., optim
#     method-dependent, FOCEI outer optimizer-dependent)
# -----------------------------------------------------------------------
.isUnboundedMethod <- function(est, control = NULL) {
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

# -----------------------------------------------------------------------
# The pre-processing hook function
# -----------------------------------------------------------------------
.preProcessBoundedTransform <- function(ui, est, data, control) {
  if (!.isUnboundedMethod(est, control)) return(NULL)

  .transforms <- .getBoundedParams(ui)
  if (length(.transforms) == 0) return(NULL)

  .newUi <- .rewriteModelWithTransforms(ui, .transforms)
  list(ui = .newUi)
}

# -----------------------------------------------------------------------
# Back-transform helpers
# -----------------------------------------------------------------------

#' Back-transform parameter history (parHist) from internal to natural scale
#'
#' @param parHist data.frame with one row per iteration, columns for parameters.
#'   Transformed columns use the internal name (e.g., \code{td1_untransformed}).
#' @param transforms list of transform specifications
#' @return modified data.frame with back-transformed values and original names
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

# -----------------------------------------------------------------------
# Post-estimation hook: back-transform bounded parameters
#
# Registered via preFinalParTableHooksAdd(). Runs after estimation
# but before the final parameter table is built (called from C++
# foceiFinalizeTables via .preFinalParTableHooksRun).
#
# Modifies env$theta, env$thetaNames, env$cov to restore natural-scale
# parameter values, original parameter names, and Jacobian-corrected
# covariance matrix.
# -----------------------------------------------------------------------
.postEstimationBoundedTransform <- function(env) {
  .ui <- env$ui
  if (is.null(.ui)) return(invisible(NULL))

  .transforms <- .ui$boundedTransforms
  if (is.null(.transforms) || length(.transforms) == 0) return(invisible(NULL))

  # --- Back-transform theta ---
  # env$theta is a data.frame with columns: lower, theta, upper, fixed
  .thetaDf <- env$theta
  if (!is.null(.thetaDf) && is.data.frame(.thetaDf)) {
    .thetaNames <- if (exists("thetaNames", envir = env)) env$thetaNames else rownames(.thetaDf)
    for (.tr in .transforms) {
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
    env$theta <- .thetaDf
    env$thetaNames <- .thetaNames
  }

  # --- Jacobian correction on covariance matrix ---
  if (exists("cov", envir = env) && !is.null(env$cov) && is.matrix(env$cov)) {
    # Compute Jacobian diagonal from natural-scale values
    # For logit: d(expit)/d(logit) = (natVal - lo) * (hi - natVal) / (hi - lo)
    # For lower_exp: d(lo + exp(x))/dx = natVal - lo
    # For upper_exp: d(hi - exp(x))/dx = -(hi - natVal) [negative sign preserved]
    .jdiag <- rep(1.0, nrow(.thetaDf))
    for (.tr in .transforms) {
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

  # --- Back-transform parameter history ---
  if (exists("parHistData", envir = env) && !is.null(env$parHistData)) {
    env$parHistData <- .backTransformParHist(env$parHistData, .transforms)
  }

  invisible(NULL)
}

# -----------------------------------------------------------------------
# Register hooks
# -----------------------------------------------------------------------
preProcessHooksAdd(".preProcessBoundedTransform", .preProcessBoundedTransform)
preFinalParTableHooksAdd(".postEstimationBoundedTransform", .postEstimationBoundedTransform)
