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
# Back-transform results after estimation
# -----------------------------------------------------------------------

#' Back-transform a theta vector from internal to natural scale
#'
#' @param theta named numeric vector in transformed space
#' @param transforms list of transform specifications
#' @return named numeric vector on natural scale
#' @noRd
.backTransformTheta <- function(theta, transforms) {
  for (.tr in transforms) {
    # The theta names use the internal name
    .w <- which(names(theta) == .tr$internalName)
    if (length(.w) != 1L) next
    .val <- theta[.w]
    theta[.w] <- switch(.tr$type,
      "logit" = rxode2::expit(.val, .tr$lower, .tr$upper),
      "lower_exp" = .tr$lower + exp(.val),
      "upper_exp" = .tr$upper - exp(.val),
      .val
    )
    names(theta)[.w] <- .tr$name  # restore original name
  }
  theta
}

#' Compute Jacobian diagonal for bounded back-transforms
#'
#' @param theta named numeric vector in transformed space
#' @param transforms list of transform specifications
#' @return numeric vector of Jacobian diagonal entries
#' @noRd
.boundedJacobianDiag <- function(theta, transforms) {
  .jdiag <- rep(1.0, length(theta))
  for (.tr in transforms) {
    .w <- which(names(theta) == .tr$internalName)
    if (length(.w) != 1L) next
    .val <- theta[.w]
    .jdiag[.w] <- switch(.tr$type,
      "logit" = {
        if (.val >= 0) {
          .t <- exp(-.val)
          (.tr$upper - .tr$lower) * .t / (1 + .t)^2
        } else {
          .t <- exp(.val)
          (.tr$upper - .tr$lower) * .t / (1 + .t)^2
        }
      },
      "lower_exp" = exp(.val),      # d/dx (lower + exp(x)) = exp(x)
      "upper_exp" = -exp(.val),     # d/dx (upper - exp(x)) = -exp(x)
      1.0
    )
  }
  .jdiag
}

#' Transform covariance matrix using Jacobian
#'
#' @param theta parameter vector in transformed space
#' @param covMat covariance matrix in transformed space
#' @param transforms list of transform specifications
#' @return covariance matrix on natural scale
#' @noRd
.backTransformCov <- function(theta, covMat, transforms) {
  .jdiag <- .boundedJacobianDiag(theta, transforms)
  .J <- diag(.jdiag)
  .J %*% covMat %*% t(.J)
}

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
# Register the hook
# -----------------------------------------------------------------------
preProcessHooksAdd(".preProcessBoundedTransform", .preProcessBoundedTransform)
