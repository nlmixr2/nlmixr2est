# fsaemInnerCov.R -- covariate-aware f-SAEM inner model.
#
# For covariate models the SAEM target absorbs the time-invariant covariate
# effect into the per-subject prior mean mprior_i (mprior varies by subject).
# The FOCEi inner centers its eta prior at 0 around a *global* theta, so to
# reproduce that target the inner is built on a phi-collapsed model where each
# mu-referenced parameter's intercept is a per-subject data column
# (phi_j = mpriorData_j[i] + eta_j): the covariate is absorbed into the mprior
# data, the etas are kept (keepEtas), and any time-varying covariate stays as a
# beta regressor.  The mprior data is refreshed each fast iteration; the model
# compiles once (distinct cache entry from the plain focei inner).

#' Build the mprior-as-data inner UI + the mprior column map for a mu-ref model.
#'
#' @return list with `ui` (parsed rxode2 ui using mprior data columns for the
#'   mu-ref intercepts), `mpriorCols` (named char: mu-ref theta -> data column),
#'   and `etaThetas` (mu-ref theta names, in muRefDataFrame/eta order).
#' @noRd
.fsaemInnerMpriorUi <- function(ui) {
  .muRef <- ui$muRefDataFrame
  .etaThetas <- .muRef$theta                       # mu-ref theta per eta (in eta order)
  .thetas <- unique(.etaThetas)
  .mpriorCols <- setNames(paste0("nlmixrMprior", seq_along(.thetas)), .thetas)

  # phi-collapsed model with etas kept (non-time-varying covariates absorbed)
  .lst <- .saemDropMuRefFromModel(ui, keepEtas = TRUE)
  .subst <- function(x) {
    if (is.name(x)) {
      .c <- as.character(x)
      if (.c %in% .thetas) return(as.name(.mpriorCols[[.c]]))
      x
    } else if (is.call(x)) {
      as.call(lapply(x, .subst))
    } else x
  }
  .modelLines <- lapply(.lst, .subst)

  # ini: keep the etas, the residual, and any theta still referenced in the
  # model (time-varying covariate betas); drop the mu-ref intercepts (now data)
  # and the absorbed non-time-varying covariate coefficients.
  .iniDf <- ui$iniDf
  .usedNames <- unique(unlist(lapply(.modelLines, all.vars)))
  .etaRows <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, , drop = FALSE]
  .thetaRows <- .iniDf[!is.na(.iniDf$ntheta), , drop = FALSE]
  .keepTheta <- .thetaRows[
    !(.thetaRows$name %in% .thetas) &                      # not a mu-ref intercept (-> data)
      (!is.na(.thetaRows$err) | .thetaRows$name %in% .usedNames), , drop = FALSE]

  .iniLines <- c(
    vapply(seq_len(nrow(.etaRows)), function(i)
      paste0(.etaRows$name[i], " ~ ", .etaRows$est[i]), character(1)),
    vapply(seq_len(nrow(.keepTheta)), function(i)
      paste0(.keepTheta$name[i], " <- ", .keepTheta$est[i]), character(1)))

  .modelTxt <- vapply(.modelLines, function(e) paste(deparse(e), collapse = " "), character(1))
  .fnTxt <- paste0(
    "function() {\n  ini({\n    ",
    paste(.iniLines, collapse = "\n    "),
    "\n  })\n  model({\n    ",
    paste(.modelTxt, collapse = "\n    "),
    "\n  })\n}")
  .fn <- eval(parse(text = .fnTxt))
  .innerUi <- rxode2::rxUiDecompress(rxode2::rxode2(.fn))
  list(ui = .innerUi, mpriorCols = .mpriorCols, etaThetas = .etaThetas)
}

#' Add/overwrite the per-subject nlmixrMprior* columns in `data`.
#'
#' `mpriorMat` is N x neta (eta/i1 order); column k for eta k is written to the
#' data column for that eta's mu-ref theta, constant within each subject.
#' @noRd
.fsaemSetMpriorData <- function(data, mpriorMat, built) {
  .idn <- if ("ID" %in% names(data)) "ID" else "id"
  .ids <- unique(data[[.idn]])
  .idx <- match(data[[.idn]], .ids)                 # per-row subject index
  for (k in seq_along(built$etaThetas)) {
    .col <- built$mpriorCols[[built$etaThetas[k]]]
    data[[.col]] <- mpriorMat[.idx, k]
  }
  data
}

#' Set up the covariate-aware f-SAEM inner (mprior-as-data model).  Compiles the
#' inner model ONCE (distinct cache entry from the plain focei inner); per-
#' iteration updates only refresh the mprior data + residual theta + omega.
#' @return list(env, built, control, neta, data0, innerTheta0)
#' @noRd
.fsaemInnerSetupCov <- function(ui, data, mpriorMat, control) {
  .built <- .fsaemInnerMpriorUi(ui)
  .neta <- nrow(ui$muRefDataFrame)
  .data <- .fsaemSetMpriorData(data, mpriorMat, .built)
  .env <- .fsaemInnerSetup(.built$ui, .data, matrix(0, nrow(mpriorMat), .neta), control)
  # inner THETA vector template (residual + any time-varying betas), in the inner
  # ui's ntheta order; refreshed each iteration.
  .innerTheta <- .built$ui$iniDf
  .innerTheta <- .innerTheta[!is.na(.innerTheta$ntheta), , drop = FALSE]
  .innerTheta <- .innerTheta[order(.innerTheta$ntheta), ]
  # non-residual inner thetas are time-varying covariate betas; map each to its
  # position in SAEM's Plambda (saemParamsToEstimate order) so it can be
  # refreshed from the live estimate each iteration.
  .betaW <- which(is.na(.innerTheta$err))
  .betaPlambda <- if (length(.betaW)) match(.innerTheta$name[.betaW], ui$saemParamsToEstimate) else integer(0)
  list(env = .env, built = .built, control = control, neta = .neta,
       data0 = data, innerTheta = .innerTheta, betaW = .betaW, betaPlambda = .betaPlambda)
}

#' Refresh the covariate-aware inner at a new estimate.  Rewrites the per-subject
#' mprior data + residual theta + omega and re-solves; the inner model is reused
#' (cache-safe, no recompile).
#' @noRd
.fsaemInnerUpdateCov <- function(setup, mpriorMat, ares, bres, plambda, omega) {
  .data <- .fsaemSetMpriorData(setup$data0, mpriorMat, setup$built)
  # inner theta: residuals from ares/bres; time-varying covariate betas from the
  # live Plambda (falling back to the ini value until Plambda is populated);
  # structural intercepts are the per-subject mprior data.
  .it <- setup$innerTheta
  .theta <- as.numeric(.it$est)
  .residW <- which(!is.na(.it$err))
  if (length(.residW)) {
    .isAdd <- .it$err[.residW] == "add"
    .theta[.residW] <- ifelse(.isAdd, ares[1], bres[1])
  }
  for (i in seq_along(setup$betaW)) {
    .pp <- setup$betaPlambda[i]
    if (!is.na(.pp) && .pp <= length(plambda) && is.finite(plambda[.pp]) && plambda[.pp] != 0) {
      .theta[setup$betaW[i]] <- plambda[.pp]
    }
  }
  setup$env <- .fsaemInnerSetup(setup$built$ui, .data, matrix(0, nrow(mpriorMat), setup$neta),
                                setup$control)
  .fsaemInnerUpdate(setup$env, .theta, omega, matrix(0, nrow(mpriorMat), setup$neta))
  invisible(setup)
}
