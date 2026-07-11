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
