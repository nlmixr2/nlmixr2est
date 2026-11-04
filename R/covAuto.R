# covMethod = "auto": analytic or finite-difference covariance, chosen from the model
# size (#1217).  inst/benchmarks/benchmark-cov-method.R measures both; the fitted log
# ratio of the analytic to the "r,s" covariance-step cost is
# .covAutoCoef[1] + .covAutoCoef[2] * log(ODE states) + .covAutoCoef[3] * log(etas) +
# .covAutoCoef[4] * log(estimated thetas, sigmas included).
.covAutoCoef <- c(0, 0, 0, 0)

#' Predicted analytic / "r,s" covariance-step cost of a model
#' @param nState ODE states
#' @param nEta etas
#' @param nTheta estimated thetas (residual parameters included)
#' @return the predicted ratio
#' @noRd
.covAutoRatio <- function(nState, nEta, nTheta) {
  exp(sum(.covAutoCoef * c(1, log(max(nState, 1L)), log(max(nEta, 1L)), log(max(nTheta, 1L)))))
}

#' Resolve covMethod = "auto" for a fit
#'
#' The analytic covariance when it is predicted to cost less than "r,s", which runs
#' otherwise; scope the analytic hook would decline before anything is computed
#' (FO, linCmt(), no etas) also takes "r,s".  The result is a control with
#' `covAuto = TRUE` and a concrete covType.
#' @param ui rxode2 ui
#' @param control foceiControl() of the fit
#' @return the control
#' @noRd
.foceiCovAuto <- function(ui, control) {
  if (!identical(control$covType, "auto")) {
    return(control)
  }
  control$covAuto <- TRUE
  .iniDf <- ui$iniDf
  .nEta <- sum(!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2 & !.iniDf$fix)
  .nTheta <- sum(!is.na(.iniDf$ntheta) & !.iniDf$fix)
  .nState <- length(.rxode2stateOdeNoOutput(ui$mv0))
  .analytic <- !isTRUE(as.logical(control$fo)) &&
    !inherits(control, c("foControl", "foiControl")) &&
    .nEta > 0L &&
    !.foceiUsesLinCmt(ui) &&
    .covAutoRatio(.nState, .nEta, .nTheta) < 1
  if (.analytic && identical(as.integer(control$covMethod), 1L)) {
    control$covType <- "analytic"
    control$covMethod <- 2L
  } else {
    control$covType <- "fd"
  }
  control
}
