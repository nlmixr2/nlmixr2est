#' Does a fit's model solve ODEs?
#'
#' A `linCmt()`-only model reports compartments but has no `d/dt()`, so it is
#' not an ODE model here.
#' @param fit nlmixr2 fit
#' @return logical
#' @noRd
.odeSettleIsOde <- function(fit) {
  .mv <- tryCatch(rxode2::rxModelVars(fit$ui), error = function(e) NULL)
  if (is.null(.mv)) {
    return(FALSE)
  }
  .norm <- .mv$model[["normModel"]]
  length(.norm) == 1L && grepl("d/dt(", .norm, fixed = TRUE)
}

#' An `rxControl()` with every ODE tolerance a factor tighter
#'
#' @param rxControl the fit's `rxControl()`
#' @param factor the tightening factor
#' @return the tightened `rxControl()`
#' @noRd
.odeSettleTighten <- function(rxControl, factor = 10) {
  for (.n in c("atol", "rtol", "atolSens", "rtolSens", "ssAtol", "ssRtol", "ssAtolSens", "ssRtolSens")) {
    if (is.numeric(rxControl[[.n]])) {
      rxControl[[.n]] <- rxControl[[.n]] / factor
    }
  }
  rxControl
}

#' The objective at the final estimates under 10x tighter ODE tolerances
#'
#' Re-solves the ETAs (from the fit's) and the objective with no outer step.
#' @param fit nlmixr2 fit
#' @return the objective, `NA` when the re-solve fails
#' @noRd
.odeSettleObjf <- function(fit) {
  .ctl <- fit$foceiControl
  .ctl$rxControl <- .odeSettleTighten(.ctl$rxControl)
  .ctl$etaMat <- if (isTRUE(.ctl$fo)) NULL else fit$etaMat
  .ctl$maxOuterIterations <- 0L
  .ctl$calcTables <- FALSE
  .ctl$covMethod <- 0L
  .ctl$compress <- FALSE
  .ctl$odeSettleCheck <- FALSE
  .ctl$print <- 0L
  .ctl <- do.call(foceiControl, .ctl)
  .newFit <- tryCatch(
    suppressMessages(suppressWarnings(
      .nlmixr2PriorGateBypass(
        nlmixr2(fit, nlme::getData(fit), "focei", control = .ctl)
      )
    )),
    error = function(e) NULL
  )
  if (is.null(.newFit)) {
    return(NA_real_)
  }
  as.double(.newFit$objective)
}

#' The `$runInfo` note for an objective that moves with the ODE tolerance
#'
#' @param shift the objective at the tighter tolerance minus the fit's
#' @param tol the threshold
#' @return the message, or `NULL` when the objective is settled
#' @noRd
.odeSettleMsg <- function(shift, tol) {
  if (!is.finite(shift)) {
    return("OFV not finite at 10x tighter ODE tolerance")
  }
  if (abs(shift) <= tol) {
    return(NULL)
  }
  sprintf("OFV moves %s at 10x tighter ODE tolerance", format(signif(shift, 3)))
}

#' Warn when a fit's objective is not settled with respect to the ODE tolerance
#'
#' Opt-in with `foceiControl(odeSettleCheck = TRUE)`; skipped for models
#' without ODEs and for fits with no outer optimization.  The tighter
#' objective is kept as `odeSettleObjf` on the fit.
#' @param fit the finalized fit
#' @return the objective shift, invisibly (`NULL` when not checked)
#' @noRd
.odeSettleCheck <- function(fit) {
  .ctl <- tryCatch(fit$control, error = function(e) NULL)
  if (!is.list(.ctl) || !isTRUE(.ctl[["odeSettleCheck"]])) {
    return(invisible(NULL))
  }
  if (identical(as.integer(.ctl[["maxOuterIterations"]]), 0L) || !.odeSettleIsOde(fit)) {
    return(invisible(NULL))
  }
  .tol <- .ctl[["odeSettleTol"]]
  if (!is.numeric(.tol) || length(.tol) != 1L) {
    .tol <- 0.1
  }
  # the nested fit resets the timing stack, so nlmixrWithTiming() cannot record this
  .pt <- proc.time()
  .objf <- .odeSettleObjf(fit)
  .env <- if (is.environment(fit)) fit else fit$env
  try(nlmixrAddTiming(.env, "odeSettle", (proc.time() - .pt)[["elapsed"]]), silent = TRUE)
  assign("odeSettleObjf", .objf, envir = .env)
  .shift <- .objf - as.double(fit$objective)
  .msg <- .odeSettleMsg(.shift, .tol)
  if (!is.null(.msg)) {
    warning(.msg, call. = FALSE)
  }
  invisible(.shift)
}
