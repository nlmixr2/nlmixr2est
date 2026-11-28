#' Warn when observed values are close to the ODE absolute tolerance
#'
#' The ODE `atol` is in the units of the states, so observations within
#' `1e3 * atol` of zero can carry integration error comparable to the data.
#' Disable with `options(nlmixr2est.atolCheck = FALSE)` or a control's
#' `atolCheck = FALSE` (`foceiControl()`, `saemControl()`).
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the ui (unmodified; this hook only warns)
#' @noRd
NULL

#' The ODE states of a model, excluding `linCmt()` compartments
#'
#' @param ui rxode2 ui
#' @return character vector of ODE state names
#' @noRd
.atolDvOdeStates <- function(ui) {
  if (is.null(ui$mvL)) {
    .state <- tryCatch(ui$mv0$state, error = function(e) character(0))
    return(as.character(.state))
  }
  .linCmtOdeDdtStates(ui$lstExpr)
}

#' The effective ODE `atol` of an estimation control
#'
#' @param control estimation control, or an `rxControl()`
#' @return the largest `atol`, or `NULL` when the control carries none
#' @noRd
.atolDvControlAtol <- function(control) {
  if (!is.list(control)) {
    return(NULL)
  }
  .rx <- if (inherits(control, "rxControl")) control else control$rxControl
  if (!is.list(.rx)) {
    return(NULL)
  }
  .atol <- suppressWarnings(as.numeric(.rx$atol))
  .atol <- .atol[is.finite(.atol)]
  if (length(.atol) == 0L) {
    return(NULL)
  }
  max(.atol)
}

#' Smallest non-zero absolute DV of the uncensored observations
#'
#' @inheritParams nlmixr2
#' @return the smallest value, or `NULL` when there is none
#' @noRd
.atolDvMinObs <- function(data) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    return(NULL)
  }
  .nms <- tolower(names(data))
  .col <- function(n) {
    .w <- which(.nms == n)
    if (length(.w) == 1L) data[[.w]] else NULL
  }
  .dv <- suppressWarnings(as.numeric(.col("dv")))
  if (length(.dv) == 0L) {
    return(NULL)
  }
  .keep <- is.finite(.dv) & .dv != 0
  .evid <- .col("evid")
  if (!is.null(.evid)) {
    .evid <- suppressWarnings(as.numeric(.evid))
    .keep <- .keep & !is.na(.evid) & .evid == 0
  } else {
    .amt <- suppressWarnings(as.numeric(.col("amt")))
    if (length(.amt) == length(.dv)) {
      .keep <- .keep & (is.na(.amt) | .amt == 0)
    }
  }
  .mdv <- suppressWarnings(as.numeric(.col("mdv")))
  if (length(.mdv) == length(.dv)) {
    .keep <- .keep & (is.na(.mdv) | .mdv == 0)
  }
  .cens <- suppressWarnings(as.numeric(.col("cens")))
  if (length(.cens) == length(.dv)) {
    .keep <- .keep & (is.na(.cens) | .cens == 0)
  }
  if (!any(.keep)) {
    return(NULL)
  }
  min(abs(.dv[.keep]))
}

.preProcessAtolDv <- function(ui, est, data, control) {
  if (isTRUE(est %in% .preProcessHooksSimulationEst)) {
    return(list(ui = ui))
  }
  if (isFALSE(getOption("nlmixr2est.atolCheck", TRUE))) {
    return(list(ui = ui))
  }
  if (is.list(control) && isFALSE(control$atolCheck)) {
    return(list(ui = ui))
  }
  if (length(.atolDvOdeStates(ui)) == 0L) {
    return(list(ui = ui))
  }
  .atol <- .atolDvControlAtol(control)
  if (is.null(.atol)) {
    return(list(ui = ui))
  }
  .minDv <- .atolDvMinObs(data)
  if (is.null(.minDv)) {
    return(list(ui = ui))
  }
  if (.minDv < 1e3 * .atol) {
    warning("DV near ODE atol; rescale units or lower atol", call. = FALSE)
  }
  list(ui = ui)
}
preProcessHooksAdd(".preProcessAtolDv", .preProcessAtolDv)
