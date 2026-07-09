# Dedicated theta-sensitivity model for est="impmap".
#
# The importance-sampling EM updates the non-mu structural thetas by a Newton
# step on the analytic gradient of the individual data log-likelihood, which
# needs d(f)/d(theta) -- the derivative of the prediction with respect to each
# theta.  This compiled model provides it as output lhs
# `rx__sens_rx_pred__BY_THETA_j___`, solved alongside the inner/pred models.
#
# It is a minimal model built in the inner-problem order (ODEs -> sensitivity
# ODEs -> prediction -> the d(f)/d(theta) outputs); no error-model transform or
# endpoint machinery is needed since we only want the prediction sensitivities.
# It works with or without ODE states (linCmt/algebraic models have empty ODE and
# sensitivity blocks and analytic prediction sensitivities).
#
# rx_pred_ is converted (rxFromSE + get) BEFORE the d(f)/d(theta) apply, because
# evaluating the sensitivity chain rule shadows base `get`/`::` in that scope.

#' Build the symengine env carrying the impmap theta-sensitivity model
#' (`..thetaSens`), whose output lhs `rx__sens_rx_pred__BY_THETA_j___` is
#' d(f)/d(theta_j).
#' @export
rxUiGet.impmapThetaSens <- function(x, ...) {
  .s <- rxUiGet.foceiThetaS(x, ...)
  # The theta sensitivity setup shadows base `get` in this scope with a symengine
  # variant, so retrieve rx_pred_ with the `$` accessor (not get()).  Assign to a
  # variable first -- the `$` symengine NSE misbehaves when nested in the call.
  .pred <- .s$`rx_pred_`
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.pred))
  .stateVars <- .rxode2stateOdeNoOutput(.s)
  .grd <- rxode2::rxExpandFEta_(.stateVars, .s$..maxTheta, 1L, isTheta = TRUE)
  .hd <- apply(.grd, 1, function(z) {
    .l <- eval(parse(text = z["calc"]))
    paste0(z["dfe"], "=", rxode2::rxFromSE(.l))
  })
  .ddt <- .s$..ddt; if (is.null(.ddt)) .ddt <- character(0)
  .sens <- .s$..sens; if (is.null(.sens)) .sens <- character(0)
  # Inner-problem order: ODEs, sensitivity ODEs, prediction, then d(f)/d(theta).
  .s$..thetaSens <- paste(c(.ddt, .sens, .prd, .hd, ""), collapse = "\n")
  .s
}
attr(rxUiGet.impmapThetaSens, "rstudio") <- emptyenv()

#' Compile the impmap theta-sensitivity model.
#'
#' @param ui rxode2 ui object
#' @return a compiled rxode2 model outputting rx__sens_rx_pred__BY_THETA_j___ per
#'   theta.
#' @noRd
.impmapThetaSensModel <- function(ui) {
  .s <- rxUiGet.impmapThetaSens(list(ui))
  nlmixr2global$toRxParam <-
    paste0(.uiGetThetaEtaParams(ui, TRUE), "\n", ui$foceiCmtPreModel, "\n")
  nlmixr2global$toRxDvidCmt <- .foceiToCmtLinesAndDvid(ui)
  .toRx(.s$..thetaSens, "compiling theta-sensitivity model...")
}
