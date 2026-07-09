# Dedicated theta-sensitivity model for est="impmap".
#
# The importance-sampling EM updates the non-mu structural thetas by a Newton
# step on the analytic gradient of the individual data log-likelihood, which
# needs d(f)/d(theta) -- the derivative of the prediction with respect to each
# theta.  This compiled model provides it as output lhs
# `rx__sens_rx_pred__BY_THETA_j___`, solved alongside the inner/pred models.
#
# Only the NON-MU STRUCTURAL thetas get a sensitivity: the mu-referenced thetas
# are updated by the EM closed form and residual-error (sigma) thetas do not enter
# the prediction (d(f)/d(sigma)=0).  Restricting to that subset also keeps the
# sensitivity-ODE state count bounded so it fits the shared inner solve buffer.
#
# It is built in the inner-problem order (ODEs, sensitivity ODEs, prediction, then
# the d(f)/d(theta) outputs); no error-model transform or endpoint machinery is
# needed.  Works with or without ODE states.
#
# rx_pred_ is retrieved with the `$` accessor via an intermediate variable
# (the theta sensitivity setup shadows base get/:: and the `$` NSE misbehaves
# when nested in a call).

#' 1-based theta indices (in the THETA_j_ / ntheta ordering) of the non-mu
#' structural thetas: not mu-referenced (intercept or covariate) and not a
#' residual-error parameter.  These are the thetas that need d(f)/d(theta).
#' @noRd
.impmapNonMuStructTheta <- function(ui) {
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta), ]
  .th <- .th[order(.th$ntheta), ]
  .muNames <- unique(ui$muRefDataFrame$theta)
  # mu covariate-coefficient thetas (handled by the regression update) -- exclude
  .covNames <- character(0)
  .mrc <- try(ui$muRefCovariateReplaceDataFrame, silent = TRUE)
  if (!inherits(.mrc, "try-error") && is.data.frame(.mrc) && "covariateParameter" %in% names(.mrc)) {
    .covNames <- unique(.mrc$covariateParameter)
  }
  .keep <- is.na(.th$err) & !(.th$name %in% .muNames) & !(.th$name %in% .covNames)
  .th$ntheta[.keep]
}

#' Build the symengine env carrying the impmap theta-sensitivity model
#' (`..thetaSens`), whose output lhs `rx__sens_rx_pred__BY_THETA_j___` (for each
#' non-mu structural theta j) is d(f)/d(theta_j).
#' @export
rxUiGet.impmapThetaSens <- function(x, ...) {
  .ui <- x[[1]]
  .idx <- .impmapNonMuStructTheta(.ui)
  if (length(.idx) == 0L) return(NULL)
  # Load the (linCmt-promoted) model into symengine, then add sensitivity ODEs
  # for the non-mu structural thetas only.
  .s <- rxUiGet.loadPruneSens(x, ...)
  if (!exists("..maxTheta", .s)) return(NULL)
  .stateVars <- .rxode2stateOdeNoOutput(.s)
  .thetaVars <- paste0("THETA_", .idx, "_")
  rxode2::.rxJacobian(.s, c(.stateVars, .thetaVars))
  rxode2::.rxSens(.s, .thetaVars)
  # rx_pred_ before the sensitivity eval shadows base get/::
  .pred <- .s$`rx_pred_`
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.pred))
  # d(f)/d(theta_j) = D(f, THETA_j) + sum_state rx__sens_state_BY_THETA_j * D(f, state)
  .hd <- vapply(.idx, function(j) {
    .terms <- paste0("D(rx_pred_, THETA_", j, "_)")
    if (length(.stateVars) > 0L) {
      .terms <- c(.terms,
                  paste0("rx__sens_", .stateVars, "_BY_THETA_", j, "___*D(rx_pred_, ",
                         .stateVars, ")"))
    }
    .l <- eval(parse(text = paste0("with(.s, ", paste(.terms, collapse = "+"), ")")))
    paste0("rx__sens_rx_pred__BY_THETA_", j, "___=", rxode2::rxFromSE(.l))
  }, character(1))
  .ddt <- .s$..ddt; if (is.null(.ddt)) .ddt <- character(0)
  .sens <- .s$..sens; if (is.null(.sens)) .sens <- character(0)
  .s$..thetaSens <- paste(c(.ddt, .sens, .prd, .hd, ""), collapse = "\n")
  .s$..thetaSensIdx <- .idx
  .s
}
attr(rxUiGet.impmapThetaSens, "rstudio") <- emptyenv()

#' Compile the impmap theta-sensitivity model.
#'
#' @param ui rxode2 ui object
#' @return a compiled rxode2 model outputting rx__sens_rx_pred__BY_THETA_j___ for
#'   each non-mu structural theta j, or NULL if there are none.
#' @noRd
.impmapThetaSensModel <- function(ui) {
  .s <- rxUiGet.impmapThetaSens(list(ui))
  if (is.null(.s)) return(NULL)
  nlmixr2global$toRxParam <-
    paste0(.uiGetThetaEtaParams(ui, TRUE), "\n", ui$foceiCmtPreModel, "\n")
  nlmixr2global$toRxDvidCmt <- .foceiToCmtLinesAndDvid(ui)
  .toRx(.s$..thetaSens, "compiling theta-sensitivity model...")
}
