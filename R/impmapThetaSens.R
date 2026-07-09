# Dedicated sensitivity model for est="impmap".
#
# The importance-sampling EM updates the non-mu thetas (structural and residual-
# error) by a Newton step on the analytic gradient of the individual data
# log-likelihood.  For a Gaussian endpoint with mean f and variance V that needs,
# per theta, both d(f)/d(theta) and d(V)/d(theta).  This compiled model provides
# them as outputs rx__sens_rx_pred__BY_THETA_j___ = d(f)/d(theta_j) and
# rx__sens_rx_r__BY_THETA_j___ = d(V)/d(theta_j), solved alongside the inner/pred
# models.
#
# Only NON-MU thetas get sensitivities: mu-referenced thetas are updated by the EM
# closed form.  Among those, structural thetas (which enter the ODE states / f)
# carry a sensitivity ODE; residual-error (sigma) thetas enter only V algebraically
# and need no state sensitivity, so restricting the sensitivity ODEs to the
# structural thetas also keeps the state count bounded for the shared solve buffer.
#
# d(f)/d(theta_j) and d(V)/d(theta_j) share the chain rule
#   d(g)/d(theta_j) = D(g, THETA_j) + sum_state rx__sens_state_BY_THETA_j * D(g, state)
# where the state-sensitivity term is present only for the structural thetas.
#
# rx_pred_ / rx_r_ are retrieved with the `$` accessor via an intermediate
# variable (the sensitivity setup shadows base get/:: and the `$` NSE misbehaves
# when nested in a call).

#' 1-based theta indices (THETA_j_ / ntheta ordering) that impmap estimates by
#' the sensitivity Newton step: not mu-referenced (intercept or covariate).  This
#' is the union of the structural thetas (\code{$struct}, which get a sensitivity
#' ODE and d(f)/d(theta)) and the residual-error thetas (\code{$sigma}, which get
#' only d(V)/d(theta)).  \code{$all} is the sorted union.
#' @noRd
.impmapEstTheta <- function(ui) {
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta), ]
  .th <- .th[order(.th$ntheta), ]
  .muNames <- unique(ui$muRefDataFrame$theta)
  .covNames <- character(0)
  .mrc <- try(ui$muRefCovariateReplaceDataFrame, silent = TRUE)
  if (!inherits(.mrc, "try-error") && is.data.frame(.mrc) && "covariateParameter" %in% names(.mrc)) {
    .covNames <- unique(.mrc$covariateParameter)
  }
  .isMu <- (.th$name %in% .muNames) | (.th$name %in% .covNames)
  .fixed <- !is.na(.th$fix) & .th$fix
  .struct <- .th$ntheta[!.isMu & !.fixed & is.na(.th$err)]
  .sigma <- .th$ntheta[!.isMu & !.fixed & !is.na(.th$err)]
  list(struct = .struct, sigma = .sigma, all = sort(unique(c(.struct, .sigma))))
}

#' @noRd
.impmapChainRule <- function(s, target, j, stateVars, structIdx) {
  .terms <- paste0("D(", target, ", THETA_", j, "_)")
  if (j %in% structIdx && length(stateVars) > 0L) {
    .terms <- c(.terms,
                paste0("rx__sens_", stateVars, "_BY_THETA_", j, "___*D(", target, ", ",
                       stateVars, ")"))
  }
  .l <- eval(parse(text = paste0("with(s, ", paste(.terms, collapse = "+"), ")")))
  rxode2::rxFromSE(.l)
}

#' Build the symengine env carrying the impmap sensitivity model
#' (\code{..thetaSens}).  For each estimated non-mu theta j it outputs
#' rx__sens_rx_pred__BY_THETA_j___ = d(f)/d(theta_j) and
#' rx__sens_rx_r__BY_THETA_j___ = d(V)/d(theta_j).
#' @export
rxUiGet.impmapThetaSens <- function(x, ...) {
  .ui <- x[[1]]
  .idx <- .impmapEstTheta(.ui)
  if (length(.idx$all) == 0L) return(NULL)
  .s <- rxUiGet.loadPruneSens(x, ...)
  if (!exists("..maxTheta", .s)) return(NULL)
  .stateVars <- .rxode2stateOdeNoOutput(.s)
  # State sensitivities only for the structural thetas.
  .thetaVars <- paste0("THETA_", .idx$struct, "_")
  if (length(.thetaVars) > 0L) {
    rxode2::.rxJacobian(.s, c(.stateVars, .thetaVars))
    rxode2::.rxSens(.s, .thetaVars)
  }
  .pred <- .s$`rx_pred_`
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.pred))
  # d(f)/d(theta_j): chain rule for structural thetas, 0 for sigma thetas.
  .dfOut <- vapply(.idx$all, function(j) {
    if (j %in% .idx$struct) {
      paste0("rx__sens_rx_pred__BY_THETA_", j, "___=",
             .impmapChainRule(.s, "rx_pred_", j, .stateVars, .idx$struct))
    } else {
      paste0("rx__sens_rx_pred__BY_THETA_", j, "___=0")
    }
  }, character(1))
  # d(V)/d(theta_j): chain rule (structural) or direct partial (sigma).
  .dvOut <- vapply(.idx$all, function(j) {
    paste0("rx__sens_rx_r__BY_THETA_", j, "___=",
           .impmapChainRule(.s, "rx_r_", j, .stateVars, .idx$struct))
  }, character(1))
  .ddt <- .s$..ddt; if (is.null(.ddt)) .ddt <- character(0)
  .sens <- .s$..sens; if (is.null(.sens)) .sens <- character(0)
  .s$..thetaSens <- paste(c(.ddt, .sens, .prd, .dfOut, .dvOut, ""), collapse = "\n")
  .s$..thetaSensIdx <- .idx$all
  .s
}
attr(rxUiGet.impmapThetaSens, "rstudio") <- emptyenv()

#' Compile the impmap sensitivity model.
#'
#' @param ui rxode2 ui object
#' @return a compiled rxode2 model outputting rx__sens_rx_pred__BY_THETA_j___ and
#'   rx__sens_rx_r__BY_THETA_j___ for each estimated non-mu theta j, or NULL if
#'   there are none.
#' @noRd
.impmapThetaSensModel <- function(ui) {
  .s <- rxUiGet.impmapThetaSens(list(ui))
  if (is.null(.s)) return(NULL)
  nlmixr2global$toRxParam <-
    paste0(.uiGetThetaEtaParams(ui, TRUE), "\n", ui$foceiCmtPreModel, "\n")
  nlmixr2global$toRxDvidCmt <- .foceiToCmtLinesAndDvid(ui)
  .toRx(.s$..thetaSens, "compiling sensitivity model...")
}
