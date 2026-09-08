## Theta sensitivities for a declared-distribution fit, built through the etas.
##
## A declared theta reaches the model ONLY through its own random effect, so
##
##   d(state)/d(theta_j) = sum_k d(state)/d(eta_k) * d(eta_k)/d(theta_j)
##
## and the first factor is already in ind->solve from the solve that produced
## the prediction, while the second is pure algebra on Q(phiU(z); args(theta))
## with the latent FIXED.  So the theta derivative is a READ AND A MULTIPLY
## against the solved buffer, not a new integration -- and the sensitivity
## system scales with the number of declared ETAS rather than with the number
## of parameters those declarations carry.
##
## This file used to also carry a peer model that scored the declared FAMILY
## against the sampled etas (a log-density table, an assembler, two compiled
## models and two odeSwap slots).  That answered the wrong question: the
## complete data is (y, z) with z the latent standard normal, log p(z) is
## theta-free, and eta is a deterministic transform -- so the family density
## never enters the M-step at all.  It is gone; `git log` has it if the family
## term is ever wanted as a penalty rather than as the objective.

#' Chain rule for a declared theta, routed through the declared ETAS
#'
#' `.impmapChainRule()` writes
#'
#'   D(target, THETA_j_) + sum_state rx__sens_<state>_BY_THETA_j___ * D(target, state)
#'
#' which needs one state-sensitivity ODE PER THETA.  A declared theta reaches
#' the model only through its own random effect, so the same derivative is
#'
#'   D(target, THETA_j_)
#'     + sum_state [ sum_k rx__sens_<state>_BY_<eta_k>__ * D(eta_k, THETA_j_) ]
#'                 * D(target, state)
#'
#' which needs one per declared ETA instead -- a fixed, small number that does
#' not grow with the declarations' parameter count.  The state sensitivities
#' are already in `ind->solve` from the solve that produced the prediction, so
#' the theta derivative is a read and a multiply rather than a new integration.
#'
#' `D(eta_k, THETA_j_)` is pure algebra on `Q(phiU(z); args(theta))` with the
#' latent z FIXED -- no solve, and it is the only place the declaration's own
#' parameters appear.
#'
#' On an analytic (`linCmt()`) model there are no ODE states and the sum is
#' empty: `D(target, THETA_j_)` alone already routes through linCmtB's own
#' parameter sensitivities (`rx__sens_central_BY_p1` and friends), which are
#' per linCmt PARAMETER, not per theta.  So this changes nothing there and
#' everything on a true ODE model -- which is where the cost was.
#'
#' @param s symengine environment carrying the model
#' @param target lhs name being differentiated, eg `"rx_pred_"`
#' @param j 1-based theta index
#' @param stateVars ODE states (empty for an analytic model)
#' @param etaNames declared random effects, in declaration order
#' @return the derivative as an rxode2 expression
#' @noRd
#' @author Matthew L. Fidler
.etaDistChainRule <- function(s, target, j, stateVars, etaNames) {
  .terms <- paste0("D(", target, ", THETA_", j, "_)")
  if (length(stateVars) > 0L && length(etaNames) > 0L) {
    for (.st in stateVars) {
      ## d(state)/d(theta_j) = sum_k d(state)/d(eta_k) * d(eta_k)/d(theta_j)
      .parts <- paste0("rx__sens_", .st, "_BY_", .etaDistSensVar(etaNames),
                       "__*D(", etaNames, ", THETA_", j, "_)")
      .terms <- c(.terms,
                  paste0("(", paste(.parts, collapse = "+"), ")*D(", target,
                         ", ", .st, ")"))
    }
  }
  .l <- eval(parse(text = paste0("with(s, ", paste(.terms, collapse = "+"), ")")))
  rxode2::rxFromSE(.l)
}

#' The symbol `.rxSens()` names a declared eta's state sensitivity by
#'
#' Kept as one function so the emitter and the reader cannot drift apart.
#'
#' @param etaNames declared random effects
#' @return character vector, same length
#' @noRd
#' @author Matthew L. Fidler
.etaDistSensVar <- function(etaNames) etaNames

#' Theta-sensitivity model for a declared-distribution fit, built through the etas
#'
#' Emits the SAME lhs names `rxUiGet.saemThetaSens()` does --
#' `rx__sens_rx_pred__BY_THETA_j___` -- so `nonMuGradPhi0()` and the
#' odeSlotThetaSens registration consume it unchanged.  Only the CONSTRUCTION
#' differs: the state sensitivities it asks the solve for are per declared ETA
#' (`.rxSens(s, etaNames)`) rather than per THETA, and the theta derivative is
#' formed from those by `.etaDistChainRule()`.
#'
#' Returns NULL -- and the caller then falls back to the ordinary
#' theta-sensitivity model -- when the model declares nothing, or when the
#' declared etas are not differentiable variables in the loaded model.  It is a
#' cheaper construction of the same quantity, never a different one, so
#' declining costs correctness nothing.
#'
#' @param x rxode2 ui, in a list (rxUiGet convention)
#' @return list with `thetaSens` (model text) and `thetaSensIdx`, matching
#'   `$saemThetaSens`'s shape, or NULL
#' @noRd
#' @author Matthew L. Fidler
#' @export
rxUiGet.etaDistThetaSens <- function(x, ...) {
  .ui <- rxode2::rxUiDecompress(x[[1]])
  .st <- .etaDistDeclGet(.ui)
  if (is.null(.st)) {
    .d <- rxode2::rxUiEtaDists(.ui)
    if (nrow(.d) == 0L) return(NULL)
    .st <- .etaDistDeclStash(.ui, .d)
    if (is.null(.st)) return(NULL)
  }
  .etaNames <- as.character(.st$name)
  if (length(.etaNames) == 0L) return(NULL)
  .idx <- .impmapEstTheta(.ui)$all
  if (length(.idx) == 0L) return(NULL)
  .s <- rxUiGet.loadPruneSens(x, ...)
  if (!exists("..maxTheta", .s)) return(NULL)
  .stateVars <- .rxode2stateOdeNoOutput(.s)
  ## State sensitivities per declared ETA, not per theta.  Skipped entirely
  ## when the model has no ODE states -- an analytic linCmt() already carries
  ## its own parameter sensitivities and needs none of this.
  if (length(.stateVars) > 0L) {
    .ok <- tryCatch({
      rxode2::.rxJacobian(.s, c(.stateVars, .etaNames))
      rxode2::.rxSens(.s, .etaNames)
      TRUE
    }, error = function(e) FALSE)
    if (!.ok) return(NULL)
  }
  .pred <- .s$`rx_pred_`
  .rvar <- .s$`rx_r_`
  .yj <- .s$`rx_yj_`; .lambda <- .s$`rx_lambda_`
  .hi <- .s$`rx_hi_`; .low <- .s$`rx_low_`
  .tbs <- c(paste0("rx_yj_~", rxode2::rxFromSE(.yj)),
            paste0("rx_lambda_~", rxode2::rxFromSE(.lambda)),
            paste0("rx_hi_~", rxode2::rxFromSE(.hi)),
            paste0("rx_low_~", rxode2::rxFromSE(.low)))
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.pred))
  .rr <- paste0("rx_r_=", rxode2::rxFromSE(.rvar))
  .dfOut <- character(0)
  for (.j in .idx) {
    .g <- tryCatch(.etaDistChainRule(.s, "rx_pred_", .j, .stateVars, .etaNames),
                   error = function(e) NULL)
    if (is.null(.g)) return(NULL)
    .dfOut <- c(.dfOut, paste0("rx__sens_rx_pred__BY_THETA_", .j, "___=", .g))
  }
  .ddt <- .s$..ddt; if (is.null(.ddt)) .ddt <- character(0)
  .sens <- .s$..sens; if (is.null(.sens)) .sens <- character(0)
  ## Lightweight return only -- never the symengine environment (see the note
  ## on rxUiGet.impmapThetaSens).
  list(thetaSens = paste(c(.ddt, .sens, .tbs, .prd, .rr, .dfOut, ""),
                         collapse = "\n"),
       thetaSensIdx = .idx)
}
attr(rxUiGet.etaDistThetaSens, "rstudio") <- emptyenv()
