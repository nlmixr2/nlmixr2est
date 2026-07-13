#'@rdname nlmixr2Est
#'@export
nlmixr2Est.mlaplace <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'mlaplace'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="mlaplaceControl")
  .mlaplaceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$mlaplaceControl <- .control
  env$est <- "mlaplace"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="mlaplace")
}
attr(nlmixr2Est.mlaplace, "iov") <- TRUE
attr(nlmixr2Est.mlaplace, "covPresent") <- TRUE
attr(nlmixr2Est.mlaplace, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.mlaplace, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.mlaplaceControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="mlaplaceControl")
}
