#'@rdname nlmixr2Est
#'@export
nlmixr2Est.mulaplace <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'mulaplace'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="mulaplaceControl")
  .mulaplaceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$mulaplaceControl <- .control
  env$est <- "mulaplace"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="mulaplace")
}
attr(nlmixr2Est.mulaplace, "iov") <- TRUE
attr(nlmixr2Est.mulaplace, "covPresent") <- TRUE
attr(nlmixr2Est.mulaplace, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.mulaplace, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.mulaplaceControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="mulaplaceControl")
}
