#'@rdname nlmixr2Est
#'@export
nlmixr2Est.ilaplace <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'ilaplace'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="ilaplaceControl")
  .ilaplaceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$ilaplaceControl <- .control
  env$est <- "ilaplace"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="ilaplace")
}
attr(nlmixr2Est.ilaplace, "iov") <- TRUE
attr(nlmixr2Est.ilaplace, "covPresent") <- TRUE
attr(nlmixr2Est.ilaplace, "unbounded") <- .foUnbounded
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.ilaplace, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.ilaplaceControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="ilaplaceControl")
}
