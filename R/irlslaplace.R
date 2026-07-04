#'@rdname nlmixr2Est
#'@export
nlmixr2Est.irlslaplace <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'irlslaplace'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="irlslaplaceControl")
  .irlslaplaceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$irlslaplaceControl <- .control
  env$est <- "irlslaplace"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="irlslaplace")
}
attr(nlmixr2Est.irlslaplace, "iov") <- TRUE
attr(nlmixr2Est.irlslaplace, "covPresent") <- TRUE
attr(nlmixr2Est.irlslaplace, "unbounded") <- .foUnbounded
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.irlslaplace, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.irlslaplaceControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="irlslaplaceControl")
}
