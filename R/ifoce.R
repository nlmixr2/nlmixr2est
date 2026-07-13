#'@rdname nlmixr2Est
#'@export
nlmixr2Est.ifoce <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'ifoce'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="ifoceControl")
  .ifoceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$ifoceControl <- .control
  env$est <- "ifoce"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="ifoce")
}
attr(nlmixr2Est.ifoce, "iov") <- TRUE
attr(nlmixr2Est.ifoce, "covPresent") <- TRUE
attr(nlmixr2Est.ifoce, "unbounded") <- .foUnbounded
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.ifoce, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
