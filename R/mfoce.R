#'@rdname nlmixr2Est
#'@export
nlmixr2Est.mfoce <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'mfoce'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="mfoceControl")
  .mfoceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$mfoceControl <- .control
  env$est <- "mfoce"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="mfoce")
}
attr(nlmixr2Est.mfoce, "iov") <- TRUE
attr(nlmixr2Est.mfoce, "covPresent") <- TRUE
attr(nlmixr2Est.mfoce, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.mfoce, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
