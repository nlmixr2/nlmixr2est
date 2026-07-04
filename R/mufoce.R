#'@rdname nlmixr2Est
#'@export
nlmixr2Est.mufoce <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'mufoce'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="mufoceControl")
  .mufoceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$mufoceControl <- .control
  env$est <- "mufoce"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="mufoce")
}
attr(nlmixr2Est.mufoce, "iov") <- TRUE
attr(nlmixr2Est.mufoce, "covPresent") <- TRUE
attr(nlmixr2Est.mufoce, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.mufoce, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
