#'@rdname nlmixr2Est
#'@export
nlmixr2Est.mfocei <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  .foceiFamilyControl(env, ..., type="mfoceiControl")
  .mfoceiControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$mfoceiControl <- .control
  env$est <- "mfocei"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="mfocei")
}
attr(nlmixr2Est.mfocei, "covPresent") <- TRUE
attr(nlmixr2Est.mfocei, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.mfocei, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
