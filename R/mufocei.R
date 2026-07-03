#'@rdname nlmixr2Est
#'@export
nlmixr2Est.mufocei <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  .foceiFamilyControl(env, ..., type="mufoceiControl")
  .mufoceiControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$mufoceiControl <- .control
  env$est <- "mufocei"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="mufocei")
}
attr(nlmixr2Est.mufocei, "covPresent") <- TRUE
attr(nlmixr2Est.mufocei, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.mufocei, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
