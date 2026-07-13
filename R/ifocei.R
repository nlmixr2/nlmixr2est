#'@rdname nlmixr2Est
#'@export
nlmixr2Est.ifocei <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  .foceiFamilyControl(env, ..., type="ifoceiControl")
  .ifoceiControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$ifoceiControl <- .control
  env$est <- "ifocei"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="ifocei")
}
attr(nlmixr2Est.ifocei, "covPresent") <- TRUE
attr(nlmixr2Est.ifocei, "unbounded") <- .foUnbounded
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.ifocei, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
