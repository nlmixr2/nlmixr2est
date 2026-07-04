#'@rdname nlmixr2Est
#'@export
nlmixr2Est.irlsfocei <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  .foceiFamilyControl(env, ..., type="irlsfoceiControl")
  .irlsfoceiControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$irlsfoceiControl <- .control
  env$est <- "irlsfocei"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="irlsfocei")
}
attr(nlmixr2Est.irlsfocei, "covPresent") <- TRUE
attr(nlmixr2Est.irlsfocei, "unbounded") <- .foUnbounded
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.irlsfocei, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
