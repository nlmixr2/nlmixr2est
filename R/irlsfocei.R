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
# Activates the existing mu2/mu3/mu4 algebraic covariate rewriting hook
# (.uiApplyMu2hook, see R/mu2.R) for complex covariate expressions -- see
# .isMuMethod(). Gated on muModel/muRefCovAlg so it only fires for this
# family, matching the "bit-identical when muModel='none'" requirement.
attr(nlmixr2Est.irlsfocei, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
