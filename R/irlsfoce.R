#'@rdname nlmixr2Est
#'@export
nlmixr2Est.irlsfoce <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'irlsfoce'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="irlsfoceControl")
  .irlsfoceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$irlsfoceControl <- .control
  env$est <- "irlsfoce"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="irlsfoce")
}
attr(nlmixr2Est.irlsfoce, "iov") <- TRUE
attr(nlmixr2Est.irlsfoce, "covPresent") <- TRUE
attr(nlmixr2Est.irlsfoce, "unbounded") <- .foUnbounded
# Activates the existing mu2/mu3/mu4 algebraic covariate rewriting hook
# (.uiApplyMu2hook, see R/mu2.R) for complex covariate expressions -- see
# .isMuMethod(). Gated on muModel/muRefCovAlg so it only fires for this
# family, matching the "bit-identical when muModel='none'" requirement.
attr(nlmixr2Est.irlsfoce, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}
