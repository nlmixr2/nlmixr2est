#'@rdname nlmixr2Est
#'@export
nlmixr2Est.muagq <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'muagq'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="muagqControl")
  .muagqControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$muagqControl <- .control
  env$est <- "muagq"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="muagq")
}
attr(nlmixr2Est.muagq, "iov") <- TRUE
attr(nlmixr2Est.muagq, "covPresent") <- TRUE
attr(nlmixr2Est.muagq, "unbounded") <- .foUnbounded
# Activates the existing mu2/mu3/mu4 algebraic covariate rewriting hook
# (.uiApplyMu2hook, see R/mu2.R) for complex covariate expressions -- see
# .isMuMethod(). Gated on muModel/muRefCovAlg so it only fires for this
# family, matching the "bit-identical when muModel='none'" requirement.
attr(nlmixr2Est.muagq, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.muagqControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="muagqControl")
}
