#'@rdname nlmixr2Est
#'@export
nlmixr2Est.irlsagq <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'irlsagq'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="irlsagqControl")
  .irlsagqControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$irlsagqControl <- .control
  env$est <- "irlsagq"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="irlsagq")
}
attr(nlmixr2Est.irlsagq, "iov") <- TRUE
attr(nlmixr2Est.irlsagq, "covPresent") <- TRUE
attr(nlmixr2Est.irlsagq, "unbounded") <- .foUnbounded
# Activates the existing mu2/mu3/mu4 algebraic covariate rewriting hook
# (.uiApplyMu2hook, see R/mu2.R) for complex covariate expressions -- see
# .isMuMethod(). Gated on muModel/muRefCovAlg so it only fires for this
# family, matching the "bit-identical when muModel='none'" requirement.
attr(nlmixr2Est.irlsagq, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.irlsagqControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="irlsagqControl")
}
