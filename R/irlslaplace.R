#'@rdname nlmixr2Est
#'@export
nlmixr2Est.irlslaplace <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'irlslaplace'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="irlslaplaceControl")
  .irlslaplaceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$irlslaplaceControl <- .control
  env$est <- "irlslaplace"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="irlslaplace")
}
attr(nlmixr2Est.irlslaplace, "iov") <- TRUE
attr(nlmixr2Est.irlslaplace, "covPresent") <- TRUE
attr(nlmixr2Est.irlslaplace, "unbounded") <- .foUnbounded
# Activates the existing mu2/mu3/mu4 algebraic covariate rewriting hook
# (.uiApplyMu2hook, see R/mu2.R) for complex covariate expressions -- see
# .isMuMethod(). Gated on muModel/muRefCovAlg so it only fires for this
# family, matching the "bit-identical when muModel='none'" requirement.
attr(nlmixr2Est.irlslaplace, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.irlslaplaceControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="irlslaplaceControl")
}
