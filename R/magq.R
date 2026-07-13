#'@rdname nlmixr2Est
#'@export
nlmixr2Est.magq <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'magq'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="magqControl")
  .magqControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$magqControl <- .control
  env$est <- "magq"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="magq")
}
attr(nlmixr2Est.magq, "iov") <- TRUE
attr(nlmixr2Est.magq, "covPresent") <- TRUE
attr(nlmixr2Est.magq, "unbounded") <- .foUnbounded
# Activates the mu2/mu3/mu4 covariate-rewriting hook (.uiApplyMu2hook,
# R/mu2.R) for this family only, gated on muModel/muRefCovAlg.
attr(nlmixr2Est.magq, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.magqControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="magqControl")
}
