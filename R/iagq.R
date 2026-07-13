#'@rdname nlmixr2Est
#'@export
nlmixr2Est.iagq <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'iagq'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type="iagqControl")
  .iagqControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$iagqControl <- .control
  env$est <- "iagq"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est="iagq")
}
attr(nlmixr2Est.iagq, "iov") <- TRUE
attr(nlmixr2Est.iagq, "covPresent") <- TRUE
attr(nlmixr2Est.iagq, "unbounded") <- .foUnbounded
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.iagq, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.iagqControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="iagqControl")
}
