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
# Activates mu2/mu3/mu4 covariate rewriting (R/mu2.R); gated on
# muModel/muRefCovAlg for bit-identical behavior when muModel="none".
attr(nlmixr2Est.irlsagq, "mu") <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

#' @export
rxUiDeparse.irlsagqControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="irlsagqControl")
}
