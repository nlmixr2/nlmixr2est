#' Control options for the posthoc estimation method
#'
#' This option is for simply getting the maximum a-prior (MAP) also
#' called the posthoc estimates
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiConrol()`
#' @param maxOuterIterations ignored, posthoc always sets this to 0.
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`, though you can set it to be `TRUE` as well.
#' @return posthocControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' posthocControl()
posthocControl <- function(sigdig = 3, ..., interaction = FALSE, maxOuterIterations = NULL) {
  .control <- foceiControl(sigdig = sigdig, ..., maxOuterIterations = 0L, interaction = interaction)
  class(.control) <- "posthocControl"
  .control
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.posthocControl <- function(control, env) assign("posthocControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.posthoc <- function(control) .getValidCtl(control, "posthocControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.posthoc <- function(x, ...) .nmObjGetControlByClass(x, "posthocControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.posthoc <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "posthocControl", assign = FALSE, set = list(maxOuterIterations = 0L))
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.posthoc <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'posthoc'", .var.name = .ui$modelName)
  .control <- env$control
  env$posthocControl <- .control
  .foceiFamilyControl(env, ..., type = "posthocControl")
  .foceiFamilyControlToFoceiControl(env, "posthocControl", set = list(maxOuterIterations = 0L))
  on.exit({
    if (exists("control", envir = .ui)) {
      rm("control", envir = .ui)
    }
  })
  rxode2::rxAssignControlValue(.ui, "maxOuterIterations", 0L)
  .ret <- .foceiFamilyReturn(env, .ui, ..., est = "posthoc")
  .ret
}
attr(nlmixr2Est.posthoc, "covPresent") <- TRUE
attr(nlmixr2Est.posthoc, "unbounded") <- FALSE
# Like "output", posthoc does not estimate anything (maxOuterIterations=0); it
# evaluates an already-specified model, so a prior in ini({}) is not something
# it could silently ignore.  See #938.
attr(nlmixr2Est.posthoc, "nlmixr2Priors") <- "all"
