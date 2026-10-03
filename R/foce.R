#' Control options for the foce estimation method
#'
#' This is the first order option without the interaction between
#' residuals and etas.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `focei` instead
#' @return foceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' foceControl()
foceControl <- function(sigdig = 3, ..., interaction = FALSE) {
  .control <- foceiControl(sigdig = sigdig, ..., interaction = FALSE)
  class(.control) <- "foceControl"
  .control
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.foceControl <- function(control, env) assign("foceControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.foce <- function(control) {
  .getValidCtl(
    control,
    "foceControl",
    convert = c("foceiControl", "foControl", "foiControl"),
    convertFun = .foceiFamilyControlAs
  )
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.foce <- function(x, ...) .nmObjGetControlByClass(x, "foceControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.foce <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "foceControl", assign = FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.foce <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'foce'", .var.name = .ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type = "foceControl")
  .foceiFamilyControlToFoceiControl(env, "foceControl")
  on.exit({
    if (exists("control", envir = .ui)) {
      rm("control", envir = .ui)
    }
  })
  env$foceControl <- .control
  env$est <- "foce"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est = "foce")
}
attr(nlmixr2Est.foce, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.foce, "iov") <- TRUE
attr(nlmixr2Est.foce, "covPresent") <- TRUE
attr(nlmixr2Est.foce, "unbounded") <- .foUnbounded
