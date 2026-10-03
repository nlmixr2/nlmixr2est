#' Control options for the focep (FOCE+) estimation method
#'
#' This is the first order conditional estimation without eta/residual
#' interaction, but keeping the live conditional residual variance R (the
#' `foce = "foce+"` option of [foceiControl()]).  It is the `foce` method with
#' `foce = "foce+"` forced; use `foce` (est = "foce") for the NONMEM-matching
#' frozen-R behavior.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `focei` instead
#' @param foce FOCE residual-variance mode; for `focepControl()` this is
#'   always `"foce+"` and cannot be changed -- use `foceControl()` for
#'   `"nonmem"`
#' @return focepControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' focepControl()
focepControl <- function(sigdig = 3, ..., interaction = FALSE, foce = "foce+") {
  .control <- foceiControl(sigdig = sigdig, ..., interaction = FALSE, foce = "foce+")
  class(.control) <- "focepControl"
  .control
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.focepControl <- function(control, env) assign("focepControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.focep <- function(control) {
  .getValidCtl(
    control,
    "focepControl",
    convert = c("foceiControl", "foControl", "foiControl"),
    convertFun = .foceiFamilyControlAs
  )
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.focep <- function(x, ...) .nmObjGetControlByClass(x, "focepControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.focep <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "focepControl", assign = FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.focep <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'focep'", .var.name = .ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type = "focepControl")
  .foceiFamilyControlToFoceiControl(env, "focepControl")
  on.exit({
    if (exists("control", envir = .ui)) {
      rm("control", envir = .ui)
    }
  })
  env$focepControl <- .control
  env$est <- "focep"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est = "focep")
}
attr(nlmixr2Est.focep, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.focep, "iov") <- TRUE
attr(nlmixr2Est.focep, "covPresent") <- TRUE
attr(nlmixr2Est.focep, "unbounded") <- .foUnbounded
