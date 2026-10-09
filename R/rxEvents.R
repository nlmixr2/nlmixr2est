## Emitting events on the rxode2 event bus (see rxode2::rxEventListen).
##
## Everything goes through these wrappers so nlmixr2est keeps working with an
## rxode2 that has no event bus: then they do nothing and scopes just evaluate
## their code.  Only the outermost operation emits, so the fits, solves and
## covariance steps run inside nlmixr2() are silent.

#' Does the installed rxode2 have the event bus?
#' @noRd
.nlmixr2EventBus <- function() {
  exists("rxEventEmit", envir = asNamespace("rxode2"), inherits = FALSE)
}

#' Call an rxode2 event-bus function by name (only when the bus exists)
#' @noRd
.nlmixr2EventFun <- function(name) {
  getExportedValue("rxode2", name)
}

#' Enter the shared operation scope
#' @noRd
.nlmixr2EventEnter <- function() {
  if (.nlmixr2EventBus()) {
    .nlmixr2EventFun(".rxEventEnter")()
  }
  invisible()
}

#' Leave the scope, emitting `event` when this was the outermost operation
#' @noRd
.nlmixr2EventExit <- function(event = NULL, ..., fun = NULL) {
  if (.nlmixr2EventBus()) {
    .nlmixr2EventFun(".rxEventExit")(event, ..., fun = fun)
  }
  invisible()
}

#' Evaluate `expr` inside the shared operation scope
#' @noRd
.nlmixr2EventScope <- function(expr) {
  if (!.nlmixr2EventBus()) {
    return(expr)
  }
  .nlmixr2EventEnter()
  on.exit(.nlmixr2EventExit(), add = TRUE)
  force(expr)
}

#' Emit an event (delivered only outside every scope)
#' @noRd
.nlmixr2EventEmit <- function(event, ..., fun = NULL) {
  if (.nlmixr2EventBus()) {
    .nlmixr2EventFun("rxEventEmit")(event, ..., fun = fun)
  }
  invisible()
}

#' The variable name of the object passed to nlmixr2(), or NULL
#'
#' magrittr's placeholder `.` (and `.x`) are not names the user gave.
#' @noRd
.nlmixr2EventObjName <- function(expr) {
  if (is.name(expr)) {
    .n <- as.character(expr)
    if (!.n %in% c(".", ".x")) {
      return(.n)
    }
  }
  NULL
}

#' Leave nlmixr2()'s scope and emit what the fit produced
#'
#' A fit emits `fitComplete`; a simulation or prediction (`est = "rxSolve"`,
#' `"simulate"`, `"predict"`) emits `solveComplete` with the input object, so
#' `simulate(fit)` stays linked to its fit.  Anything else (e.g. an rxUi from
#' `nlmixr2(model)`) emits nothing.
#' @noRd
.nlmixr2EventExitFit <- function(result, object, call, objName, est) {
  if (inherits(result, "nlmixr2FitCore")) {
    .nlmixr2EventExit(
      "fitComplete",
      fit = result,
      object = object,
      call = call,
      objName = objName,
      source = "fit",
      fun = "nlmixr2"
    )
  } else if (inherits(result, "rxSolve")) {
    .nlmixr2EventExit(
      "solveComplete",
      result = result,
      object = object,
      call = call,
      kind = if (is.character(est) && length(est) == 1L) est else "rxSolve",
      fun = "nlmixr2"
    )
  } else {
    .nlmixr2EventExit()
  }
}

#' Leave an updating function's scope and emit fitUpdate
#'
#' Emits nothing when the function returned its input unchanged (e.g.
#' `addCwres()` on a fit that already has CWRES) unless `force = TRUE`
#' (`setOfv()` changes the fit in place and returns the same object).
#'
#' @param result The function's return value.
#' @param original The fit passed in.
#' @param name The variable name of `original`, or NULL.
#' @param what What changed: "cwres", "npde", "table", "ofv".
#' @param inPlace Whether the user's object now holds the result.
#' @noRd
.nlmixr2EventExitUpdate <- function(result, original, name, what, inPlace, force = FALSE) {
  if (
    !inherits(result, "nlmixr2FitCore") ||
      (!force && identical(result, original))
  ) {
    return(.nlmixr2EventExit())
  }
  .nlmixr2EventExit(
    "fitUpdate",
    fit = result,
    original = original,
    name = if (isTRUE(inPlace)) name else NULL,
    what = what,
    inPlace = isTRUE(inPlace),
    fun = NULL
  )
}

#' Leave a simulating function's scope and emit solveComplete for its result
#' @noRd
.nlmixr2EventExitSolve <- function(result, object, call, kind) {
  if (is.null(result)) {
    return(.nlmixr2EventExit())
  }
  .nlmixr2EventExit("solveComplete", result = result, object = object, call = call, kind = kind, fun = kind)
}

## Fields whose access may compute a deferred objective function
.nmObjObjectiveArgs <- c("logLik", "value", "obf", "ofv", "objf", "OBJF", "objective", "AIC", "BIC")

#' Leave `$`'s scope; emit fitUpdate when the objective was just computed
#' @noRd
.nlmixr2EventExitObjective <- function(fit, env, ofv0) {
  .ofv1 <- if (is.environment(env)) get0("objective", envir = env, inherits = FALSE) else NULL
  if (length(ofv0) == 1L && is.na(ofv0) && length(.ofv1) == 1L && !is.na(.ofv1)) {
    .nlmixr2EventExit("fitUpdate", fit = fit, original = fit, name = NULL, what = "ofv", inPlace = TRUE, fun = NULL)
  } else {
    .nlmixr2EventExit()
  }
}
