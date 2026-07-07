# Convenience "*f" estimation methods: each is exactly its base method with
# foceiControl(fast = TRUE) as the default (the analytic Almquist-2015 outer
# gradient + Eq-48 warm-start).  They are thin delegates -- the control validator
# forces fast=TRUE through the base control constructor (so the derivative-free
# outerOpt downgrade still runs), and the estimator dispatches to the base method
# (which reports the base `est`).  Method attributes mirror the base family.

#' Force `fast = TRUE` through a base control constructor (respecting the
#' derivative-free `outerOpt` downgrade, which runs inside the constructor).
#' @param control the `getValidNlmixrControl` wrapper list
#' @param ctlFun the base control constructor (e.g. `foceiControl`)
#' @return a base control object with `fast = TRUE`
#' @noRd
.foceiFastCtl <- function(control, ctlFun) {
  .ctl <- control[[1]]
  .l <- if (is.null(.ctl)) list() else unclass(.ctl)
  .l$fast <- TRUE
  do.call(ctlFun, .l)
}

# mu-referenced families share this activation predicate (see R/mu2.R).
.foceiFastMuAttr <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

## ---- focei / foce / focep (no mu) ------------------------------------------

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.foceif <- function(control) .foceiFastCtl(control, foceiControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.foceif <- function(env, ...) nlmixr2Est.focei(env, ...)
attr(nlmixr2Est.foceif, "covPresent") <- TRUE
attr(nlmixr2Est.foceif, "unbounded") <- .foUnbounded
attr(nlmixr2Est.foceif, "iov") <- TRUE

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.focef <- function(control) .foceiFastCtl(control, foceControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.focef <- function(env, ...) nlmixr2Est.foce(env, ...)
attr(nlmixr2Est.focef, "covPresent") <- TRUE
attr(nlmixr2Est.focef, "unbounded") <- .foUnbounded
attr(nlmixr2Est.focef, "iov") <- TRUE

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.focepf <- function(control) .foceiFastCtl(control, focepControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.focepf <- function(env, ...) nlmixr2Est.focep(env, ...)
attr(nlmixr2Est.focepf, "covPresent") <- TRUE
attr(nlmixr2Est.focepf, "unbounded") <- .foUnbounded
attr(nlmixr2Est.focepf, "iov") <- TRUE

## ---- mufoce / mufocep / mufocei (mu-referenced) ----------------------------

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mufoceif <- function(control) .foceiFastCtl(control, mufoceiControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mufoceif <- function(env, ...) nlmixr2Est.mufocei(env, ...)
attr(nlmixr2Est.mufoceif, "covPresent") <- TRUE
attr(nlmixr2Est.mufoceif, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mufoceif, "iov") <- TRUE
attr(nlmixr2Est.mufoceif, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mufocef <- function(control) .foceiFastCtl(control, mufoceControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mufocef <- function(env, ...) nlmixr2Est.mufoce(env, ...)
attr(nlmixr2Est.mufocef, "covPresent") <- TRUE
attr(nlmixr2Est.mufocef, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mufocef, "iov") <- TRUE
attr(nlmixr2Est.mufocef, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mufocepf <- function(control) .foceiFastCtl(control, mufocepControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mufocepf <- function(env, ...) nlmixr2Est.mufocep(env, ...)
attr(nlmixr2Est.mufocepf, "covPresent") <- TRUE
attr(nlmixr2Est.mufocepf, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mufocepf, "iov") <- TRUE
attr(nlmixr2Est.mufocepf, "mu") <- .foceiFastMuAttr

## ---- irlsfoce / irlsfocep / irlsfocei (mu-referenced, irls) ----------------

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.irlsfoceif <- function(control) .foceiFastCtl(control, irlsfoceiControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.irlsfoceif <- function(env, ...) nlmixr2Est.irlsfocei(env, ...)
attr(nlmixr2Est.irlsfoceif, "covPresent") <- TRUE
attr(nlmixr2Est.irlsfoceif, "unbounded") <- .foUnbounded
attr(nlmixr2Est.irlsfoceif, "iov") <- TRUE
attr(nlmixr2Est.irlsfoceif, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.irlsfocef <- function(control) .foceiFastCtl(control, irlsfoceControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.irlsfocef <- function(env, ...) nlmixr2Est.irlsfoce(env, ...)
attr(nlmixr2Est.irlsfocef, "covPresent") <- TRUE
attr(nlmixr2Est.irlsfocef, "unbounded") <- .foUnbounded
attr(nlmixr2Est.irlsfocef, "iov") <- TRUE
attr(nlmixr2Est.irlsfocef, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.irlsfocepf <- function(control) .foceiFastCtl(control, irlsfocepControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.irlsfocepf <- function(env, ...) nlmixr2Est.irlsfocep(env, ...)
attr(nlmixr2Est.irlsfocepf, "covPresent") <- TRUE
attr(nlmixr2Est.irlsfocepf, "unbounded") <- .foUnbounded
attr(nlmixr2Est.irlsfocepf, "iov") <- TRUE
attr(nlmixr2Est.irlsfocepf, "mu") <- .foceiFastMuAttr
