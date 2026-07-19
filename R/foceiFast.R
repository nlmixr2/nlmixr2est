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
  # a defaulted (not user-chosen) outer optimizer re-defaults under fast=TRUE
  # (nlminb -> lbfgsb3c); an explicit outerOpt is kept as given
  if (isTRUE(.l$outerOptDefault)) {
    .l$outerOpt <- NULL
    .l$outerOptTxt <- NULL
    .l$outerOptFun <- NULL
    .l$outerOptDefault <- NULL
  }
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

## ---- mfoce / mfocep / mfocei (mu-referenced) ----------------------------

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfoceif <- function(control) .foceiFastCtl(control, mfoceiControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mfoceif <- function(env, ...) nlmixr2Est.mfocei(env, ...)
attr(nlmixr2Est.mfoceif, "covPresent") <- TRUE
attr(nlmixr2Est.mfoceif, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mfoceif, "iov") <- TRUE
attr(nlmixr2Est.mfoceif, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfocef <- function(control) .foceiFastCtl(control, mfoceControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mfocef <- function(env, ...) nlmixr2Est.mfoce(env, ...)
attr(nlmixr2Est.mfocef, "covPresent") <- TRUE
attr(nlmixr2Est.mfocef, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mfocef, "iov") <- TRUE
attr(nlmixr2Est.mfocef, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfocepf <- function(control) .foceiFastCtl(control, mfocepControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.mfocepf <- function(env, ...) nlmixr2Est.mfocep(env, ...)
attr(nlmixr2Est.mfocepf, "covPresent") <- TRUE
attr(nlmixr2Est.mfocepf, "unbounded") <- .foUnbounded
attr(nlmixr2Est.mfocepf, "iov") <- TRUE
attr(nlmixr2Est.mfocepf, "mu") <- .foceiFastMuAttr

## ---- ifoce / ifocep / ifocei (mu-referenced, irls) ----------------

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifoceif <- function(control) .foceiFastCtl(control, ifoceiControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.ifoceif <- function(env, ...) nlmixr2Est.ifocei(env, ...)
attr(nlmixr2Est.ifoceif, "covPresent") <- TRUE
attr(nlmixr2Est.ifoceif, "unbounded") <- .foUnbounded
attr(nlmixr2Est.ifoceif, "iov") <- TRUE
attr(nlmixr2Est.ifoceif, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifocef <- function(control) .foceiFastCtl(control, ifoceControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.ifocef <- function(env, ...) nlmixr2Est.ifoce(env, ...)
attr(nlmixr2Est.ifocef, "covPresent") <- TRUE
attr(nlmixr2Est.ifocef, "unbounded") <- .foUnbounded
attr(nlmixr2Est.ifocef, "iov") <- TRUE
attr(nlmixr2Est.ifocef, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifocepf <- function(control) .foceiFastCtl(control, ifocepControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.ifocepf <- function(env, ...) nlmixr2Est.ifocep(env, ...)
attr(nlmixr2Est.ifocepf, "covPresent") <- TRUE
attr(nlmixr2Est.ifocepf, "unbounded") <- .foUnbounded
attr(nlmixr2Est.ifocepf, "iov") <- TRUE
attr(nlmixr2Est.ifocepf, "mu") <- .foceiFastMuAttr

## ---- agq / magq / iagq (adaptive Gaussian quadrature) --------------------

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.agqf <- function(control) .foceiFastCtl(control, agqControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.agqf <- function(env, ...) nlmixr2Est.agq(env, ...)
attr(nlmixr2Est.agqf, "covPresent") <- TRUE
attr(nlmixr2Est.agqf, "unbounded") <- .foUnbounded
attr(nlmixr2Est.agqf, "iov") <- TRUE

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.magqf <- function(control) .foceiFastCtl(control, magqControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.magqf <- function(env, ...) nlmixr2Est.magq(env, ...)
attr(nlmixr2Est.magqf, "covPresent") <- TRUE
attr(nlmixr2Est.magqf, "unbounded") <- .foUnbounded
attr(nlmixr2Est.magqf, "iov") <- TRUE
attr(nlmixr2Est.magqf, "mu") <- .foceiFastMuAttr

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.iagqf <- function(control) .foceiFastCtl(control, iagqControl)
#' @rdname nlmixr2Est
#' @export
nlmixr2Est.iagqf <- function(env, ...) nlmixr2Est.iagq(env, ...)
attr(nlmixr2Est.iagqf, "covPresent") <- TRUE
attr(nlmixr2Est.iagqf, "unbounded") <- .foUnbounded
attr(nlmixr2Est.iagqf, "iov") <- TRUE
attr(nlmixr2Est.iagqf, "mu") <- .foceiFastMuAttr
