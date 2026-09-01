#' Compile a generated model.
#'
#' A thin wrapper over `rxode2::rxode2()`, kept because every generated model in the
#' package goes through one place.
#'
#' It used to do much more.  rxode2 named a model's `.c`/`.so` from the parsed model
#' text alone, while the emitted C also depends on the event-sensitivity code generated
#' afterwards -- so two builds of one text whose event-sensitivity code differed landed
#' on a single `.so`, the later build won for both, and because entry points resolve BY
#' NAME (`R_GetCCallable`) a model object bound to the earlier one silently began
#' executing the replacement.  Measured: an augmented FOCEi sensitivity model built at
#' 193856 bytes was replaced by a 167792-byte build under the same name, after which its
#' `calc_lhs` computed 4 of the 29 columns it declares; every analytic outer gradient
#' then came back non-finite and the fit silently finished on finite differences.
#'
#' This package worked around that by naming every artifact `role_eventSens_md5` and
#' building into `tempdir()/nlmixr2estSens` rather than the rxode2 cache.  rxode2 now
#' folds the event-sensitivity code into its own cache key (nlmixr2/rxode2#1171), so the
#' variants separate upstream and the models belong back in the rxode2 cache directory
#' under rxode2's naming.  Requires rxode2 >= 5.1.6.
#'
#' @param model model text (or anything `rxode2::rxModelVars()` accepts)
#' @param role short role tag, e.g. `"rxInner"`.  No longer shapes the compiled
#'   artifact -- kept so the call sites still read as self-documenting.
#' @param ... passed to `rxode2::rxode2()` (e.g. `eventSens`)
#' @return the compiled rxode2 model
#' @noRd
.nlmixr2estRxode2 <- function(model, role, ...) {
  ## `role` no longer shapes the artifact, and the srcref strip that used to live here
  ## moved to rxode2 (rxStripModelSrc), applied while it builds the model.  Both need
  ## rxode2 >= 5.1.6 (DESCRIPTION).
  .ret <- rxode2::rxode2(model, ...)
  ## Record which event-sensitivity mode built this model.  The focei model bundle is
  ## cached as rxNorm() TEXT and rehydrated with rxode2() (.foceiModelCacheInflate),
  ## and eventSens is NOT recoverable from a compiled model -- rxModelVars() does not
  ## carry it -- so without this tag a rehydrated sensitivity model is rebuilt in "fd"
  ## mode and every dosing-parameter (f/alag/rate/dur) sensitivity silently becomes
  ## zero.  Only set when the caller passed one: rxode2() resolves an absent eventSens
  ## through getOption("rxode2.eventSens"), and `missing(eventSens)` changes what it
  ## builds, so "not passed" has to replay as "not passed".
  .es <- list(...)$eventSens
  if (!is.null(.es)) attr(.ret, "nlmixr2estEventSens") <- .es
  .ret
}
