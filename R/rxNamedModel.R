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
  ## `role` is kept in the signature for the call sites but no longer shapes the
  ## artifact.  It used to, together with a private build directory: rxode2 keyed the
  ## compiled .c/.so on the PARSED MODEL alone while the emitted C also depends on the
  ## event-sensitivity code generated afterwards, so two builds of one text whose
  ## eventSensCode differed shared one .so and the later build won for both -- silently,
  ## as an augmented model declaring 29 lhs whose calc_lhs computed 4.  nlmixr2est
  ## worked around it by naming every artifact `role_eventSens_md5` and building into
  ## tempdir()/nlmixr2estSens instead of the rxode2 cache.
  ##
  ## rxode2 now folds eventSensCode into its own cache key (nlmixr2/rxode2#1171), so
  ## the variants separate upstream and the models belong back in the rxode2 cache
  ## directory under rxode2's naming.  Requires rxode2 >= 5.1.6 (DESCRIPTION).
  ##
  ## The srcref strip that used to live here moved to rxode2 (rxStripModelSrc), which
  ## now applies it while building the model -- it belongs with the sweep rxode2
  ## already runs over the model env and cmpMgr, and every consumer needs it, not
  ## just this one.  Requires rxode2 >= 5.1.6 (DESCRIPTION).
  rxode2::rxode2(model, ...)
}
