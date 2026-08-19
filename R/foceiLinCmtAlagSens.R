# A modeled alag()/f() on a linCmt() compartment moves/scales the DOSE that
# feeds the closed-form solution -- it never appears among linCmtB()'s
# algebraic arguments (p1/v1/ka/...), so ordinary symengine differentiation of
# rx_pred_ (and anything built from it, e.g. rx_r_) is structurally blind to
# it: d(rx_pred_)/d(eta_lag) comes out 0 even though the true gradient is not.
# rxode2/rxode2#1235 added linCmtB(which1=-3) precisely for this: the
# derivative of the linear-system solution wrt a delay applied to every dose.
# This file adds the missing chain-rule term back into the FOCEi eta gradient
# (rx__sens_rx_pred__BY_ETA_*___ / rx__sens_rx_r__BY_ETA_*___), issue #920.
#
# Bioavailability needs no new linCmtB() mode: the system is linear in dose
# amount, so d(pred)/dF = pred/F exactly; only the dF/d(eta) chain-rule factor
# is new.
#
# Scope: the ETA-driven inner-model gradient only (rxUiGet.foceiHdEta /
# rxUiGet.foceiEnv, used by FOCEI/FOCE/EBE). A THETA that drives alag/F without
# going through an ETA needs no correction here (it does not enter the eta
# gradient). The combined eta+theta build (`combSens`, #958) and the fast=TRUE
# analytic 2nd-order Hessian (`.foceiAddHdEta2`) are NOT covered -- known gap,
# left as a follow-up; `fast=TRUE`'s analytic AUGMENTED outer gradient already
# declines linCmt() models entirely (`foceiGradAnalytic.R`), so it is unaffected.
#
# rxode2/rxode2#1236 (infusions): CONFIRMED to reproduce here, not just in the
# plain rxSolve() path the upstream issue was filed against -- a subject whose
# regimen infuses into the lagged/F'd linCmt() compartment reads `NA` for
# `linCmtB(which1=-3)` even when evaluated through this package's own
# generated model (`ind_solve()`-driven `calc_lhs`, same as FOCEi uses).
# Whether a subject's regimen contains an infusion is event-table data this
# build-time step never sees (the compiled model is cached/reused across
# datasets), so it cannot be gated per-subject here. Gated instead on
# `foceiControl(eventSens=)` -- already the documented opt-out for the
# analogous ODE jump-sensitivity feature -- so `eventSens="fd"` is the escape
# hatch for a model that infuses into a lagged/F'd linCmt() compartment.
#
# rxode2/rxode2#1237 (mixed-route dosing) -- CONFIRMED, and BROADER than "two
# explicitly different alag()/f() declarations": per rxode2's own doc comment
# on `linCmtB` (src/linCmt.cpp), which1=-3 "requires that EVERY dose reaching
# the linear system carr[y] the same alag()" -- an oral depot with a modeled
# alag() PLUS an unlagged IV bolus/infusion into central in the SAME regimen
# is just as out of scope as two differently-lagged compartments (an unlagged
# dose is a lag of 0, still a different delay). Confirmed by direct
# comparison to central differences: a mixed depot+central regimen returns a
# gradient biased by a roughly constant, nonzero amount, not NA -- worse than
# #1236 in that it fails silently rather than loudly. `.rxFoceiLinCmtEventRows()`
# only sees the MODEL's alag()/f() declarations, not the DATA's dosing
# routes, so it cannot detect this case; it is the same class of build-time
# blind spot as #1236 and shares its `eventSens="fd"` escape hatch. This
# matters most for f() (bioavailability), which is routinely estimated from
# exactly this kind of paired IV+oral study design.

#' Discover (linCmt compartment, driving ETA/THETA, symengine expr) for a
#' modeled alag()/f() on a linCmt() compartment
#'
#' @param s FOCEi symengine env (post `.sensEtaOrTheta`/`rxUiGet.foceiEtaS`)
#' @param lin Physical linCmt compartment names (`rxode2::.rxLinCmt()`,
#'   sensitivity-compartment names already stripped)
#' @param kind `"lag"` or `"f"` (the symengine symbol prefix; both `alag()` and
#'   `lag()` land on `rx_lag_<cmt>_`, both `F()` and `f()` on `rx_f_<cmt>_`)
#' @param etaVars ETA_i_/THETA_j_ names in gradient-column order
#' @return named list (by compartment) of `list(sym=, drivers=)`; empty list
#'   when the model has no such event-timing dependency
#' @noRd
.rxFoceiLinCmtEventRows <- function(s, lin, kind, etaVars) {
  .rows <- list()
  for (.st in lin) {
    .nm <- paste0("rx_", kind, "_", .st, "_")
    if (!exists(.nm, envir = s, inherits = FALSE)) next
    .sym <- get(.nm, envir = s, inherits = FALSE)
    .free <- tryCatch(vapply(symengine::free_symbols(.sym), as.character, character(1)),
                       error = function(e) character(0))
    .drv <- intersect(.free, etaVars)
    if (length(.drv) == 0L) next
    .rows[[.st]] <- list(sym = .sym, drivers = .drv)
  }
  .rows
}

#' Build the `linCmtB(which1=-3, ...)` moving-boundary sensitivity call from
#' the model's own structural `rx_pred_` call
#'
#' Reuses every argument of the structural linCmtB() call (pointer, time,
#' linCmt/ncmt/oral0, trans, p1..ka) so the shape always matches what
#' `rx_pred_` already solves; only `which1`/`which2` change.  `which2`
#' mirrors the structural encoding: `-1` (reported concentration) maps to the
#' new mode's `-3`; an amount-in-compartment `which1` (structural `-2` case)
#' maps to that same compartment index in the new mode's `which2` slot
#' (rxode2/rxode2#1235's `which2 >= 0` case).
#'
#' @param predArgs `symengine::get_args()` of the structural `rx_pred_` call
#' @return symengine `linCmtB()` FunctionSymbol call
#' @noRd
.rxFoceiLinCmtWhich3Call <- function(predArgs) {
  .args <- predArgs
  .origWhich1 <- .args[[6]]
  .origWhich2 <- .args[[7]]
  .isConc <- isTRUE(all.equal(as.numeric(.origWhich2), -1))
  .args[[6]] <- symengine::S("-3.0")
  .args[[7]] <- if (.isConc) symengine::S("-3.0") else .origWhich1
  do.call(symengine::Function("linCmtB"), as.list(.args))
}

#' Extra d(rx_pred_)/d(eta) terms from a modeled alag()/f() on a linCmt()
#' compartment (the moving-boundary contribution ordinary symengine
#' differentiation of `rx_pred_` cannot see)
#'
#' Two upstream preconditions gate this: `foceiControl(eventSens="fd")` skips
#' the correction entirely (rxode2/rxode2#1236 -- an infused dose into the
#' lagged/F'd compartment reads `NA`; this build-time step has no per-subject
#' event data to gate on more precisely), and more than one lagged linCmt()
#' compartment drops the alag correction, keeping only F() if present
#' (rxode2/rxode2#1237 -- `which1=-3` is one delay shared by every dose, not
#' a per-compartment one).
#'
#' @param x list(rxUi)
#' @param s symengine env with `rx_pred_` and (when present) `rx_lag_<cmt>_`/
#'   `rx_f_<cmt>_` loaded
#' @param etaVars ETA_i_/THETA_j_ names in gradient-column order
#' @return named list, param -> symengine addition to `d(rx_pred_)/d(param)`;
#'   `NULL` when nothing needs correcting
#' @noRd
.rxFoceiLinCmtEventPredExtra <- function(x, s, etaVars) {
  .ui <- x[[1]]
  if (rxode2::.rxLinNcmt(.ui)["numLin"] <= 0L) return(NULL)
  # rxode2/rxode2#1236: which1=-3 reads NA for an infused dose into the
  # lagged/F'd compartment; this build-time step has no per-subject event
  # data to gate on, so reuse the existing analytic-vs-fd opt-out instead.
  if (identical(rxode2::rxGetControl(.ui, "eventSens", "jump"), "fd")) return(NULL)
  if (!exists("rx_pred_", envir = s, inherits = FALSE)) return(NULL)
  .pred <- get("rx_pred_", envir = s, inherits = FALSE)
  if (!identical(tryCatch(symengine::get_name(.pred), error = function(e) ""), "linCmtB")) {
    return(NULL)
  }
  .predArgs <- symengine::get_args(.pred)
  if (length(.predArgs) != 15L) return(NULL)
  .lin <- rxode2::.rxLinCmt(.ui)
  .lin <- grep("^rx__sens_", .lin, value = TRUE, invert = TRUE)
  if (length(.lin) == 0L) return(NULL)
  .lagRows <- .rxFoceiLinCmtEventRows(s, .lin, "lag", etaVars)
  .fRows <- .rxFoceiLinCmtEventRows(s, .lin, "f", etaVars)
  # rxode2/rxode2#1237: which1=-3 is the derivative wrt ONE delay shared by
  # EVERY dose reaching the linear system, not a per-compartment one -- more
  # than one MODELED alag() is only the detectable half of that; a regimen
  # that also doses an unlagged compartment (a delay of 0) is the same
  # violation and is NOT detectable here (that is data, not model,
  # structure -- see the file banner). Skip when the model-visible half is
  # present; leave the existing (structural-only, still-incomplete) gradient
  # unchanged rather than emit a call that answers a different question.
  if (length(.lagRows) > 1L) .lagRows <- list()
  # d(pred)/dF = pred/F assumes ALL of `pred` comes from the F-scaled
  # compartment's dose; more than one MODELED f() is the detectable half of
  # that (same data-invisible caveat as alag() above -- an unscaled dose
  # elsewhere in the regimen breaks the identity too).
  if (length(.fRows) > 1L) .fRows <- list()
  if (length(.lagRows) == 0L && length(.fRows) == 0L) return(NULL)
  .extra <- stats::setNames(vector("list", length(etaVars)), etaVars)
  .accum <- function(param, term) {
    .extra[[param]] <<- if (is.null(.extra[[param]])) term else .extra[[param]] + term
  }
  if (length(.lagRows) == 1L) {
    .lagCall <- .rxFoceiLinCmtWhich3Call(.predArgs)
    .row <- .lagRows[[1]]
    for (.p in .row$drivers) {
      .d <- symengine::D(.row$sym, symengine::S(.p))
      if (paste(.d) %in% c("0", "0.0")) next
      .accum(.p, .d * .lagCall)
    }
  }
  if (length(.fRows) == 1L) {
    .row <- .fRows[[1]]
    for (.p in .row$drivers) {
      .d <- symengine::D(.row$sym, symengine::S(.p))
      if (paste(.d) %in% c("0", "0.0")) next
      .accum(.p, .d * (.pred / .row$sym))
    }
  }
  .extra <- .extra[!vapply(.extra, is.null, logical(1))]
  if (length(.extra) == 0L) return(NULL)
  .extra
}

#' Chain the linCmt() alag/f() moving-boundary correction through a
#' downstream quantity that embeds `rx_pred_`'s linCmtB() call (e.g. `rx_r_`)
#'
#' `rx_r_` is built by substituting `rx_pred_`'s full expression into its own
#' formula (e.g. `(prop.sd*rx_pred_)^2`), so the embedded linCmtB() call is a
#' literal subexpression, not a reference to the `rx_pred_` symbol --
#' `D(rx_r_, rx_pred_)` is identically 0 for a symengine model built this way.
#' `subs()` finds the structurally-identical subexpression by value, so
#' substituting it for a placeholder symbol lets `D()` recover
#' d(target)/d(pred); substituting the placeholder back gives the chain-rule
#' factor entirely in terms of the model's own symbols.  A no-op (identity)
#' when `target` IS `pred` (`rx_pred_` itself, i.e. `rxUiGet.foceiHdEta`).
#'
#' @param target symengine expression that may embed `pred`'s linCmtB() call
#' @param pred `rx_pred_` symengine expression (the embedded call to find)
#' @param extraPred named list from `.rxFoceiLinCmtEventPredExtra()`
#' @return named list, param -> d(target)/d(param) addition (nonzero entries
#'   only); `NULL` when `target` does not depend on `pred` at all
#' @noRd
.rxFoceiLinCmtEventChain <- function(target, pred, extraPred) {
  if (is.null(extraPred) || length(extraPred) == 0L) return(NULL)
  .ph <- symengine::S("rx__linCmtEventPh__")
  .sub <- symengine::subs(target, pred, .ph)
  .dTdPred <- symengine::D(.sub, .ph)
  if (paste(.dTdPred) %in% c("0", "0.0")) return(NULL)
  .dTdPred <- symengine::subs(.dTdPred, .ph, pred)
  lapply(extraPred, function(.e) .dTdPred * .e)
}
