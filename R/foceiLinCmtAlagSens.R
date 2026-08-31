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
# The correction is emitted PER COMPARTMENT, using rxode2's per-origin
# decomposition of the linCmt() amounts (rxode2/rxode2#1119 part B):
# `linCmtB(which1=-9)` is the derivative wrt a delay on ONE compartment's
# doses, and `linCmtB(which1=-10)` returns the part of the prediction that
# arrived through one compartment, so d(pred)/dF_q = A^(q)/F_q. That is what
# makes a mixed-route regimen -- a modeled alag()/f() on an oral depot
# alongside an unlagged, unscaled IV dose into central, exactly the paired
# design bioavailability is usually estimated from -- come out right.
# `which1=-3` (one delay shared by every dose) and `pred/F` (all of `pred`
# from the F-scaled dose) are only correct when the whole regimen feeds one
# compartment, and fail SILENTLY, with a roughly constant bias, when it does
# not (rxode2/rxode2#1237); neither is used here any more, and the
# one-lagged-compartment / one-F'd-compartment restrictions they forced are
# gone with them.
#
# Remaining gap (rxode2/rxode2#1236's tail): a subject whose regimen reaches
# the lagged compartment with a STEADY-STATE INFUSION still reads `NA` -- the
# SS infusion's rate is not carried past the SS solve, so -dA/dt is not
# recoverable there. Whether a subject's regimen contains one is event-table
# data this build-time step never sees (the compiled model is cached/reused
# across datasets), so it cannot be gated per-subject here. Gated instead on
# `foceiControl(eventSens=)` -- already the documented opt-out for the
# analogous ODE jump-sensitivity feature -- so `eventSens="fd"` is the escape
# hatch for such a model.

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

# rxode2's which2 packing for the per-origin modes (`q*8 + out`, with `out`
# the compartment whose amount is wanted or 7 for the reported
# concentration); see RX_LINCMT_ORIGIN_W2 / RX_LINCMT_ORIGIN_CONC in
# rxode2's src/linCmt.cpp.
.rxFoceiLinCmtOriginW2 <- 8L
.rxFoceiLinCmtOriginConc <- 7L

#' Build a per-origin `linCmtB()` call from the model's own structural
#' `rx_pred_` call
#'
#' Reuses every argument of the structural linCmtB() call (pointer, time,
#' linCmt/ncmt/oral0, trans, p1..ka) so the shape always matches what
#' `rx_pred_` already solves; only `which1`/`which2` change.  The output the
#' structural call asks for carries over: a reported concentration
#' (structural `which2 = -1`) becomes `RX_LINCMT_ORIGIN_CONC`, an
#' amount-in-compartment (structural `which1`) stays that compartment index.
#'
#' @param predArgs `symengine::get_args()` of the structural `rx_pred_` call
#' @param which1 `-9` (per-compartment dose-time) or `-10` (per-compartment
#'   amounts, the bioavailability chain-rule factor)
#' @param q linCmt() block index of the origin compartment (0-based, depot
#'   first when oral)
#' @return symengine `linCmtB()` FunctionSymbol call
#' @noRd
.rxFoceiLinCmtOriginCall <- function(predArgs, which1, q) {
  .args <- predArgs
  .origWhich1 <- .args[[6]]
  .origWhich2 <- .args[[7]]
  .isConc <- isTRUE(all.equal(as.numeric(.origWhich2), -1))
  .out <- if (.isConc) .rxFoceiLinCmtOriginConc else as.numeric(.origWhich1)
  .args[[6]] <- symengine::S(sprintf("%.1f", which1))
  .args[[7]] <- symengine::S(sprintf("%.1f", q * .rxFoceiLinCmtOriginW2 + .out))
  do.call(symengine::Function("linCmtB"), as.list(.args))
}

#' Extra d(rx_pred_)/d(eta) terms from a modeled alag()/f() on a linCmt()
#' compartment (the moving-boundary contribution ordinary symengine
#' differentiation of `rx_pred_` cannot see)
#'
#' Every lagged and every F'd linCmt() compartment contributes its own term,
#' read off rxode2's per-origin decomposition of the amounts, so a model may
#' lag or scale as many of them as it likes.  One upstream precondition gates
#' this: `foceiControl(eventSens="fd")` skips the correction entirely (a
#' steady-state infusion into the lagged compartment reads `NA`, and this
#' build-time step has no per-subject event data to gate on more precisely).
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
  # A steady-state infusion into the lagged compartment still reads NA; this
  # build-time step has no per-subject event data to gate on, so reuse the
  # existing analytic-vs-fd opt-out instead.
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
  if (length(.lagRows) == 0L && length(.fRows) == 0L) return(NULL)
  # .lin is the linCmt() block order (depot first when oral), which is the
  # 0-based origin index the which2 packing takes.
  .qOf <- stats::setNames(seq_along(.lin) - 1L, .lin)
  .extra <- stats::setNames(vector("list", length(etaVars)), etaVars)
  .accum <- function(param, term) {
    .extra[[param]] <<- if (is.null(.extra[[param]])) term else .extra[[param]] + term
  }
  # d(pred)/d(alag_q) = the dose-time sensitivity of compartment q's own
  # doses; d(pred)/dF_q = A^(q)/F_q, the part of pred that came through q.
  # Both are per compartment, so every lagged/F'd compartment contributes.
  .addRows <- function(rows, which1, factor) {
    for (.cmt in names(rows)) {
      .q <- .qOf[[.cmt]]
      .call <- .rxFoceiLinCmtOriginCall(.predArgs, which1, .q)
      .row <- rows[[.cmt]]
      for (.p in .row$drivers) {
        .d <- symengine::D(.row$sym, symengine::S(.p))
        if (paste(.d) %in% c("0", "0.0")) next
        .accum(.p, .d * factor(.call, .row$sym))
      }
    }
  }
  .addRows(.lagRows, -9, function(call, sym) call)
  .addRows(.fRows, -10, function(call, sym) call / sym)
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
