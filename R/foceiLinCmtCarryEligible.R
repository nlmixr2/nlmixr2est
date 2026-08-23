# Carry-eligibility detection for the linCmt() sensitivity-carry (phase 3b.1).
#
# A time-varying covariate on a linCmt() parameter makes the analytic
# d(pred)/d(eta) wrong today (getAlastAD reconstructs the carried Alast
# sensitivity assuming theta is constant across the subject).  The fix (the
# which1=-5/-6/-7 carry recurrence) is only VALIDATED for formulas where the
# eta enters separably: theta = h(cov)*exp(eta) ("mult", d(f)/d(eta) == f)
# or theta = h(cov) + c*eta ("add", d(f)/d(eta) covariate-free).  This file
# only DETECTS eligible (linCmt parameter slot, eta) pairs; nothing is
# substituted here.  Every ambiguous case must come back NOT eligible: a
# false negative keeps today's (already-shipping) behavior, a false positive
# would silently corrupt a gradient.
#
# Conservative rules beyond the shape check, each biased to FALSE:
# - an eta appearing in more than one linCmtB() argument slot is dropped
#   (whole-expression substitution in 3b.3 would have to cover every channel
#   at once; only the single-slot case is validated)
# - a slot containing more than one eta is dropped
# - IOV/occasion etas (iniDf condition != "id") are dropped
# - a slot whose formula references time (t) directly is dropped
# - rx_pred_ must be a bare 15-argument linCmtB() call (matching #920's
#   precondition in foceiLinCmtAlagSens.R); anything else returns no pairs
#
# Data-confirmed variation (varying NA/TRUE/FALSE) and the "linear"
# interpolation error live with the data helpers in foceiLinCmtCarryData.R;
# a CONFIRMED varying covariate under linear interpolation is a clear error,
# not a silent discretization.

.rxFoceiLinCmtCarrySlotNames <- c("p1", "v1", "p2", "p3", "p4", "p5", "ka")

#' Free symbol names of a symengine expression (empty on error/constant)
#' @noRd
.rxFoceiCarryFreeSyms <- function(expr) {
  tryCatch(vapply(symengine::free_symbols(expr), as.character, character(1)),
           error = function(e) character(0))
}

#' Is a symengine expression identically zero?
#' @noRd
.rxFoceiCarryIsZero <- function(expr) {
  paste(expr) %in% c("0", "0.0")
}

#' Detect carry-eligible (linCmt parameter slot, eta) pairs
#'
#' @param x list(rxUi) (matching `.rxFoceiLinCmtEventPredExtra`'s convention)
#' @param s FOCEi symengine env (post `rxUiGet.foceiEtaS`) holding `rx_pred_`
#' @param etaVars ETA_i_ names in gradient-column order
#' @param data optional event data; when given, `varying` is confirmed
#'   per-subject instead of reported as candidate (`NA`)
#' @param interpolation covariate interpolation method; "linear" combined
#'   with a confirmed-varying covariate on an eligible slot is an error
#' @return data.frame, one row per eligible pair: `slot` (1-7 into
#'   p1/v1/p2/p3/p4/p5/ka), `slotName`, `eta` (ETA_n_), `etaName` (original),
#'   `covs` (comma-joined), `shape` ("mult"/"add"), `formula`/`dEtaFormula`
#'   (rxode2 text of the substituted slot expression and its eta derivative),
#'   `varying` (NA without data, else TRUE/FALSE); zero rows when nothing is
#'   eligible
#' @noRd
.rxFoceiLinCmtCarryEligible <- function(x, s, etaVars, data = NULL,
                                        interpolation = c("locf", "nocb", "midpoint", "linear"),
                                        render = TRUE) {
  interpolation <- match.arg(interpolation)
  .ui <- x[[1]]
  .empty <- .rxFoceiCarryEmpty() # nolint: object_usage_linter.
  .predArgs <- .rxFoceiCarryPredArgs(.ui, s)
  if (is.null(.predArgs)) return(.empty)
  .allCovs <- .ui$allCovs
  if (length(.allCovs) == 0L) return(.empty)
  .iniDf <- .ui$iniDf
  .etaDf <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, , drop = FALSE]
  # per-slot free symbols (slots 9-15 of the linCmtB call are p1..ka)
  .slotExpr <- lapply(1:7, function(k) .predArgs[[k + 8L]])
  .slotFree <- lapply(.slotExpr, .rxFoceiCarryFreeSyms)
  .ret <- .empty
  for (.e in seq_along(etaVars)) {
    .row <- .rxFoceiCarryEligibleEta(.e, etaVars, .etaDf, .allCovs, .slotExpr,
                                      .slotFree, data, interpolation, render)
    if (!is.null(.row)) .ret <- rbind(.ret, .row)
  }
  .ret
}

#' The 15 `linCmtB()` arguments of `rx_pred_`, or NULL when the model is not
#' a pure linCmt() model with a bare linCmtB() prediction
#' @noRd
.rxFoceiCarryPredArgs <- function(ui, s) {
  if (rxode2::.rxLinNcmt(ui)["numLin"] <= 0L) return(NULL)
  # mixed ODE + linCmt() models evaluate linCmtB() inside the integrator as
  # well as in the lhs pass; the once-per-row carry stepping is only
  # validated for pure linCmt() models
  if (length(.rxode2stateOdeNoOutput(s)) > 0L) return(NULL)
  if (!exists("rx_pred_", envir = s, inherits = FALSE)) return(NULL)
  .pred <- get("rx_pred_", envir = s, inherits = FALSE)
  if (!identical(tryCatch(symengine::get_name(.pred), error = function(e) ""),
                 "linCmtB")) {
    return(NULL)
  }
  .predArgs <- symengine::get_args(.pred)
  if (length(.predArgs) != 15L) return(NULL)
  .predArgs
}

#' Slot index the eta occupies when it qualifies structurally, else NULL
#' @noRd
.rxFoceiCarryEtaSlot <- function(e, etaVars, etaDf, slotFree) {
  .eta <- etaVars[e]
  .inSlot <- which(vapply(slotFree, function(f) .eta %in% f, logical(1)))
  # eta must live in exactly one slot; multi-slot substitution is unvalidated
  if (length(.inSlot) != 1L) return(NULL)
  .free <- slotFree[[.inSlot]]
  # exactly one eta in that slot, no direct time dependence
  if (sum(etaVars %in% .free) != 1L || "t" %in% .free) return(NULL)
  # IOV/occasion etas carry per-occasion values; unvalidated
  .wEta <- which(etaDf$neta1 == e)
  if (length(.wEta) != 1L) return(NULL)
  .cond <- etaDf$condition[.wEta]
  if (!is.na(.cond) && !identical(.cond, "id")) return(NULL)
  .inSlot
}

#' "mult" / "add" when the eta enters separably, else NULL
#' @noRd
.rxFoceiCarryShape <- function(expr, d, allCovs, etaVars) {
  if (.rxFoceiCarryIsZero(symengine::expand(d - expr))) return("mult")
  .dFree <- .rxFoceiCarryFreeSyms(d)
  if (length(intersect(.dFree, c(allCovs, etaVars, "t"))) == 0L) return("add")
  NULL
}

#' One eligibility row for eta `e`, or NULL
#' @noRd
.rxFoceiCarryEligibleEta <- function(e, etaVars, etaDf, allCovs, slotExpr,
                                     slotFree, data, interpolation, render) {
  .k <- .rxFoceiCarryEtaSlot(e, etaVars, etaDf, slotFree)
  if (is.null(.k)) return(NULL)
  .covs <- intersect(slotFree[[.k]], allCovs)
  if (length(.covs) == 0L) return(NULL)
  .eta <- etaVars[e]
  .expr <- slotExpr[[.k]]
  .d <- tryCatch(symengine::D(.expr, symengine::S(.eta)),
                 error = function(err) NULL)
  if (is.null(.d) || .rxFoceiCarryIsZero(.d)) return(NULL)
  .shape <- .rxFoceiCarryShape(.expr, .d, allCovs, etaVars)
  if (is.null(.shape)) return(NULL)
  .varying <- .rxFoceiCarryVarying(.covs, data) # nolint: object_usage_linter.
  if (identical(interpolation, "linear") && isTRUE(.varying)) {
    stop("time-varying covariate '", paste(.covs, collapse = "', '"),
         "' on a linCmt() parameter needs 'locf', 'nocb' or 'midpoint' ",
         "interpolation; 'linear' cannot be represented by the linCmt() solution",
         call. = FALSE)
  }
  # rxFromSE() is substitute()-based: it deparses a non-character
  # argument's EXPRESSION, so it must always be handed the repr STRING
  # (paste(<Basic>)), never a Basic-yielding call.  render=FALSE keeps
  # the raw symengine repr for callers that render later themselves.
  .fTxt <- paste(.expr)
  .dTxt <- paste(.d)
  if (render) {
    .fTxt <- rxode2::rxFromSE(.fTxt)
    .dTxt <- rxode2::rxFromSE(.dTxt)
  }
  data.frame(slot = .k,
             slotName = .rxFoceiLinCmtCarrySlotNames[.k],
             eta = .eta,
             etaName = etaDf$name[which(etaDf$neta1 == e)],
             covs = paste(.covs, collapse = ","),
             shape = .shape,
             formula = .fTxt,
             dEtaFormula = .dTxt,
             varying = .varying,
             stringsAsFactors = FALSE)
}
