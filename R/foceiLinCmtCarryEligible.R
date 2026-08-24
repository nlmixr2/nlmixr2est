# Carry-eligibility detection for the linCmt() sensitivity-carry (phase 3b.1).
#
# A time-varying covariate on a linCmt() parameter makes the analytic
# d(pred)/d(eta) wrong today (getAlastAD reconstructs the carried Alast
# sensitivity assuming theta is constant across the subject).  The fix (the
# which1=-5/-6/-7 carry recurrence) is only VALIDATED for formulas where the
# eta enters separably: theta = h(cov)*exp(eta) ("mult", d(f)/d(eta) == f)
# or theta = h(cov) + c*eta ("add", d(f)/d(eta) covariate-free).  This file
# only DETECTS eligible (eta, channels) pairs; nothing is substituted here.
# Every ambiguous case must come back NOT eligible: a false negative keeps
# today's (already-shipping) behavior, a false positive would silently
# corrupt a gradient.
#
# An eta reaches the state through up to three channels: one linCmt()
# parameter slot (p1..ka), a modeled f() and a modeled alag() on the dosed
# compartment (foceiLinCmtCarryEvent.R).  A pair is made when ANY channel
# needs the carry (a covariate on the slot or on F, or a lag on a
# time-varying kernel) and EVERY channel the eta drives is representable;
# otherwise the eta keeps the status quo gradient (#920's row-local terms).
#
# Conservative rules beyond the shape check, each biased to FALSE:
# - an eta appearing in more than one linCmtB() argument slot is dropped
# - a slot containing more than one eta is dropped
# - IOV/occasion etas (iniDf condition != "id") are dropped
# - a slot whose formula references time (t) directly is dropped
# - a covariate-driven alag() is dropped (its boundary terms need the
#   entering dose's own d(lag)/d(eta))
# - rx_pred_ must contain exactly one distinct 15-argument linCmtB() value
#   call: either rx_pred_ IS that call (the add/prop endpoints) or the call
#   sits inside a larger expression (an ll() endpoint, #1004), which is then
#   factored through the symbol rx_lcConc_ so the carry supplies
#   d(rx_lcConc_)/d(eta) and symengine the outer chain rule; two different
#   calls, or none, return no pairs
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

#' Detect carry-eligible (eta, channels) pairs
#'
#' @param x list(rxUi) (matching `.rxFoceiLinCmtEventPredExtra`'s convention)
#' @param s FOCEi symengine env (post `rxUiGet.foceiEtaS`) holding `rx_pred_`
#' @param etaVars ETA_i_ names in gradient-column order
#' @param data optional event data; when given, `varying` is confirmed
#'   per-subject instead of reported as candidate (`NA`)
#' @param interpolation covariate interpolation method; "linear" combined
#'   with a confirmed-varying covariate on an eligible pair is an error
#' @return data.frame, one row per eligible pair: `slot` (1-7 into
#'   p1/v1/p2/p3/p4/p5/ka, `NA` for an event-only pair), `slotName`, `eta`
#'   (ETA_n_), `etaName` (original), `covs` (comma-joined), `shape`
#'   ("mult"/"add"/`NA`), `formula`/`dEtaFormula` (rxode2 text of the slot
#'   expression and its eta derivative, `NA` without a slot), `varying`
#'   (NA without data, else TRUE/FALSE), `fD` (d(ln F)/d(eta) text or `NA`),
#'   `fCov` (F references a covariate), `fCmt`, `lagD` (d(lag)/d(eta) text
#'   or `NA`), `lagCmt`; zero rows when nothing is eligible
#' @noRd
.rxFoceiLinCmtCarryEligible <- function(x, s, etaVars, data = NULL,
                                        interpolation = c("locf", "nocb", "midpoint", "linear"),
                                        render = TRUE) {
  interpolation <- match.arg(interpolation)
  .ui <- x[[1]]
  .empty <- .rxFoceiCarryEmpty() # nolint: object_usage_linter.
  .allCovs <- .ui$allCovs
  if (length(.allCovs) == 0L) return(.empty)
  .predArgs <- .rxFoceiCarryPredArgs(.ui, s) # nolint: object_usage_linter.
  if (is.null(.predArgs)) return(.empty)
  .iniDf <- .ui$iniDf
  .etaDf <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, , drop = FALSE]
  # per-slot free symbols (slots 9-15 of the linCmtB call are p1..ka)
  .slotExpr <- lapply(1:7, function(k) .predArgs[[k + 8L]])
  .slotFree <- lapply(.slotExpr, .rxFoceiCarryFreeSyms)
  .mods <- .rxFoceiCarryEventMods(.ui, s, etaVars) # nolint: object_usage_linter.
  .ret <- .empty
  for (.e in seq_along(etaVars)) {
    .row <- .rxFoceiCarryEligibleEta(.e, etaVars, .etaDf, .allCovs, .slotExpr,
                                      .slotFree, .mods, data, interpolation, render)
    if (!is.null(.row)) .ret <- rbind(.ret, .row)
  }
  .ret
}

# (the structural-call location itself lives in foceiLinCmtCarryPredCall.R)

#' Is eta `e` a plain between-subject eta (not IOV/occasion)?
#' @noRd
.rxFoceiCarryEtaIdOnly <- function(e, etaDf) {
  .wEta <- which(etaDf$neta1 == e)
  if (length(.wEta) != 1L) return(FALSE)
  .cond <- etaDf$condition[.wEta]
  is.na(.cond) || identical(.cond, "id")
}

#' The eta's slot channel: list(k, g, shape, covs) for a single qualifying
#' slot, list() when the eta is in no slot, NULL when it is in a slot the
#' carry cannot represent
#' @noRd
.rxFoceiCarryEtaSlotInfo <- function(eta, etaVars, allCovs, slotExpr, slotFree) {
  .inSlot <- which(vapply(slotFree, function(f) eta %in% f, logical(1)))
  if (length(.inSlot) == 0L) return(list())
  # eta must live in exactly one slot; multi-slot substitution is unvalidated
  if (length(.inSlot) != 1L) return(NULL)
  .free <- slotFree[[.inSlot]]
  # exactly one eta in that slot, no direct time dependence
  if (sum(etaVars %in% .free) != 1L || "t" %in% .free) return(NULL)
  .expr <- slotExpr[[.inSlot]]
  .d <- tryCatch(symengine::D(.expr, symengine::S(eta)), error = function(err) NULL)
  if (is.null(.d) || .rxFoceiCarryIsZero(.d)) return(NULL)
  .shape <- .rxFoceiCarryShape(.expr, .d, allCovs, etaVars)
  if (is.null(.shape)) return(NULL)
  list(k = .inSlot, g = .d, expr = .expr, shape = .shape,
       covs = intersect(.free, allCovs))
}

#' "mult" / "add" when the eta enters separably, else NULL
#' @noRd
.rxFoceiCarryShape <- function(expr, d, allCovs, etaVars) {
  if (.rxFoceiCarryIsZero(symengine::expand(d - expr))) return("mult")
  .dFree <- .rxFoceiCarryFreeSyms(d)
  if (length(intersect(.dFree, c(allCovs, etaVars, "t"))) == 0L) return("add")
  NULL
}

#' Render a symengine repr to model text (or pass it through)
#' @noRd
.rxFoceiCarryTxt <- function(repr, render) {
  if (is.null(repr) || is.na(repr)) return(NA_character_)
  if (render) rxode2::rxFromSE(repr) else repr
}

#' One eligibility row for eta `e`, or NULL
#' @noRd
.rxFoceiCarryEligibleEta <- function(e, etaVars, etaDf, allCovs, slotExpr,
                                     slotFree, mods, data, interpolation, render) {
  .eta <- etaVars[e]
  if (!.rxFoceiCarryEtaIdOnly(e, etaDf)) return(NULL)
  .jump <- .rxFoceiCarryEtaJump(.eta, mods, allCovs) # nolint: object_usage_linter.
  if (!isTRUE(.jump$ok)) return(NULL)
  .slot <- .rxFoceiCarryEtaSlotInfo(.eta, etaVars, allCovs, slotExpr, slotFree)
  if (is.null(.slot)) return(NULL)
  .hasSlot <- length(.slot) > 0L
  .why <- .rxFoceiCarryWhy(if (.hasSlot) .slot else NULL, .jump, mods, slotFree, allCovs)
  if (length(.why) == 0L) return(NULL)
  .varying <- .rxFoceiCarryVarying(.why, data) # nolint: object_usage_linter.
  if (identical(interpolation, "linear") && isTRUE(.varying)) {
    stop("time-varying covariate '", paste(.why, collapse = "', '"),
         "' on a linCmt() parameter needs 'locf', 'nocb' or 'midpoint' ",
         "interpolation; 'linear' cannot be represented by the linCmt() solution",
         call. = FALSE)
  }
  .rxFoceiCarryPairRow(.eta, etaDf$name[which(etaDf$neta1 == e)], # nolint: object_usage_linter.
                       if (.hasSlot) .slot else NULL, .jump, mods, .why, .varying, render)
}

#' The covariates that make this eta's gradient history-dependent: those on
#' its slot, those in dlnF, and -- for a lag channel -- every covariate on
#' the kernel; empty when the status quo gradient is already exact
#' @noRd
.rxFoceiCarryWhy <- function(slot, jump, mods, slotFree, allCovs) {
  .why <- if (is.null(slot)) character(0) else slot$covs
  if (isTRUE(jump$fCov)) {
    .why <- c(.why, intersect(.rxFoceiCarryFreeSyms(mods$f$sym), allCovs))
  }
  if (!is.null(jump$lagD) && .rxFoceiCarryKernelHasCov(slotFree, allCovs)) { # nolint: object_usage_linter.
    .why <- c(.why, unlist(lapply(slotFree, intersect, allCovs)))
  }
  unique(.why)
}
