# Carry-eligibility detection for the linCmt() sensitivity-carry (phase 3b.1).
#
# A time-varying covariate on a linCmt() parameter makes the analytic
# d(pred)/d(eta) wrong today (getAlastAD reconstructs the carried Alast
# sensitivity assuming theta is constant across the subject).  The fix (the
# which1=-5/-6/-7 carry recurrence) is only VALIDATED for formulas where the
# eta enters separably: theta = h(cov)*exp(eta) ("mult", d(f)/d(eta) == f)
# or theta = h(cov) + c*eta ("add", d(f)/d(eta) covariate-free).  This file
# only DETECTS eligible (linCmt parameter slot, eta) pairs; nothing is
# substituted yet.  Every ambiguous case must come back NOT eligible: a
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
# "Time-varying" is a data property: with no data, a covariate in
# ui$allCovs is only a CANDIDATE (varying = NA); with data, actual
# within-subject variation is confirmed (TRUE/FALSE).  3b.3 must require
# isTRUE(varying) before substituting.  "linear" covariate interpolation
# cannot be represented by linCmt()'s one-sample-per-row evaluation, so a
# CONFIRMED varying covariate under linear interpolation is a clear error,
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

#' Does a covariate vary within any subject in the data?
#'
#' @return `NA` when the data has no usable ID/covariate column, otherwise
#'   `TRUE`/`FALSE`
#' @noRd
.rxFoceiCarryCovVaries <- function(data, cov) {
  .n <- names(data)
  .idCol <- .n[tolower(.n) == "id"]
  if (length(.idCol) != 1L || !(cov %in% .n)) return(NA)
  any(vapply(split(data[[cov]], data[[.idCol]]),
             function(v) {
               v <- v[!is.na(v)]
               length(unique(v)) > 1L
             }, logical(1)))
}

#' Empty carry-eligibility result (fixed column layout)
#' @noRd
.rxFoceiCarryEmpty <- function() {
  data.frame(slot = integer(0), slotName = character(0),
             eta = character(0), etaName = character(0),
             covs = character(0), shape = character(0),
             formula = character(0), dEtaFormula = character(0),
             varying = logical(0), stringsAsFactors = FALSE)
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
                                        interpolation = c("locf", "nocb", "midpoint", "linear")) {
  interpolation <- match.arg(interpolation)
  .ui <- x[[1]]
  .empty <- .rxFoceiCarryEmpty()
  if (rxode2::.rxLinNcmt(.ui)["numLin"] <= 0L) return(.empty)
  if (!exists("rx_pred_", envir = s, inherits = FALSE)) return(.empty)
  .pred <- get("rx_pred_", envir = s, inherits = FALSE)
  if (!identical(tryCatch(symengine::get_name(.pred), error = function(e) ""),
                 "linCmtB")) {
    return(.empty)
  }
  .predArgs <- symengine::get_args(.pred)
  if (length(.predArgs) != 15L) return(.empty)
  .allCovs <- .ui$allCovs
  if (length(.allCovs) == 0L) return(.empty)
  .iniDf <- .ui$iniDf
  .etaDf <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, , drop = FALSE]
  # per-slot free symbols (slots 9-15 of the linCmtB call are p1..ka)
  .slotExpr <- lapply(1:7, function(k) .predArgs[[k + 8L]])
  .slotFree <- lapply(.slotExpr, .rxFoceiCarryFreeSyms)
  .ret <- .empty
  for (.e in seq_along(etaVars)) {
    .eta <- etaVars[.e]
    .inSlot <- which(vapply(.slotFree, function(f) .eta %in% f, logical(1)))
    # eta must live in exactly one slot; multi-slot substitution is unvalidated
    if (length(.inSlot) != 1L) next
    .k <- .inSlot
    .free <- .slotFree[[.k]]
    # exactly one eta in that slot
    if (sum(etaVars %in% .free) != 1L) next
    # no direct time dependence
    if ("t" %in% .free) next
    # IOV/occasion etas carry per-occasion values; unvalidated
    .wEta <- which(.etaDf$neta1 == .e)
    if (length(.wEta) != 1L) next
    .cond <- .etaDf$condition[.wEta]
    if (!is.na(.cond) && !identical(.cond, "id")) next
    .covs <- intersect(.free, .allCovs)
    if (length(.covs) == 0L) next
    .expr <- .slotExpr[[.k]]
    .d <- tryCatch(symengine::D(.expr, symengine::S(.eta)),
                   error = function(err) NULL)
    if (is.null(.d) || .rxFoceiCarryIsZero(.d)) next
    .dFree <- .rxFoceiCarryFreeSyms(.d)
    .shape <- NULL
    if (.rxFoceiCarryIsZero(symengine::expand(.d - .expr))) {
      .shape <- "mult"
    } else if (length(intersect(.dFree, c(.allCovs, etaVars, "t"))) == 0L) {
      .shape <- "add"
    }
    if (is.null(.shape)) next
    .varying <- NA
    if (!is.null(data)) {
      .v <- vapply(.covs, function(cv) .rxFoceiCarryCovVaries(data, cv), logical(1))
      .varying <- if (all(is.na(.v))) NA else isTRUE(any(.v, na.rm = TRUE))
    }
    if (identical(interpolation, "linear") && isTRUE(.varying)) {
      stop("time-varying covariate '", paste(.covs, collapse = "', '"),
           "' on a linCmt() parameter needs 'locf', 'nocb' or 'midpoint' ",
           "interpolation; 'linear' cannot be represented by the linCmt() solution",
           call. = FALSE)
    }
    .ret <- rbind(.ret,
                  data.frame(slot = .k,
                             slotName = .rxFoceiLinCmtCarrySlotNames[.k],
                             eta = .eta,
                             etaName = .etaDf$name[.wEta],
                             covs = paste(.covs, collapse = ","),
                             shape = .shape,
                             formula = rxode2::rxFromSE(.expr),
                             dEtaFormula = rxode2::rxFromSE(.d),
                             varying = .varying,
                             stringsAsFactors = FALSE))
  }
  .ret
}

#' Carry-eligible pairs straight from a model UI (test/entry convenience)
#'
#' @param ui rxode2 UI model
#' @inheritParams .rxFoceiLinCmtCarryEligible
#' @return see `.rxFoceiLinCmtCarryEligible`
#' @noRd
.foceiLinCmtCarryPairs <- function(ui, data = NULL,
                                   interpolation = c("locf", "nocb", "midpoint", "linear")) {
  .ui <- rxode2::assertRxUi(ui)
  .s <- .ui$foceiEtaS
  .etaVars <- paste0("ETA_", seq_len(.s$..maxEta), "_")
  .rxFoceiLinCmtCarryEligible(list(.ui), .s, .etaVars, data = data,
                              interpolation = match.arg(interpolation))
}
