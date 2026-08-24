# Data-dependent confirmation for the linCmt() sensitivity carry, the
# eligibility row layout, the fit-time data checks the jump channels need,
# and the UI-facing entry point.  "Time-varying" is a data property: with no data a
# covariate in ui$allCovs is only a CANDIDATE (varying = NA); with data,
# actual within-subject variation is confirmed (TRUE/FALSE).  The symbolic
# detection itself lives in foceiLinCmtCarryEligible.R.

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

#' NA without data, else whether any of the covariates varies within a subject
#' @noRd
.rxFoceiCarryVarying <- function(covs, data) {
  if (is.null(data)) return(NA)
  .v <- vapply(covs, function(cv) .rxFoceiCarryCovVaries(data, cv), logical(1))
  if (all(is.na(.v))) NA else isTRUE(any(.v, na.rm = TRUE))
}

#' Empty carry-eligibility result (fixed column layout)
#' @noRd
.rxFoceiCarryEmpty <- function() {
  data.frame(slot = integer(0), slotName = character(0),
             eta = character(0), etaName = character(0),
             covs = character(0), shape = character(0),
             formula = character(0), dEtaFormula = character(0),
             varying = logical(0),
             fD = character(0), fCov = logical(0), fCmt = character(0),
             lagD = character(0), lagCmt = character(0),
             stringsAsFactors = FALSE)
}

#' Session memo for the data-independent candidate pairs, keyed by the
#' focei model digest (which already covers the model text, iniDf,
#' linCmtSensCarry, the rxode2 carry capability and covsInterpolation --
#' everything the symbolic detection depends on).  Only the `data = NULL`
#' candidate result is memoized; every data-dependent conclusion (varying,
#' ss/evid, jump-regimen checks) is re-derived per fit.
#' @noRd
.carryPairsMemo <- new.env(parent = emptyenv())

#' Drop every memoized entry (tests; a session never needs this itself)
#' @noRd
.foceiLinCmtCarryMemoClear <- function() {
  rm(list = ls(envir = .carryPairsMemo, all.names = FALSE), envir = .carryPairsMemo)
  invisible(NULL)
}

#' Memo hit/miss counters (mechanism observable for tests)
#' @noRd
.foceiLinCmtCarryMemoStats <- function(reset = FALSE) {
  .st <- c(hits = get0("hits", envir = .carryPairsMemo, ifnotfound = 0L),
           misses = get0("misses", envir = .carryPairsMemo, ifnotfound = 0L))
  if (reset) {
    assign("hits", 0L, envir = .carryPairsMemo)
    assign("misses", 0L, envir = .carryPairsMemo)
  }
  .st
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
  # cheap UI-level exits before ui$foceiEtaS builds a full symengine
  # environment (~0.25 s per call): no linCmt() or no covariate means no
  # pair can exist, and the empty result is identical
  if (rxode2::.rxLinNcmt(.ui)["numLin"] <= 0L) return(.rxFoceiCarryEmpty())
  if (length(.ui$allCovs) == 0L) return(.rxFoceiCarryEmpty())
  interpolation <- match.arg(interpolation)
  # the data-independent candidate result is a pure function of the model
  # digest; skip the foceiEtaS symengine build on a repeat fit of the same
  # model.  Data-dependent checks never reach this memo (data = NULL only).
  .key <- if (is.null(data) &&
              !identical(Sys.getenv("NLMIXR2EST_CARRY_MEMO"), "off")) {
    tryCatch(rxUiGet.foceiModelDigest(list(.ui)), error = function(e) NULL)
  }
  if (!is.null(.key) && exists(.key, envir = .carryPairsMemo, inherits = FALSE)) {
    assign("hits", .foceiLinCmtCarryMemoStats()[["hits"]] + 1L, envir = .carryPairsMemo)
    return(get(.key, envir = .carryPairsMemo, inherits = FALSE))
  }
  .s <- .ui$foceiEtaS
  .etaVars <- paste0("ETA_", seq_len(.s$..maxEta), "_")
  .ret <- .rxFoceiLinCmtCarryEligible(list(.ui), .s, .etaVars, data = data, # nolint: object_usage_linter.
                                      interpolation = interpolation)
  if (!is.null(.key)) {
    assign("misses", .foceiLinCmtCarryMemoStats()[["misses"]] + 1L, envir = .carryPairsMemo)
    # bound the memo (a session rarely fits more than a handful of models)
    .keys <- setdiff(ls(envir = .carryPairsMemo, all.names = FALSE), c("hits", "misses"))
    if (length(.keys) >= 64L) {
      rm(list = .keys, envir = .carryPairsMemo)
    }
    assign(.key, .ret, envir = .carryPairsMemo)
  }
  .ret
}

#' One eligibility data.frame row (rxFromSE() is substitute()-based: it
#' deparses a non-character argument's EXPRESSION, so it must always be
#' handed the repr STRING; render=FALSE keeps the raw repr for callers that
#' render later themselves)
#' @noRd
.rxFoceiCarryPairRow <- function(eta, etaName, slot, jump, mods, why, varying, render) {
  .txt <- function(x) .rxFoceiCarryTxt(x, render) # nolint: object_usage_linter.
  data.frame(slot = if (is.null(slot)) NA_integer_ else slot$k,
             slotName = if (is.null(slot)) NA_character_ else .rxFoceiLinCmtCarrySlotNames[slot$k], # nolint: object_usage_linter.
             eta = eta,
             etaName = etaName,
             covs = paste(why, collapse = ","),
             shape = if (is.null(slot)) NA_character_ else slot$shape,
             formula = if (is.null(slot)) NA_character_ else .txt(paste(slot$expr)),
             dEtaFormula = if (is.null(slot)) NA_character_ else .txt(paste(slot$g)),
             varying = varying,
             fD = .txt(jump$fD),
             fCov = isTRUE(jump$fCov),
             fCmt = if (is.null(jump$fD)) NA_character_ else mods$f$cmt,
             lagD = .txt(jump$lagD),
             lagCmt = if (is.null(jump$lagD)) NA_character_ else mods$lag$cmt,
             stringsAsFactors = FALSE)
}

#' Expected data CMT values for a linCmt() compartment name
#' @noRd
.rxFoceiCarryCmtValues <- function(cmt, oral0) {
  if (identical(cmt, "depot")) return(list(num = 1L, chr = "depot"))
  if (identical(cmt, "central")) return(list(num = oral0 + 1L, chr = "central"))
  list(num = NA_integer_, chr = cmt)
}

#' Do every dose's CMT values point at the modified compartment?
#' @noRd
.rxFoceiCarryDoseCmtOk <- function(dose, cmt, oral0) {
  if (!("CMT" %in% names(dose))) return(TRUE)
  .v <- dose[["CMT"]]
  .ok <- .rxFoceiCarryCmtValues(cmt, oral0)
  if (is.numeric(.v)) return(all(.v == .ok$num))
  all(as.character(.v) == .ok$chr)
}

#' Does the regimen contain an infusion (RATE/DUR)?
#' @noRd
.rxFoceiCarryHasInfusion <- function(dose) {
  .inf <- FALSE
  if ("RATE" %in% names(dose)) .inf <- .inf || any(dose[["RATE"]] != 0, na.rm = TRUE)
  if ("DUR" %in% names(dose)) .inf <- .inf || any(dose[["DUR"]] != 0, na.rm = TRUE)
  .inf
}

#' Dose rows (EVID 1/4/101) of an event data.frame with upper-case names,
#' or NULL when there is nothing to check
#' @noRd
.rxFoceiCarryDoseRows <- function(data) {
  .d <- data
  names(.d) <- toupper(names(.d))
  if (!("EVID" %in% names(.d))) return(NULL)
  .dose <- .d[.d[["EVID"]] %in% c(1L, 4L, 101L), , drop = FALSE]
  if (nrow(.dose) == 0L) return(NULL)
  .dose
}

#' Fit-time data check for jump pairs: every dose must enter the modified
#' compartment, and a lag pair or a covariate-driven F pair needs a
#' bolus-only regimen.  Returns NULL when fine, else a short reason.
#' @noRd
.rxFoceiCarryJumpDataProblem <- function(pairs, data, oral0) {
  .jump <- pairs[!is.na(pairs$fD) | !is.na(pairs$lagD), , drop = FALSE]
  if (nrow(.jump) == 0L) return(NULL)
  .dose <- .rxFoceiCarryDoseRows(data)
  if (is.null(.dose)) return(NULL)
  .cmts <- unique(c(.jump$fCmt, .jump$lagCmt))
  .cmts <- .cmts[!is.na(.cmts)]
  if (length(.cmts) > 1L) return("f() and alag() on different compartments")
  if (!.rxFoceiCarryDoseCmtOk(.dose, .cmts, oral0)) {
    return("doses outside the f()/alag() compartment")
  }
  .needBolus <- any(!is.na(.jump$lagD)) || any(isTRUE(.jump$fCov))
  if (.needBolus && .rxFoceiCarryHasInfusion(.dose)) {
    return("infusions with alag()/covariate f()")
  }
  NULL
}
