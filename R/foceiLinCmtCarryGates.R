# Build-time gates for the linCmt() sensitivity carry: which of the detected
# carry-eligible pairs (foceiLinCmtCarryEligible.R) actually get their
# generated gradient line substituted (foceiLinCmtCarry.R).

#' Does the loaded rxode2 have the carry sentinels (which1=-5/-6/-7)?
#'
#' Detected via the linCmtCarryLiveTest export added with them; an rxode2
#' without the sentinels would mis-dispatch which1=-5 into the structural
#' solve, so nothing may be emitted without this.
#' @noRd
.rxFoceiLinCmtCarryCapable <- function() {
  exists("linCmtCarryLiveTest", envir = asNamespace("rxode2"), inherits = FALSE)
}

#' Does the loaded rxode2 have the fast-path pin (which1=-8) an event jump
#' needs?  Detected via linCmtCarrySentinelMax(), added with it.
#' @noRd
.rxFoceiLinCmtCarryJumpCapable <- function() {
  .ns <- asNamespace("rxode2")
  if (!exists("linCmtCarrySentinelMax", envir = .ns, inherits = FALSE)) {
    return(FALSE)
  }
  isTRUE(get("linCmtCarrySentinelMax", envir = .ns)() >= 8L)
}

#' Carry columns the loaded rxode2 has
#'
#' rxode2's `RX_LINCMT_CARRY_MAXPAIRS`, read through its
#' `linCmtCarryMaxPairs()`.  8 when rxode2 has no such function, or when it
#' fails (an rxode2lincmt without `.linCmtCarryMaxPairs()` indexes 8 columns
#' whatever rxode2 has) or gives no count; this runs inside the model
#' digest, so it never errors.  Remove the fallback with
#' test-dependency-workarounds.R's carry-cap block.
#' @return integer
#' @noRd
.rxFoceiLinCmtCarryMaxPairs <- function() {
  .ns <- asNamespace("rxode2")
  if (!exists("linCmtCarryMaxPairs", envir = .ns, inherits = FALSE)) {
    return(8L)
  }
  .cap <- try(get("linCmtCarryMaxPairs", envir = .ns)(), silent = TRUE)
  if (inherits(.cap, "try-error") || length(.cap) != 1L || !is.numeric(.cap) || is.na(.cap) || .cap < 1) {
    return(8L)
  }
  as.integer(.cap)
}

#' The carry-cap warning: the first eta past the cap and how many more
#'
#' Kept to one `$runInfo` line (under 75 characters, no method prefix).
#' @param etas the etas that keep the standard gradient
#' @param cap the carry columns the loaded rxode2 has
#' @noRd
.rxFoceiLinCmtCarryCapMsg <- function(etas, cap) {
  .msg <- if (length(etas) == 1L) {
    sprintf("%s has an approximate gradient (%d carry columns)", etas, cap)
  } else {
    sprintf(
      "%s and %d more have an approximate gradient (%d carry columns)",
      etas[1L],
      length(etas) - 1L,
      cap
    )
  }
  if (nchar(.msg) < 75L) {
    return(.msg)
  }
  # a long eta name: the count alone
  sprintf(
    "%d %s an approximate gradient (%d carry columns)",
    length(etas),
    if (length(etas) == 1L) "eta has" else "etas have",
    cap
  )
}

#' Record a carry-cap warning for the model being built
#'
#' The warning is given while the focei model builds; rxUiGet.foceiModel()
#' stores what was recorded with the cached bundle and gives it again when a
#' later fit reads the bundle from the cache.
#' @param msg the warning text
#' @noRd
.rxFoceiLinCmtCarryCapWarn <- function(msg) {
  nlmixr2global$carryCapWarnings <- c(nlmixr2global$carryCapWarnings, msg)
  warning(msg, call. = FALSE)
}

#' Carry-eligible pairs, gated for actual codegen substitution
#'
#' Applies every build-time gate on top of `.rxFoceiLinCmtCarryEligible()`:
#' rxode2 capability (the jump pin for pairs with an f()/alag() channel),
#' the `linCmtSensCarry` control ("auto"/"none"), interpolation, the
#' parameterization (rxode2's own linToOde micro-constant builder must know
#' the trans, which is where the volume and system matrix come from), and
#' the carry-column cap.  Detection is data-independent (candidate based)
#' so the compiled-model cache stays coherent across datasets; a
#' constant-in-data covariate still produces a correct (identical) gradient
#' through the carry, just at extra cost (the 3b.4 runtime fast path is the
#' optimization for that).
#'
#' @param x list(rxUi)
#' @param s FOCEi symengine env holding `rx_pred_`
#' @param etaVars ETA_i_ names in gradient-column order
#' @param extraPred unused (kept for the call site); a pair's own f()/alag()
#'   channel replaces the #920 row-local term for that eta
#' @return eligibility data.frame (zero rows -> NULL), or NULL
#' @noRd
.rxFoceiLinCmtCarryPairsForBuild <- function(x, s, etaVars, extraPred = NULL) {
  if (!.rxFoceiLinCmtCarryBuildEnabled(x[[1]])) {
    return(NULL)
  }
  .pairs <- .rxFoceiLinCmtCarryDetect(x, s, etaVars)
  if (is.null(.pairs) || nrow(.pairs) == 0L) {
    return(NULL)
  }
  if (!.rxFoceiLinCmtCarryModelOk(s)) {
    return(NULL)
  }
  .pairs <- .rxFoceiLinCmtCarryJumpFilter(.pairs)
  if (is.null(.pairs)) {
    return(NULL)
  }
  .pairs <- .rxFoceiLinCmtCarryCapPairs(.pairs)
  .rxFoceiLinCmtCarryCheckCap(.pairs)
  attr(.pairs, "carryCovInterval") <- .rxFoceiLinCmtCarryCovInterval(
    x[[1]]$allCovs,
    .rxFoceiLinCmtCarryInterp(x[[1]], x[[1]]$allCovs)
  )
  .pairs
}

#' Drop the pairs a jump channel or multCov shape needs the fast-path pin
#' (-8) for when the loaded rxode2 lacks it
#' @param pairs the eligible pairs
#' @return the pairs that can be carried, or NULL when none can
#' @noRd
.rxFoceiLinCmtCarryJumpFilter <- function(pairs) {
  if (is.null(pairs) || nrow(pairs) == 0L) {
    return(NULL)
  }
  .jump <- !is.na(pairs$fD) | !is.na(pairs$lagD) | pairs$shape %in% "multCov"
  if (any(.jump) && !.rxFoceiLinCmtCarryJumpCapable()) {
    pairs <- pairs[!.jump, , drop = FALSE]
    if (nrow(pairs) == 0L) {
      return(NULL)
    }
  }
  pairs
}

#' Is the carry switched on for this build (capability, control, interpolation)?
#' @noRd
.rxFoceiLinCmtCarryBuildEnabled <- function(ui) {
  # covsInterpolation rides inside control$rxControl; a bare rxGetControl()
  # lookup never sees it
  .rxFoceiLinCmtCarryOn(
    rxode2::rxGetControl(ui, "linCmtSensCarry", "auto"),
    rxode2::rxGetControl(ui, "rxControl", NULL)$covsInterpolation
  )
}

#' Is the carry switched on for a `linCmtSensCarry` value and covariate
#' interpolation (0=linear, 1=locf, 2=nocb, 3=midpoint)?
#' @param carry the `linCmtSensCarry` control
#' @param interp the `covsInterpolation` control
#' @noRd
.rxFoceiLinCmtCarryOn <- function(carry, interp) {
  if (!.rxFoceiLinCmtCarryCapable()) {
    return(FALSE)
  }
  if (identical(carry, "none")) {
    return(FALSE)
  }
  !(identical(suppressWarnings(as.integer(interp)), 0L) || identical(interp, "linear"))
}

#' Run the eligibility detection without rendering; a failure warns and
#' keeps the status quo gradient
#' @noRd
.rxFoceiLinCmtCarryDetect <- function(x, s, etaVars) {
  # render=FALSE: no rxFromSE may run before rxUiGet.foceiHdEta's own
  # D(rx_pred_, ETA_n_) evaluations (shared-symengine-state corruption);
  # the emission renders in-loop via rxFromSE(S(<repr>)) instead.
  # (nolint: lintr resolves cross-file helpers against the installed package)
  tryCatch(
    .rxFoceiLinCmtCarryEligible(
      x,
      s,
      etaVars,
      data = NULL,
      render = FALSE
    ),
    error = function(e) {
      warning("linCmt() carry detection failed; standard gradient used", call. = FALSE)
      NULL
    }
  )
}

#' ncmt / oral0 / trans of the structural linCmtB() call in `rx_pred_`
#' @noRd
.rxFoceiLinCmtCarryShape <- function(s) {
  .a <- .rxFoceiCarryPredCall(s)$args
  .num <- function(i) suppressWarnings(as.numeric(paste(.a[[i]])))
  list(ncmt = as.integer(.num(4)), oral0 = as.integer(.num(5)), trans = .num(8))
}

#' The parameterization must be one rxode2's own linToOde translation knows
#' (that is where Vc and the system matrix come from)
#' @noRd
.rxFoceiLinCmtCarryModelOk <- function(s) {
  .sh <- .rxFoceiLinCmtCarryShape(s)
  if (any(is.na(unlist(.sh)))) {
    return(FALSE)
  }
  !is.null(.rxFoceiCarryMicro(.sh$ncmt, .sh$oral0, .sh$trans))
}

#' Carry columns the pairs need
#'
#' One per pair, one tracker per slot (`.rxFoceiLinCmtCarryTrackers()`), plus
#' one shared amounts tracker when any pair has a jump channel and one shared
#' lag tracker when any pair has an alag() channel
#' @noRd
.rxFoceiLinCmtCarryColumns <- function(pairs) {
  nrow(pairs) +
    .rxFoceiLinCmtCarryTrackers(pairs$slot)$n +
    as.integer(any(!is.na(pairs$fD) | !is.na(pairs$lagD))) +
    as.integer(any(!is.na(pairs$lagD)))
}

#' Drop the pairs past the carry's column cap, the newest kinds first
#'
#' A fit never gets here past the cap: `.uiApplyLinCmtCarryCap()` solves
#' such a model as ODEs.  A build that does (a direct `ui$foceiModel`, or a
#' model `linToOde()` cannot translate) drops the multCov pairs (an IOV
#' rewrite's occasion etas), last first, and then the other pairs that share
#' a slot.  The dropped etas keep the standard gradient, which is wrong for a
#' within-subject covariate and makes the objective depend on eta order
#' (#1318), so the warning says it is approximate; a model still past the cap
#' stops in `.rxFoceiLinCmtCarryCheckCap()`.
#' @param pairs the eligible pairs
#' @return `pairs` without the dropped rows
#' @noRd
.rxFoceiLinCmtCarryCapPairs <- function(pairs) {
  .cap <- .rxFoceiLinCmtCarryMaxPairs()
  if (.rxFoceiLinCmtCarryColumns(pairs) <= .cap) {
    return(pairs)
  }
  .shared <- !is.na(pairs$slot) &
    (duplicated(pairs$slot) | duplicated(pairs$slot, fromLast = TRUE))
  .order <- c(
    rev(which(pairs$shape %in% "multCov")),
    rev(which(.shared & !(pairs$shape %in% "multCov")))
  )
  .keep <- rep(TRUE, nrow(pairs))
  for (.w in .order) {
    if (.rxFoceiLinCmtCarryColumns(pairs[.keep, , drop = FALSE]) <= .cap) {
      break
    }
    .keep[.w] <- FALSE
  }
  if (any(!.keep) && .rxFoceiLinCmtCarryColumns(pairs[.keep, , drop = FALSE]) <= .cap) {
    .rxFoceiLinCmtCarryCapWarn(.rxFoceiLinCmtCarryCapMsg(pairs$etaName[!.keep], .cap))
    pairs <- pairs[.keep, , drop = FALSE]
  }
  pairs
}

#' @noRd
.rxFoceiLinCmtCarryCheckCap <- function(pairs) {
  .need <- .rxFoceiLinCmtCarryColumns(pairs)
  .cap <- .rxFoceiLinCmtCarryMaxPairs()
  if (.need > .cap) {
    stop(
      "the carry-eligible (linCmt parameter, eta) pairs need ",
      .need,
      " carry columns (",
      .cap,
      " available); reduce the model or use ",
      "foceiControl(linCmtSensCarry=\"none\")",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
