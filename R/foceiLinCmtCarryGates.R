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
  # a jump channel and a multCov pair both need the fast-path pin (-8)
  .jump <- !is.na(.pairs$fD) | !is.na(.pairs$lagD) | .pairs$shape %in% "multCov"
  if (any(.jump) && !.rxFoceiLinCmtCarryJumpCapable()) {
    .pairs <- .pairs[!.jump, , drop = FALSE]
    if (nrow(.pairs) == 0L) {
      return(NULL)
    }
  }
  .pairs <- .rxFoceiLinCmtCarryCapPairs(.pairs)
  .rxFoceiLinCmtCarryCheckCap(.pairs)
  attr(.pairs, "carryCovInterval") <- .rxFoceiLinCmtCarryCovInterval(
    x[[1]]$allCovs,
    .rxFoceiLinCmtCarryInterp(x[[1]], x[[1]]$allCovs)
  )
  .pairs
}

#' Is the carry switched on for this build (capability, control, interpolation)?
#' @noRd
.rxFoceiLinCmtCarryBuildEnabled <- function(ui) {
  if (!.rxFoceiLinCmtCarryCapable()) {
    return(FALSE)
  }
  if (identical(rxode2::rxGetControl(ui, "linCmtSensCarry", "auto"), "none")) {
    return(FALSE)
  }
  # covsInterpolation rides inside control$rxControl (0=linear, 1=locf,
  # 2=nocb, 3=midpoint); a bare rxGetControl() lookup never sees it
  .interp <- rxode2::rxGetControl(ui, "rxControl", NULL)$covsInterpolation
  !(identical(as.integer(.interp), 0L) || identical(.interp, "linear"))
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
#' rxode2 has RX_LINCMT_CARRY_MAXPAIRS = 8 carry columns.  Past them, the
#' multCov pairs (an IOV rewrite's occasion etas), last first, and then the
#' other pairs that share a slot keep the standard gradient, with one
#' message, so a model past the cap is carried no worse than without them; a
#' model still past it stops in `.rxFoceiLinCmtCarryCheckCap()`.
#' @param pairs the eligible pairs
#' @return `pairs` without the dropped rows
#' @noRd
.rxFoceiLinCmtCarryCapPairs <- function(pairs) {
  if (.rxFoceiLinCmtCarryColumns(pairs) <= 8L) {
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
    if (.rxFoceiLinCmtCarryColumns(pairs[.keep, , drop = FALSE]) <= 8L) {
      break
    }
    .keep[.w] <- FALSE
  }
  if (any(!.keep) && .rxFoceiLinCmtCarryColumns(pairs[.keep, , drop = FALSE]) <= 8L) {
    .minfo(paste0(
      "linCmt() carry: ",
      paste(pairs$etaName[!.keep], collapse = ", "),
      " keep the standard gradient (carry column limit)"
    ))
    pairs <- pairs[.keep, , drop = FALSE]
  }
  pairs
}

#' @noRd
.rxFoceiLinCmtCarryCheckCap <- function(pairs) {
  .need <- .rxFoceiLinCmtCarryColumns(pairs)
  if (.need > 8L) {
    stop(
      "the carry-eligible (linCmt parameter, eta) pairs need ",
      .need,
      " carry columns (8 available); reduce the model or use ",
      "foceiControl(linCmtSensCarry=\"none\")",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
