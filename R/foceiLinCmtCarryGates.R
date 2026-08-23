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
  exists("linCmtCarryLiveTest", envir = asNamespace("rxode2"),
         inherits = FALSE)
}

#' Carry-eligible pairs, gated for actual codegen substitution
#'
#' Applies every build-time gate on top of `.rxFoceiLinCmtCarryEligible()`:
#' rxode2 capability, the `linCmtSensCarry` control ("auto"/"none"),
#' interpolation, trans (the final composition needs `conc = A_central/v1`,
#' verified for trans 1 and 2), the #920 alag/f interaction, and the
#' two-columns-per-pair cap.  Detection is data-independent (candidate
#' based) so the compiled-model cache stays coherent across datasets; a
#' constant-in-data covariate still produces a correct (identical) gradient
#' through the carry, just at extra cost (the 3b.4 runtime fast path is the
#' planned optimization).
#'
#' @param x list(rxUi)
#' @param s FOCEi symengine env holding `rx_pred_`
#' @param etaVars ETA_i_ names in gradient-column order
#' @param extraPred `.rxFoceiLinCmtEventPredExtra()` result (or NULL)
#' @return eligibility data.frame (zero rows -> NULL), or NULL
#' @noRd
.rxFoceiLinCmtCarryPairsForBuild <- function(x, s, etaVars, extraPred = NULL) {
  if (!.rxFoceiLinCmtCarryBuildEnabled(x[[1]])) return(NULL)
  .pairs <- .rxFoceiLinCmtCarryDetect(x, s, etaVars)
  if (is.null(.pairs) || nrow(.pairs) == 0L) return(NULL)
  if (!.rxFoceiLinCmtCarryTransOk(s)) return(NULL)
  if (!is.null(extraPred)) {
    .pairs <- .pairs[!(.pairs$eta %in% names(extraPred)), , drop = FALSE]
    if (nrow(.pairs) == 0L) return(NULL)
  }
  .rxFoceiLinCmtCarryCheckCap(.pairs)
  .pairs
}

#' Is the carry switched on for this build (capability, control, interpolation)?
#' @noRd
.rxFoceiLinCmtCarryBuildEnabled <- function(ui) {
  if (!.rxFoceiLinCmtCarryCapable()) return(FALSE)
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
  tryCatch(.rxFoceiLinCmtCarryEligible(x, s, etaVars, data = NULL, # nolint: object_usage_linter.
                                       render = FALSE),
           error = function(e) {
             warning("linCmt() carry detection failed; standard gradient used",
                     call. = FALSE)
             NULL
           })
}

#' The v1 direct term is only derived for trans 1/2
#' @noRd
.rxFoceiLinCmtCarryTransOk <- function(s) {
  .pred <- get("rx_pred_", envir = s, inherits = FALSE)
  .trans <- suppressWarnings(as.numeric(paste(symengine::get_args(.pred)[[8]])))
  isTRUE(.trans %in% c(1, 2))
}

#' @noRd
.rxFoceiLinCmtCarryCheckCap <- function(pairs) {
  .max <- 4L # two carry columns per pair; RX_LINCMT_CARRY_MAXPAIRS = 8
  if (nrow(pairs) > .max) {
    stop("more than ", .max, " carry-eligible (linCmt parameter, eta) pairs (",
         nrow(pairs), "); reduce the model or use foceiControl(linCmtSensCarry=\"none\")",
         call. = FALSE)
  }
  invisible(TRUE)
}
