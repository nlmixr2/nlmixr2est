# Codegen substitution for the linCmt() sensitivity carry (phase 3b.3).
#
# For a carry-eligible (linCmt parameter slot, eta) pair (see
# foceiLinCmtCarryEligible.R) the ordinary symengine expansion of
# d(rx_pred_)/d(eta) is WRONG whenever a covariate makes that parameter vary
# within a subject: linCmtB()'s production Jacobian reconstructs the carried
# Alast sensitivity assuming theta is constant across the subject.  This
# file replaces the naive generated line for such a pair with lines that
# drive the exact carry recurrence through linCmtB()'s which1=-5/-6/-7
# sentinels (rxode2's per-subject ind->linCmtCarryT storage):
#
#   s_i = M_i s_{i-1} + g_i J_local_i(:, k)
#
# using ONLY quantities the generated model can read back per row:
#   - M_i s_{i-1}: the which1=-5 advance (applies this row's transition
#     matrix to every carry column at once; emitted once per row).
#   - J_local_i(:, k) = J_i(:, k) - M_i J_{i-1}(:, k): the production
#     cumulative amount Jacobian satisfies this identity BY CONSTRUCTION of
#     its constant-theta reconstruction, so the local Jacobian is recovered
#     from the restorable J_i (which1=row, which2=col) and a second
#     dedicated "tracker" carry column that holds J_{i-1}(:, k) -- after the
#     -5 advance the tracker reads back as exactly M_i J_{i-1}(:, k).
#   - g_i = d(theta_slot)/d(eta): the pair's dEtaFormula, an ordinary local
#     expression.
#
# Each eligible pair therefore consumes TWO carry columns (the carry itself
# and its tracker), so at most RX_LINCMT_CARRY_MAXPAIRS/2 = 4 pairs fit; the
# cap fails loudly at model build.  calc_lhs fires exactly once per event
# row in solve order (dose rows included, output filtering notwithstanding),
# which is what makes lhs-driven stepping sound; the -5 advance derives its
# interval from its own previous invocation time (rxode2 side) because
# ind->tprior is stale in the post-solve lhs pass.
#
# Known limitations (documented, all bias to the status quo):
# - Steady-state (ss) dose rows: the carry recurrence does not model the
#   SS fixed-point reset, so a regimen with ss doses on a carry-eligible
#   model can still be wrong -- same build-time-cannot-see-data class as
#   rxode2#1236/#1237; foceiControl(linCmtSensCarry="none") is the opt-out.
# - covsInterpolation="linear": cannot be represented by linCmt()'s
#   one-sample-per-row evaluation; the build skips substitution, and the
#   fit path (.foceiFamilyReturn) errors when the data confirms the
#   covariate actually varies within a subject.
# - An eta that also drives a modeled alag()/f() correction (#920) keeps
#   the status quo path (the substitution would drop that extra term).
# - The analytic 2nd-order inner Hessian (fast=TRUE ll()) has no
#   second-order carry; models with eligible pairs keep the Shi21
#   finite-difference Hessian fallback (gated in .foceiMaybeAddHdEta2).

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
  # the emission renders in-loop via rxFromSE(S(<repr>)) instead
  tryCatch(.rxFoceiLinCmtCarryEligible(x, s, etaVars, data = NULL,
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

#' Render one linCmtB() model-text call for the carry lines
#' @noRd
.rxFoceiLinCmtCarryCall <- function(pfx, which1, which2, trans, thetas) {
  paste0("linCmtB(", pfx, ",", which1, ",", which2, ",", trans, ",",
         paste(thetas, collapse = ","), ")")
}

#' Substituted d(rx_pred_)/d(eta) lines for one carry-eligible pair
#'
#' Replaces the naive generated gradient line for pair `w` with the full
#' carry block: the once-per-row -5 advance (first pair only), the pair's
#' g/J/tracker/carry/tracker-update intermediates (all `~`, evaluation
#' order is load-bearing), and the final `dfe=` composition
#' `s_central/v1` plus, for a v1-slot pair, the row-local direct term
#' `-g*rx_pred_/v1` (conc = A_central/v1 depends on v1 directly, not just
#' through the amounts; verified against FD in the 3b.2 bench).
#'
#' @param pairs data.frame from `.rxFoceiLinCmtCarryPairsForBuild()`
#' @param w 1-based row of `pairs` being emitted
#' @param s symengine env holding `rx_pred_`
#' @param dfe generated gradient lhs name (rx__sens_rx_pred__BY_ETA_n___)
#' @return multi-line model text (newline-joined) replacing the naive line
#' @noRd
.rxFoceiLinCmtCarryEmit <- function(pairs, w, s, dfe) {
  .pred <- get("rx_pred_", envir = s, inherits = FALSE)
  .a <- symengine::get_args(.pred)
  # rxFromSE() deparses a non-character argument's expression
  # (substitute()-based), so hand it the repr string
  .txt <- vapply(seq_along(.a), function(i) {
    .r <- paste(.a[[i]])
    rxode2::rxFromSE(.r)
  }, character(1))
  .ncmt <- as.integer(as.numeric(.txt[4]))
  .oral0 <- as.integer(as.numeric(.txt[5]))
  .m <- .ncmt + .oral0
  .central <- .oral0
  .trans <- .txt[8]
  # args 1-5 (ptr, t, linCmt, ncmt, oral0) are shared by every call
  .pfx <- paste(.txt[1:5], collapse = ",")
  .thetas <- .txt[9:15]
  .zero <- rep("0", 7)
  .nP <- nrow(pairs)
  .p <- w - 1L # 0-based pair column; tracker column is .nP + .p
  .slot <- pairs$slot[w]
  .kcol <- if (.slot == 7L) 2L*.ncmt else .slot - 1L
  .g <- paste0("rx_lcCarryG", .p, "_")
  .l <- character(0)
  if (w == 1L) {
    # once-per-row advance of every carry column (M depends only on this
    # row's theta, so one call serves every pair AND every tracker)
    .l <- c(.l, paste0("rx_lcCarryAdv_~",
                       .rxFoceiLinCmtCarryCall(.pfx, -5L, 0L, .trans, .thetas)))
  }
  # dEtaFormula is stored as symengine repr (render=FALSE); render to model
  # text here.  rxFromSE() deparses its argument EXPRESSION, so it must see
  # a bare symbol bound to the string, never an indexing expression.
  .gRepr <- pairs$dEtaFormula[w]
  .gTxt <- rxode2::rxFromSE(.gRepr)
  .l <- c(.l, paste0(.g, "~", .gTxt))
  .sC <- ""
  for (.r in seq_len(.m) - 1L) {
    .j <- paste0("rx_lcCarryJ", .p, "r", .r, "_")
    .pv <- paste0("rx_lcCarryP", .p, "r", .r, "_")
    .sv <- paste0("rx_lcCarryS", .p, "r", .r, "_")
    .uv <- paste0("rx_lcCarryU", .p, "r", .r, "_")
    .dloc <- paste0("(", .j, "-", .pv, ")")
    .z7 <- .zero
    .l <- c(.l,
            paste0(.j, "~", .rxFoceiLinCmtCarryCall(.pfx, .r, .kcol, .trans, .z7)),
            paste0(.pv, "~", .rxFoceiLinCmtCarryCall(.pfx, -6L, .r + .m*(.nP + .p),
                                                     .trans, .z7)))
    .z7[3] <- paste0(.g, "*", .dloc) # -7's added value rides in the p2 slot
    .l <- c(.l, paste0(.sv, "~", .rxFoceiLinCmtCarryCall(.pfx, -7L, .r + .m*.p,
                                                         .trans, .z7)))
    .z7[3] <- .dloc
    .l <- c(.l, paste0(.uv, "~", .rxFoceiLinCmtCarryCall(.pfx, -7L,
                                                         .r + .m*(.nP + .p),
                                                         .trans, .z7)))
    if (.r == .central) .sC <- .sv
  }
  .v1 <- paste0("(", .txt[10], ")")
  .fin <- paste0(dfe, "=", .sC, "/", .v1)
  if (.slot == 2L) {
    # conc = A_central/v1: the amounts carry misses the direct d(1/v1)/d(eta)
    .fin <- paste0(.fin, "-(", .g, ")*rx_pred_/", .v1)
  }
  paste(c(.l, .fin), collapse = "\n")
}
