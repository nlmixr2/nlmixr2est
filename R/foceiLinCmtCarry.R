# Codegen substitution for the linCmt() sensitivity carry (phase 3b.3).
#
# For a carry-eligible (eta, channels) pair (see foceiLinCmtCarryEligible.R,
# gated in foceiLinCmtCarryGates.R) the ordinary symengine expansion of
# d(rx_pred_)/d(eta) is WRONG whenever a covariate makes a linCmt()
# parameter -- or a modeled f() -- vary within a subject: linCmtB()'s
# production Jacobian reconstructs the carried Alast sensitivity assuming
# theta is constant across the subject, and #920's f()/alag() terms are
# row-local under the same assumption.  This file replaces the naive
# generated line for such a pair with lines that drive the exact carry
# recurrence through linCmtB()'s which1=-5/-6/-7 sentinels (rxode2's
# per-subject ind->linCmtCarryT storage):
#
#   s_i = M_i s_{i-1} + g_i J_local_i(:, k) + jump_i
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
#   - jump_i: the f()/alag() contribution of the input that entered during
#     the row's interval, built from an amounts tracker the same way
#     (foceiLinCmtCarryEvent.R); it does not telescope, so the row is
#     pinned to the full advance (which1=-8).
#
# Each eligible pair therefore consumes TWO carry columns (the carry itself
# and its tracker), plus one shared amounts column when any pair has a jump
# channel and one shared lag column when any has an alag() channel; the cap
# fails loudly at model build.  calc_lhs fires exactly once per event row in
# solve order (dose rows included, a lagged dose on its own shifted row,
# output filtering notwithstanding), which is what makes lhs-driven stepping
# sound; the -5 advance derives its interval from its own previous
# invocation time (rxode2 side) because ind->tprior is stale in the
# post-solve lhs pass.
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
# - Jump channels need every dose to enter the modified compartment, and
#   an alag() or covariate-driven f() channel a bolus-only regimen (checked
#   against the data in .foceiFamilyReturn); modeled rate()/dur() are not
#   channels at all.
# - The analytic 2nd-order inner Hessian (fast=TRUE ll()) has no
#   second-order carry; models with eligible pairs keep the Shi21
#   finite-difference Hessian fallback (gated in .foceiMaybeAddHdEta2).

#' Render one linCmtB() model-text call for the carry lines
#' @noRd
.rxFoceiLinCmtCarryCall <- function(pfx, which1, which2, trans, thetas) {
  paste0("linCmtB(", pfx, ",", which1, ",", which2, ",", trans, ",",
         paste(thetas, collapse = ","), ")")
}

#' Shared per-model context for the carry emission (parsed once per pair
#' row; cheap)
#' @noRd
.rxFoceiLinCmtCarryCtx <- function(pairs, s) {
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
  .nP <- nrow(pairs)
  .anyJump <- any(!is.na(pairs$fD) | !is.na(pairs$lagD))
  .anyLag <- any(!is.na(pairs$lagD))
  .m <- .ncmt + .oral0
  list(ncmt = .ncmt, oral0 = .oral0, m = .m, central = .oral0,
       trans = .txt[8], pfx = paste(.txt[1:5], collapse = ","),
       thetas = .txt[9:15], zero = rep("0", 7), nP = .nP,
       anyJump = .anyJump, anyLag = .anyLag,
       aCol = 2L * .nP, lCol = 2L * .nP + 1L,
       slotExpr = lapply(1:7, function(k) .a[[k + 8L]]),
       rows = seq_len(.m) - 1L)
}

#' One pair's per-row lines: g, the per-compartment J/tracker reads, the
#' -7 carry add (kernel + jump terms) and the tracker update
#' @noRd
.rxFoceiLinCmtCarryPairLines <- function(cx, pairs, w) {
  .p <- w - 1L # 0-based pair column; tracker column is nP + p
  .slot <- pairs$slot[w]
  .hasSlot <- !is.na(.slot)
  .l <- character(0)
  .g <- paste0("rx_lcCarryG", .p, "_")
  .f <- paste0("rx_lcCarryF", .p, "_")
  .lg <- paste0("rx_lcCarryLg", .p, "_")
  # the stored reprs (render=FALSE) are rendered here; rxFromSE() deparses
  # its argument EXPRESSION, so it must see a bare symbol bound to the string
  if (.hasSlot) {
    .gRepr <- pairs$dEtaFormula[w]
    .l <- c(.l, paste0(.g, "~", rxode2::rxFromSE(.gRepr)))
  }
  if (!is.na(pairs$fD[w])) {
    .fRepr <- pairs$fD[w]
    .l <- c(.l, paste0(.f, "~", rxode2::rxFromSE(.fRepr)))
  }
  if (!is.na(pairs$lagD[w])) {
    .lRepr <- pairs$lagD[w]
    .l <- c(.l, paste0(.lg, "~", rxode2::rxFromSE(.lRepr)))
  }
  .kcol <- if (isTRUE(.slot == 7L)) 2L * cx$ncmt else .slot - 1L
  for (.r in cx$rows) {
    .terms <- character(0)
    .z7 <- cx$zero
    if (.hasSlot) {
      .j <- paste0("rx_lcCarryJ", .p, "r", .r, "_")
      .pv <- paste0("rx_lcCarryP", .p, "r", .r, "_")
      .dloc <- paste0("(", .j, "-", .pv, ")")
      .l <- c(.l,
              paste0(.j, "~", .rxFoceiLinCmtCarryCall(cx$pfx, .r, .kcol, cx$trans, .z7)),
              paste0(.pv, "~", .rxFoceiLinCmtCarryCall(cx$pfx, -6L, .r + cx$m * (cx$nP + .p),
                                                       cx$trans, .z7)))
      .terms <- c(.terms, paste0(.g, "*", .dloc))
    }
    if (!is.na(pairs$fD[w])) {
      .terms <- c(.terms, paste0("rx_lcCarryD", .r, "_*", .f))
    }
    if (!is.na(pairs$lagD[w])) {
      .terms <- c(.terms, paste0("rx_lcCarryKP", .r, "_*", .lg,
                                 "*(rx_lcCarryDel_-rx_lcCarryDelP_)"))
    }
    .z7[3] <- paste(.terms, collapse = "+") # -7's added value rides in the p2 slot
    .l <- c(.l, paste0("rx_lcCarryS", .p, "r", .r, "_~",
                       .rxFoceiLinCmtCarryCall(cx$pfx, -7L, .r + cx$m * .p, cx$trans, .z7)))
    if (.hasSlot) {
      .z7[3] <- .dloc
      .l <- c(.l, paste0("rx_lcCarryU", .p, "r", .r, "_~",
                         .rxFoceiLinCmtCarryCall(cx$pfx, -7L, .r + cx$m * (cx$nP + .p),
                                                 cx$trans, .z7)))
    }
  }
  .l
}

#' The pair's final `dfe=` composition: s_central / Vc plus, when the
#' slot enters the observation scaling, the row-local direct term
#' -(dVc/dslot) g rx_pred_ / Vc (Vc from rxode2's own micro-constant
#' builder for this trans)
#' @noRd
.rxFoceiLinCmtCarryFinal <- function(cx, pairs, w, dfe) {
  .p <- w - 1L
  .sC <- paste0("rx_lcCarryS", .p, "r", cx$central, "_")
  .vc <- .rxFoceiCarryVc(cx$ncmt, cx$oral0, as.numeric(cx$trans), pairs$slot[w], # nolint: object_usage_linter.
                         cx$slotExpr)
  .vcRepr <- paste(.vc$vc)
  .vcTxt <- paste0("(", rxode2::rxFromSE(.vcRepr), ")")
  .fin <- paste0(dfe, "=", .sC, "/", .vcTxt)
  if (!is.null(.vc$dVc)) {
    .dRepr <- paste(.vc$dVc)
    .dTxt <- rxode2::rxFromSE(.dRepr)
    # lead with the symbol: rxode2's expression optimizer folds "-(1)*x"
    # to "-x", and "-" + "-x" would read as a C decrement
    .fin <- paste0(.fin, "-rx_lcCarryG", .p, "_*(", .dTxt, ")*rx_pred_/", .vcTxt)
  }
  .fin
}

#' Substituted d(rx_pred_)/d(eta) lines for one carry-eligible pair
#'
#' Replaces the naive generated gradient line for pair `w` with the full
#' carry block: the once-per-row prelude (first pair only), the pair's own
#' intermediates (all `~`, evaluation order is load-bearing), the tracker
#' epilogue (last pair only) and the final `dfe=` composition.
#'
#' @param pairs data.frame from `.rxFoceiLinCmtCarryPairsForBuild()`
#' @param w 1-based row of `pairs` being emitted
#' @param s symengine env holding `rx_pred_`
#' @param dfe generated gradient lhs name (rx__sens_rx_pred__BY_ETA_n___)
#' @return multi-line model text (newline-joined) replacing the naive line
#' @noRd
.rxFoceiLinCmtCarryEmit <- function(pairs, w, s, dfe) {
  .cx <- .rxFoceiLinCmtCarryCtx(pairs, s)
  .l <- character(0)
  if (w == 1L) .l <- c(.l, .rxFoceiLinCmtCarryPrelude(.cx)) # nolint: object_usage_linter.
  .l <- c(.l, .rxFoceiLinCmtCarryPairLines(.cx, pairs, w))
  if (w == .cx$nP) .l <- c(.l, .rxFoceiLinCmtCarryEpilogue(.cx)) # nolint: object_usage_linter.
  paste(c(.l, .rxFoceiLinCmtCarryFinal(.cx, pairs, w, dfe)), collapse = "\n")
}
