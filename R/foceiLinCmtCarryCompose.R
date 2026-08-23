# The per-model context and the final gradient composition of the linCmt()
# carry emission (foceiLinCmtCarry.R holds the per-row lines).  The context
# is parsed once per pair row from rx_pred_'s structural call
# (foceiLinCmtCarryPredCall.R); the composition turns the carried
# compartment sensitivity into d(rx_pred_)/d(eta) -- directly for a bare
# prediction, through the outer chain rule around rx_lcConc_ for a wrapped
# ll() likelihood (#1004).

#' Shared per-model context for the carry emission (parsed once per pair
#' row; cheap)
#' @noRd
.rxFoceiLinCmtCarryCtx <- function(pairs, s) {
  .pc <- .rxFoceiCarryPredCall(s) # nolint: object_usage_linter.
  .a <- .pc$args
  # rxFromSE() deparses a non-character argument's expression
  # (substitute()-based), so hand it the repr string
  .txt <- vapply(seq_along(.a), function(i) {
    .r <- paste(.a[[i]])
    rxode2::rxFromSE(.r)
  }, character(1))
  .callRepr <- paste(.pc$call)
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
       rows = seq_len(.m) - 1L,
       # an ll() endpoint embeds the concentration call in a larger
       # expression (#1004): the carry differentiates the concentration,
       # read back as rx_lcConc_, and symengine supplies the outer chain rule
       bare = .pc$bare, predSym = .pc$predSym,
       conc = if (.pc$bare) "rx_pred_" else "rx_lcConc_",
       concLine = if (.pc$bare) character(0) else
         paste0("rx_lcConc_~", rxode2::rxFromSE(.callRepr)))
}

#' The pair's final `dfe=` composition: s_central / Vc plus, when the
#' slot enters the observation scaling, the row-local direct term
#' -(dVc/dslot) g rx_pred_ / Vc (Vc from rxode2's own micro-constant
#' builder for this trans); a wrapped prediction multiplies the outer
#' derivative in and adds the direct eta term of the likelihood
#' @noRd
.rxFoceiLinCmtCarryFinal <- function(cx, pairs, w, dfe) {
  .p <- w - 1L
  .sC <- paste0("rx_lcCarryS", .p, "r", cx$central, "_")
  .vc <- .rxFoceiCarryVc(cx$ncmt, cx$oral0, as.numeric(cx$trans), pairs$slot[w], # nolint: object_usage_linter.
                         cx$slotExpr)
  .vcRepr <- paste(.vc$vc)
  .vcTxt <- paste0("(", rxode2::rxFromSE(.vcRepr), ")")
  .conc <- paste0(.sC, "/", .vcTxt)
  if (!is.null(.vc$dVc)) {
    .dRepr <- paste(.vc$dVc)
    .dTxt <- rxode2::rxFromSE(.dRepr)
    # lead with the symbol: rxode2's expression optimizer folds "-(1)*x"
    # to "-x", and "-" + "-x" would read as a C decrement
    .conc <- paste0(.conc, "-rx_lcCarryG", .p, "_*(", .dTxt, ")*", cx$conc, "/", .vcTxt)
  }
  if (cx$bare) return(paste0(dfe, "=", .conc))
  # ll() endpoint: d(pred)/d(eta) = d(pred)/d(conc) * d(conc)/d(eta) plus
  # whatever eta dependence the expression has with the concentration held
  # fixed (symengine sees rx_lcConc_ as a free symbol for that)
  .outer <- paste(symengine::D(cx$predSym, symengine::S("rx_lcConc_")))
  .direct <- symengine::D(cx$predSym, symengine::S(pairs$eta[w]))
  .fin <- paste0(dfe, "=(", rxode2::rxFromSE(.outer), ")*(", .conc, ")")
  if (!.rxFoceiCarryIsZero(.direct)) { # nolint: object_usage_linter.
    .directRepr <- paste(.direct)
    .fin <- paste0(.fin, "+(", rxode2::rxFromSE(.directRepr), ")")
  }
  .fin
}
