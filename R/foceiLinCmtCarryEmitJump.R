# The once-per-row prelude and epilogue of the linCmt() carry emission
# (foceiLinCmtCarry.R): the fast-path pin, the -5 advance, and -- when a
# pair has an f()/alag() channel -- the amounts/lag trackers and the
# dose-entry indicators the jump contributions read (foceiLinCmtCarryEvent.R
# for the math).

#' Lines emitted once per row before any pair (with the first pair)
#' @noRd
.rxFoceiLinCmtCarryPrelude <- function(cx) { # nolint: object_usage_linter.
  .l <- character(0)
  if (cx$anyJump) {
    # pin the full -5 advance: a jump contribution does not telescope
    .l <- c(.l, paste0("rx_lcCarryPin_~",
                       .rxFoceiLinCmtCarryCall(cx$pfx, -8L, 0L, cx$trans, cx$zero))) # nolint: object_usage_linter.
  }
  # once-per-row advance of every carry column (M depends only on this
  # row's theta, so one call serves every pair AND every tracker)
  .l <- c(.l, paste0("rx_lcCarryAdv_~",
                     .rxFoceiLinCmtCarryCall(cx$pfx, -5L, 0L, cx$trans, cx$thetas))) # nolint: object_usage_linter.
  if (!cx$anyJump) return(.l)
  for (.r in cx$rows) {
    .l <- c(.l,
            # amounts at this row (post-dose) and the advanced amounts
            # tracker M_r A_{r-1}; their difference is the interval's input
            paste0("rx_lcCarryA", .r, "_~",
                   .rxFoceiLinCmtCarryCall(cx$pfx, .r, -2L, cx$trans, cx$zero)), # nolint: object_usage_linter.
            paste0("rx_lcCarryPA", .r, "_~",
                   .rxFoceiLinCmtCarryCall(cx$pfx, -6L, .r + cx$m * cx$aCol, cx$trans, cx$zero)), # nolint: object_usage_linter.
            paste0("rx_lcCarryD", .r, "_~(rx_lcCarryA", .r, "_-rx_lcCarryPA", .r, "_)"))
  }
  if (!cx$anyLag) return(.l)
  for (.r in cx$rows) {
    .l <- c(.l, paste0("rx_lcCarryL", .r, "_~",
                       .rxFoceiLinCmtCarryCall(cx$pfx, -6L, .r + cx$m * cx$lCol, cx$trans, cx$zero))) # nolint: object_usage_linter.
  }
  # K_r P_r: the kernel's right-hand side applied to the advanced amounts
  .kp <- .rxFoceiCarryRhsTxt(cx$ncmt, cx$oral0, as.numeric(cx$trans), cx$slotExpr, # nolint: object_usage_linter.
                             paste0("rx_lcCarryPA", cx$rows, "_"))
  for (.r in cx$rows) .l <- c(.l, paste0("rx_lcCarryKP", .r, "_~", .kp[.r + 1L]))
  .absD <- paste(paste0("abs(rx_lcCarryD", cx$rows, "_)"), collapse = "+")
  .absL <- paste(paste0("abs(rx_lcCarryL", cx$rows, "_)"), collapse = "+")
  .scale <- paste(paste0("abs(rx_lcCarryA", cx$rows, "_)+abs(rx_lcCarryPA", cx$rows, "_)"),
                  collapse = "+")
  c(.l,
    # did an input enter at this row / at the previous row?  D is a
    # difference of two solves of the same state, so compare it against
    # the amounts' own scale, not against zero
    paste0("rx_lcCarryDel_~ifelse((", .absD, ")>1e-8*(", .scale, ")+1e-300,1,0)"),
    paste0("rx_lcCarryDelP_~ifelse((", .absL, ")>1e-8*(", .scale, ")+1e-300,1,0)"))
}

#' Lines emitted once per row after every pair (with the last pair): keep
#' the amounts / lag trackers at this row's values
#' @noRd
.rxFoceiLinCmtCarryEpilogue <- function(cx) { # nolint: object_usage_linter.
  if (!cx$anyJump) return(character(0))
  .l <- character(0)
  for (.r in cx$rows) {
    .z7 <- cx$zero
    .z7[3] <- paste0("rx_lcCarryD", .r, "_")
    .l <- c(.l, paste0("rx_lcCarryUA", .r, "_~",
                       .rxFoceiLinCmtCarryCall(cx$pfx, -7L, .r + cx$m * cx$aCol, cx$trans, .z7))) # nolint: object_usage_linter.
    if (cx$anyLag) {
      .z7[3] <- paste0("(rx_lcCarryD", .r, "_-rx_lcCarryL", .r, "_)")
      .l <- c(.l, paste0("rx_lcCarryUL", .r, "_~",
                         .rxFoceiLinCmtCarryCall(cx$pfx, -7L, .r + cx$m * cx$lCol, cx$trans, .z7))) # nolint: object_usage_linter.
    }
  }
  .l
}

