# Event-modifier (jump) sensitivities and parameterization support for the
# linCmt() sensitivity carry.
#
# The carry recurrence (foceiLinCmtCarry.R) covers every way an eta reaches
# the state THROUGH THE KERNEL PARAMETERS between rows.  An eta can also
# reach the state through an EVENT: a modeled f() scales the dose that
# enters a compartment, a modeled alag() moves the time it enters.  Both are
# invisible to the kernel Jacobian, and the row-local #920 terms
# (foceiLinCmtAlagSens.R) are exact only while theta is constant across the
# subject.  This file supplies the per-row jump contributions the carry
# needs instead, built from quantities the generated model can read back:
#
#   D_r = A_r - M_r A_{r-1}   the input that entered during row r's interval
#                             (amounts read minus the advanced amounts
#                             tracker), zero on a row with no dose
#   f():   s_r += D_r * d(ln F)/d(eta)       (every input is proportional
#                                            to F, bolus and infusion alike)
#   alag(): s_r += K_r P_r * (delta_r - delta_{r-1}) * d(lag)/d(eta)
#          with P_r = M_r A_{r-1} and K_r the row's system matrix: the
#          interval lengths on either side of a lagged dose row both move
#          with the lag, and the entry term -K Phi x0 is what those two
#          boundary terms sum to after the next advance (bolus only)
#
# Parameterizations: the reported concentration is A_central / Vc(theta)
# for every trans (rxode2's getVc), so the observation-scaling direct term
# is -(dVc/dslot) g rx_pred_ / Vc with Vc taken from rxode2's own
# linToOde micro-constant builders -- the same source the system matrix
# K comes from.  Nothing here is hand-derived.

.rxFoceiCarrySlotPh <- paste0("rx__lcSlot", 1:7, "__")

#' rxode2's micro constants (k, k12, k21, k13, k31, v, ka) for one
#' (ncmt, oral0, trans), as symengine expressions in the slot placeholders;
#' NULL for a translation rxode2 itself does not support
#' @noRd
.rxFoceiCarryMicro <- function(ncmt, oral0, trans) {
  .args <- list(ncmt = ncmt, oral0 = oral0, trans = trans)
  for (.i in 1:7) .args[[.rxFoceiLinCmtCarrySlotNames[.i]]] <- as.name(.rxFoceiCarrySlotPh[.i]) # nolint: object_usage_linter.
  .build <- utils::getFromNamespace(paste0(".linToOdeBuildMicro", ncmt), "rxode2")
  .m <- tryCatch(.build(.args), error = function(e) NULL)
  if (is.null(.m)) {
    return(NULL)
  }
  .keep <- intersect(c("k", "k12", "k21", "k13", "k31", "v", "ka"), names(.m))
  .ret <- lapply(.keep, function(n) {
    symengine::S(paste(deparse(.m[[n]], width.cutoff = 500L), collapse = ""))
  })
  names(.ret) <- .keep
  .ret
}

#' Substitute the slot placeholders with the model's slot expressions
#' @noRd
.rxFoceiCarrySubsSlots <- function(expr, slotExpr) {
  .free <- .rxFoceiCarryFreeSyms(expr) # nolint: object_usage_linter.
  for (.i in 1:7) {
    if (.rxFoceiCarrySlotPh[.i] %in% .free) {
      expr <- symengine::subs(expr, symengine::S(.rxFoceiCarrySlotPh[.i]), slotExpr[[.i]])
    }
  }
  expr
}

#' Vc and dVc/dslot (symengine, slots substituted) for a pair's slot;
#' dVc is NULL when the slot does not enter the observation scaling, the
#' whole result NULL for an unsupported translation
#' @noRd
.rxFoceiCarryVc <- function(ncmt, oral0, trans, slot, slotExpr) {
  .m <- .rxFoceiCarryMicro(ncmt, oral0, trans)
  if (is.null(.m) || is.null(.m$v)) {
    return(NULL)
  }
  .vc <- .m$v
  .dvc <- NULL
  if (!is.na(slot)) {
    .d <- symengine::D(.vc, symengine::S(.rxFoceiCarrySlotPh[slot]))
    if (!.rxFoceiCarryIsZero(.d)) .dvc <- .rxFoceiCarrySubsSlots(.d, slotExpr) # nolint: object_usage_linter.
  }
  list(vc = .rxFoceiCarrySubsSlots(.vc, slotExpr), dVc = .dvc)
}

#' The m rows of K x (the kernel's right-hand side applied to the vector
#' named by `xNames`), as rendered model text; NULL if unsupported
#' @noRd
.rxFoceiCarryRhsTxt <- function(ncmt, oral0, trans, slotExpr, xNames) {
  .m <- .rxFoceiCarryMicro(ncmt, oral0, trans)
  if (is.null(.m)) {
    return(NULL)
  }
  .S <- symengine::S
  .x <- lapply(xNames, .S)
  .c <- oral0 + 1L # central index (1-based) in the kernel's row order
  .k <- .m$k
  .k12 <- if (ncmt >= 2L) .m$k12 else .S("0")
  .k21 <- if (ncmt >= 2L) .m$k21 else .S("0")
  .k13 <- if (ncmt >= 3L) .m$k13 else .S("0")
  .k31 <- if (ncmt >= 3L) .m$k31 else .S("0")
  .rows <- vector("list", ncmt + oral0)
  .central <- -(.k + .k12 + .k13) * .x[[.c]]
  if (ncmt >= 2L) .central <- .central + .k21 * .x[[.c + 1L]]
  if (ncmt >= 3L) .central <- .central + .k31 * .x[[.c + 2L]]
  if (oral0 == 1L) {
    .ka <- .S(.rxFoceiCarrySlotPh[7])
    .rows[[1L]] <- -.ka * .x[[1L]]
    .central <- .central + .ka * .x[[1L]]
  }
  .rows[[.c]] <- .central
  if (ncmt >= 2L) .rows[[.c + 1L]] <- .k12 * .x[[.c]] - .k21 * .x[[.c + 1L]]
  if (ncmt >= 3L) .rows[[.c + 2L]] <- .k13 * .x[[.c]] - .k31 * .x[[.c + 2L]]
  vapply(.rows, function(r) {
    .r <- paste(.rxFoceiCarrySubsSlots(r, slotExpr))
    rxode2::rxFromSE(.r)
  }, character(1))
}

#' Model-level event modifiers on the linCmt() compartments: the single
#' modeled f() and/or alag() (the #920 preconditions -- more than one of a
#' kind is unsupported), each as list(cmt, sym, drivers)
#' @noRd
.rxFoceiCarryEventMods <- function(ui, s, etaVars) {
  .lin <- rxode2::.rxLinCmt(ui)
  .lin <- grep("^rx__sens_", .lin, value = TRUE, invert = TRUE)
  if (length(.lin) == 0L) {
    return(list(f = NULL, lag = NULL))
  }
  .one <- function(kind) {
    .rows <- .rxFoceiLinCmtEventRows(s, .lin, kind, etaVars) # nolint: object_usage_linter.
    if (length(.rows) != 1L) {
      return(NULL)
    }
    list(cmt = names(.rows), sym = .rows[[1]]$sym, drivers = .rows[[1]]$drivers)
  }
  list(f = .one("f"), lag = .one("lag"))
}

#' Per-eta jump information: list(fD, fCov, lagD, ok); fD/lagD are symengine
#' reprs (NULL when the eta does not drive that modifier); ok=FALSE means the
#' eta drives a modifier in a way the carry cannot represent (bias to the
#' status quo for that eta)
#' @noRd
.rxFoceiCarryEtaJump <- function(eta, mods, allCovs) {
  .fD <- NULL
  .fCov <- FALSE
  .lagD <- NULL
  if (!is.null(mods$f) && eta %in% mods$f$drivers) {
    .d <- symengine::D(mods$f$sym, symengine::S(eta))
    if (!.rxFoceiCarryIsZero(.d)) { # nolint: object_usage_linter.
      .dln <- .d / mods$f$sym
      .fD <- paste(.dln)
      # every input is proportional to F, so d(pred)/d(eta) = pred * dlnF
      # is exact -- across a covariate change too -- unless dlnF itself
      # changes with the covariate (the non-multiplicative shapes)
      .fCov <- length(intersect(.rxFoceiCarryFreeSyms(.dln), allCovs)) > 0L # nolint: object_usage_linter.
    }
  }
  if (!is.null(mods$lag) && eta %in% mods$lag$drivers) {
    .d <- symengine::D(mods$lag$sym, symengine::S(eta))
    if (!.rxFoceiCarryIsZero(.d)) { # nolint: object_usage_linter.
      # a covariate-driven lag gives every dose its own d(lag)/d(eta); the
      # boundary terms need the entering dose's value, which a later row
      # cannot recover
      if (length(intersect(.rxFoceiCarryFreeSyms(mods$lag$sym), allCovs)) > 0L) { # nolint: object_usage_linter.
        return(list(ok = FALSE))
      }
      .lagD <- paste(.d)
    }
  }
  list(fD = .fD, fCov = .fCov, lagD = .lagD, ok = TRUE)
}

#' Does any linCmt() slot of the model reference a covariate (a time-varying
#' kernel makes the cumulative #920 alag term wrong even for a lag that has
#' no covariate of its own)?
#' @noRd
.rxFoceiCarryKernelHasCov <- function(slotFree, allCovs) {
  any(vapply(slotFree, function(f) length(intersect(f, allCovs)) > 0L, logical(1)))
}
