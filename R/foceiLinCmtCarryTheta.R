# Theta-side linCmt() sensitivity carry for the nlm family and nls (#1003).
#
# The carry (foceiLinCmtCarry.R) is direction-agnostic: a pair is a linCmt()
# slot (and/or event channel) times a direction symbol, and rxode2's
# which1=-5/-6/-7/-8 sentinels never look at what the direction is.  The
# FOCEi-family hook drives it with ETA_n_ directions; the population
# methods that request d(pred)/d(THETA_k_) symbolically (rxUiGet.nlmHdTheta,
# rxUiGet.nlsHdTheta) suffer the same conflation for a theta on a
# covariate-driven parameter and drive it here with THETA_k_ directions.
#
# Unlike the FOCEi inner model, these methods' rx_pred_ is not a bare
# linCmtB() call: the nlm family wraps the concentration in a
# log-likelihood (-llikNorm(DV, linCmtB(...), sd)).  The single linCmtB()
# call is therefore factored out as its own intermediate, rx_lcConc_, and
#
#   d(rx_pred_)/d(THETA_k_) = d(outer)/d(rx_lcConc_) * d(rx_lcConc_)/d(THETA_k_)
#                             + d(outer)/d(THETA_k_)|rx_lcConc_ fixed
#
# where outer is rx_pred_ with the call replaced by the symbol, the first
# factor is an ordinary symengine partial, the carried concentration
# sensitivity comes from the unchanged carry emission, and the last term is
# the direct dependence (a theta that also sits in the residual sd).
#
# Out of scope (documented): the FOCEi outer theta gradient (the analytic
# outer path already downgrades for linCmt()), and SAEM -- its linearized
# FIM (calc.COV) perturbs phi and re-solves the values, which is exact
# under a time-varying covariate, and saem.cpp consumes no symbolic
# sensitivity.  All intermediates are `~` so the population methods'
# positional lhs layout (rx_pred_ then the theta columns) is unchanged.

#' Every linCmtB() call inside a symengine expression
#' @noRd
.rxCarryFindLinCmtB <- function(expr, acc = list()) {
  .nm <- tryCatch(symengine::get_name(expr), error = function(e) "")
  if (identical(.nm, "linCmtB")) return(c(acc, list(expr)))
  .a <- tryCatch(symengine::get_args(expr), error = function(e) NULL)
  if (!is.null(.a)) {
    for (.i in seq_along(.a)) acc <- .rxCarryFindLinCmtB(.a[[.i]], acc)
  }
  acc
}

#' Factor the single 15-argument linCmtB() call out of rx_pred_
#'
#' @return list(call, args, outer) -- the call, its arguments, and rx_pred_
#'   with the call replaced by the symbol rx_lcConc_ -- or NULL when
#'   rx_pred_ does not contain exactly one structural call
#' @noRd
.rxCarryFactorPred <- function(s) {
  if (!exists("rx_pred_", envir = s, inherits = FALSE)) return(NULL)
  .pred <- get("rx_pred_", envir = s, inherits = FALSE)
  .calls <- .rxCarryFindLinCmtB(.pred)
  if (length(.calls) != 1L) return(NULL)
  .args <- symengine::get_args(.calls[[1]])
  if (length(.args) != 15L) return(NULL)
  # the structural (value) call has which1 = which2 = -1
  .w <- suppressWarnings(as.numeric(c(paste(.args[[6]]), paste(.args[[7]]))))
  if (any(is.na(.w)) || any(.w != -1)) return(NULL)
  .outer <- symengine::subs(.pred, .calls[[1]], symengine::S("rx_lcConc_"))
  list(call = .calls[[1]], args = .args, outer = .outer)
}

#' The iniDf name of THETA_k_
#' @noRd
.rxCarryThetaName <- function(ui, k) {
  .iniDf <- ui$iniDf
  .w <- which(!is.na(.iniDf$ntheta) & .iniDf$ntheta == k)
  if (length(.w) == 1L) .iniDf$name[.w] else paste0("THETA_", k, "_")
}

#' Carry-eligible (theta, channels) pairs for the population methods
#'
#' Same rules as `.rxFoceiLinCmtCarryEligible()` with THETA_k_ as the
#' direction symbols (a slot must hold exactly one estimated theta and enter
#' it separably; event channels as for etas); the IOV check does not apply.
#' @noRd
.rxCarryThetaEligible <- function(x, s, thetaVars, fp, render = FALSE) {
  .ui <- x[[1]]
  .empty <- .rxFoceiCarryEmpty() # nolint: object_usage_linter.
  .allCovs <- .ui$allCovs
  if (length(.allCovs) == 0L) return(.empty)
  .slotExpr <- lapply(1:7, function(k) fp$args[[k + 8L]])
  .slotFree <- lapply(.slotExpr, .rxFoceiCarryFreeSyms) # nolint: object_usage_linter.
  .mods <- .rxFoceiCarryEventMods(.ui, s, thetaVars) # nolint: object_usage_linter.
  .ret <- .empty
  for (.k in seq_along(thetaVars)) {
    .th <- thetaVars[.k]
    .jump <- .rxFoceiCarryEtaJump(.th, .mods, .allCovs) # nolint: object_usage_linter.
    if (!isTRUE(.jump$ok)) next
    .slot <- .rxFoceiCarryEtaSlotInfo(.th, thetaVars, .allCovs, .slotExpr, .slotFree) # nolint: object_usage_linter.
    if (is.null(.slot)) next
    .hasSlot <- length(.slot) > 0L
    .why <- .rxFoceiCarryWhy(if (.hasSlot) .slot else NULL, .jump, .mods, .slotFree, .allCovs) # nolint: object_usage_linter.
    if (length(.why) == 0L) next
    .ret <- rbind(.ret, .rxFoceiCarryPairRow(.th, .rxCarryThetaName(.ui, .k), # nolint: object_usage_linter.
                                             if (.hasSlot) .slot else NULL,
                                             .jump, .mods, .why, NA, render))
  }
  .ret
}

#' ncmt / oral0 / trans of a factored linCmtB() call
#' @noRd
.rxCarryCallShape <- function(fp) {
  .num <- function(i) suppressWarnings(as.numeric(paste(fp$args[[i]])))
  list(ncmt = as.integer(.num(4)), oral0 = as.integer(.num(5)), trans = .num(8))
}

#' Theta pairs gated for substitution in a population method's gradient
#'
#' @param x list(rxUi)
#' @param s symengine env after `.sensEtaOrTheta(theta = TRUE)`
#' @param thetaVars THETA_k_ names in gradient-column order
#' @return NULL, or list(pairs, fp) for `.rxCarryThetaEmit()`
#' @noRd
.rxCarryThetaPairsForBuild <- function(x, s, thetaVars) {
  .fp <- .rxCarryThetaModelOk(x[[1]], s)
  if (is.null(.fp)) return(NULL)
  .pairs <- .rxCarryThetaDetect(x, s, thetaVars, .fp)
  if (is.null(.pairs)) return(NULL)
  .rxFoceiLinCmtCarryCheckCap(.pairs) # nolint: object_usage_linter.
  list(pairs = .pairs, fp = .fp)
}

#' The factored call when the model/control can carry at all, else NULL
#' @noRd
.rxCarryThetaModelOk <- function(ui, s) {
  if (!.rxFoceiLinCmtCarryBuildEnabled(ui)) return(NULL) # nolint: object_usage_linter.
  if (rxode2::.rxLinNcmt(ui)["numLin"] <= 0L) return(NULL)
  # the once-per-row stepping is only validated for pure linCmt() models
  if (length(.rxode2stateOdeNoOutput(s)) > 0L) return(NULL)
  .fp <- .rxCarryFactorPred(s)
  if (is.null(.fp)) return(NULL)
  .sh <- .rxCarryCallShape(.fp)
  if (any(is.na(unlist(.sh))) ||
        is.null(.rxFoceiCarryMicro(.sh$ncmt, .sh$oral0, .sh$trans))) { # nolint: object_usage_linter.
    return(NULL)
  }
  .fp
}

#' Detected theta pairs, slot channels only; NULL when none qualify
#' @noRd
.rxCarryThetaDetect <- function(x, s, thetaVars, fp) {
  .pairs <- tryCatch(.rxCarryThetaEligible(x, s, thetaVars, fp, render = FALSE),
                     error = function(e) {
                       warning("linCmt() theta carry detection failed; standard gradient used",
                               call. = FALSE)
                       NULL
                     })
  if (is.null(.pairs) || nrow(.pairs) == 0L) return(NULL)
  # slot channels only: a theta on f()/alag() is a dosing parameter the nlm
  # family already routes through its own eventSens machinery
  # (rxUiGet.nlmEnv's eventTheta), and that interaction is unvalidated
  .pairs <- .pairs[is.na(.pairs$fD) & is.na(.pairs$lagD), , drop = FALSE]
  if (nrow(.pairs) == 0L) return(NULL)
  .pairs
}

#' Substituted d(rx_pred_)/d(THETA_k_) lines for one theta pair
#'
#' Runs the unchanged carry emission against the factored call (so it
#' produces the carried concentration sensitivity into an intermediate),
#' then composes the outer chain rule.  When `predMinusDv` is FALSE the
#' consumer expects -d(pred)/d(theta) (rxExpandFEta_'s pred = 2 form).
#' @noRd
.rxCarryThetaEmit <- function(pairs, w, s, dfe, fp, predMinusDv = TRUE) {
  .p <- w - 1L
  .proxy <- new.env(parent = emptyenv())
  assign("rx_pred_", fp$call, envir = .proxy)
  .tmp <- paste0("rx_lcCarryC", .p, "_")
  .lines <- strsplit(.rxFoceiLinCmtCarryEmit(pairs, w, .proxy, .tmp), "\n", fixed = TRUE)[[1]] # nolint: object_usage_linter.
  .i <- length(.lines)
  # the emission's final line assigns the concentration sensitivity to
  # `tmp=`; make it an intermediate and point its direct term at the
  # factored concentration rather than at rx_pred_
  .lines[.i] <- gsub("rx_pred_", "rx_lcConc_",
                     sub(paste0("^", .tmp, "="), paste0(.tmp, "~"), .lines[.i]),
                     fixed = TRUE)
  if (w == 1L) {
    .callRepr <- paste(fp$call)
    .lines <- c(paste0("rx_lcConc_~", rxode2::rxFromSE(.callRepr)), .lines)
  }
  .S <- symengine::S
  .dcRepr <- paste(symengine::D(fp$outer, .S("rx_lcConc_")))
  .out <- paste0("(", rxode2::rxFromSE(.dcRepr), ")*", .tmp)
  .dt <- symengine::D(fp$outer, .S(pairs$eta[w]))
  if (!.rxFoceiCarryIsZero(.dt)) { # nolint: object_usage_linter.
    .dtRepr <- paste(.dt)
    .out <- paste0(.out, "+(", rxode2::rxFromSE(.dtRepr), ")")
  }
  if (!isTRUE(predMinusDv)) .out <- paste0("-(", .out, ")")
  paste(c(.lines, paste0(dfe, "=", .out)), collapse = "\n")
}
