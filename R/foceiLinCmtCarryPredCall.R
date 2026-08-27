# Locating the structural linCmtB() value call of rx_pred_ for the linCmt()
# sensitivity carry.  rx_pred_ is either the call itself (the add/prop
# endpoints) or a generalized ll() likelihood embedding it (#1004); in both
# cases the concentration the carry differentiates is that single call, and
# a wrapped prediction is factored through the symbol rx_lcConc_ so the
# emission (foceiLinCmtCarryCompose.R) can hand symengine the outer chain
# rule.

#' Every `linCmtB()` call inside a symengine expression (depth-first)
#' @noRd
.rxFoceiCarryFindLinCmtB <- function(expr) {
  .calls <- list()
  .walk <- function(x) {
    .nm <- tryCatch(symengine::get_name(x), error = function(e) NULL)
    if (identical(.nm, "linCmtB")) {
      .calls[[length(.calls) + 1L]] <<- x
      return(invisible())
    }
    .a <- tryCatch(symengine::get_args(x), error = function(e) NULL)
    if (is.null(.a)) {
      return(invisible())
    }
    for (.i in seq_along(.a)) .walk(.a[[.i]])
    invisible()
  }
  .walk(expr)
  .calls
}

#' The structural `linCmtB()` value call of `rx_pred_`
#'
#' `rx_pred_` is either the call itself (add/prop endpoints) or an ll()
#' expression embedding it; in both cases the concentration the carry
#' differentiates is that single call.  NULL when there is no call, more
#' than one distinct call, or the call is not the 15-argument value form.
#'
#' @return list(call = symengine FunctionSymbol, args = its arguments,
#'   bare = is `rx_pred_` the call itself, predSym = `rx_pred_` with the
#'   call replaced by the symbol `rx_lcConc_`)
#' @noRd
.rxFoceiCarryPredCall <- function(s) {
  if (!exists("rx_pred_", envir = s, inherits = FALSE)) {
    return(NULL)
  }
  .pred <- get("rx_pred_", envir = s, inherits = FALSE)
  .calls <- .rxFoceiCarryFindLinCmtB(.pred)
  if (length(.calls) == 0L) {
    return(NULL)
  }
  .reprs <- vapply(.calls, paste, character(1))
  if (length(unique(.reprs)) != 1L) {
    return(NULL)
  }
  .call <- .calls[[1L]]
  .args <- symengine::get_args(.call)
  if (length(.args) != 15L) {
    return(NULL)
  }
  # the value form (which1 = -1); a sensitivity read inside rx_pred_ is not
  # something the carry can stand in for
  if (!isTRUE(all.equal(suppressWarnings(as.numeric(paste(.args[[6]]))), -1))) {
    return(NULL)
  }
  .bare <- identical(paste(.pred), .reprs[1L])
  list(
    call = .call, args = .args, bare = .bare,
    predSym = if (.bare) {
      symengine::S("rx_lcConc_")
    } else {
      symengine::subs(.pred, .call, symengine::S("rx_lcConc_"))
    }
  )
}

#' The 15 `linCmtB()` arguments of `rx_pred_`'s structural call, or NULL
#' when the model is not a pure linCmt() model with one such call
#' @noRd
.rxFoceiCarryPredArgs <- function(ui, s) {
  if (rxode2::.rxLinNcmt(ui)["numLin"] <= 0L) {
    return(NULL)
  }
  # mixed ODE + linCmt() models evaluate linCmtB() inside the integrator as
  # well as in the lhs pass; the once-per-row carry stepping is only
  # validated for pure linCmt() models
  if (length(.rxode2stateOdeNoOutput(s)) > 0L) {
    return(NULL)
  }
  .pc <- .rxFoceiCarryPredCall(s)
  if (is.null(.pc)) {
    return(NULL)
  }
  .pc$args
}
