#' Normal endpoints with a modeled (non-theta) residual error component
#'
#' rxode2 stores a model variable used as an error argument (`a` in
#' `a <- add.sd*exp(eta.sd); cp ~ add(a)`) by name in the `predDf` argument
#' columns; a plain theta leaves them `NA`.
#'
#' @param ui rxode2 ui
#' @return character vector of endpoint conditions
#' @noRd
.saemModeledResidualCond <- function(ui) {
  .pred <- ui$predDf
  if (is.null(.pred) || length(.pred$cond) == 0L) return(character(0))
  .cols <- intersect(c("a", "b", "c", "d", "e", "f", "lambda"), names(.pred))
  if (length(.cols) == 0L) return(character(0))
  .modeled <- vapply(seq_along(.pred$cond), function(i) {
    .pred$distribution[i] == "norm" &&
      any(!is.na(unlist(.pred[i, .cols, drop = TRUE])))
  }, logical(1), USE.NAMES = FALSE)
  as.character(.pred$cond[.modeled])
}

#' Append `+ dnorm()` to an error line, before any `| condition`
#'
#' @param line error model line, e.g. `cp ~ add(a) | cp`
#' @return the line with `dnorm()` added
#' @noRd
.saemAddDnormToErrLine <- function(line) {
  .rhs <- line[[3]]
  if (is.call(.rhs) && identical(.rhs[[1]], as.name("|"))) {
    .rhs[[2]] <- call("+", .rhs[[2]], quote(dnorm()))
  } else {
    .rhs <- call("+", .rhs, quote(dnorm()))
  }
  line[[3]] <- .rhs
  line
}

#' Fit a modeled residual error component as its `dnorm()` likelihood in saem
#'
#' The closed-form residual M-step can only estimate a constant residual
#' parameter, so `cp ~ add(a)` with a modeled `a` is rewritten to the
#' equivalent `cp ~ add(a) + dnorm()`.
#'
#' @param ui rxode2 ui
#' @param est estimation method
#' @param data dataset (unused)
#' @param control control (unused)
#' @return list with the rewritten ui, or `NULL` when nothing changes
#' @noRd
.preProcessSaemModeledResid <- function(ui, est, data, control) {
  if (!identical(est, "saem")) return(NULL)
  .conds <- .saemModeledResidualCond(ui)
  if (length(.conds) == 0L) return(NULL)
  .pred <- ui$predDf
  # dnorm() scores the transformed DV without the lambda-dependent Jacobian
  .lam <- .pred$cond %in% .conds & grepl("boxCox|yeoJohnson", paste(.pred$transform))
  if (any(.lam)) {
    stop("saem cannot fit a modeled residual error with a boxCox()/yeoJohnson() transform ('",
         paste(.pred$cond[.lam], collapse = "', '"), "')", call. = FALSE)
  }
  for (.cond in .conds) {
    .new <- .saemAddDnormToErrLine(ui$lstExpr[[.pred$line[.pred$cond == .cond]]])
    ui <- eval(bquote(rxode2::model(ui, .(.new))))
    warning(sprintf("modeled residual error for '%s'; fit as dnorm() likelihood", .cond),
            call. = FALSE)
  }
  list(ui = ui)
}

preProcessHooksAdd(".preProcessSaemModeledResid", .preProcessSaemModeledResid)
