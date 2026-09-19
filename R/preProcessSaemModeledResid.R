#' Normal endpoints with a modeled (non-theta) residual error component
#'
#' rxode2 stores a model variable used as an error argument (`a` in
#' `a <- add.sd*exp(eta.sd); cp ~ add(a)`) by name in the `predDf` argument
#' columns; a plain theta leaves them `NA`.
#'
#' @param ui rxode2 ui
#' @param dist endpoint distributions to consider
#' @return character vector of endpoint conditions
#' @noRd
.saemModeledResidualCond <- function(ui, dist = "norm") {
  .pred <- ui$predDf
  if (is.null(.pred) || length(.pred$cond) == 0L) {
    return(character(0))
  }
  # `f` is propF()/powF()'s prediction variable, always a name, not a residual parameter
  .cols <- intersect(c("a", "b", "c", "d", "e", "lambda"), names(.pred))
  if (length(.cols) == 0L) {
    return(character(0))
  }
  .modeled <- vapply(
    seq_along(.pred$cond),
    function(i) {
      .pred$distribution[i] %in% dist && any(!is.na(unlist(.pred[i, .cols, drop = TRUE])))
    },
    logical(1),
    USE.NAMES = FALSE
  )
  as.character(.pred$cond[.modeled])
}

#' Whether a modeled residual error component depends on an eta
#'
#' Follows the model's assignments back from each error argument variable.
#'
#' @param ui rxode2 ui
#' @param conds endpoint conditions from `.saemModeledResidualCond()`
#' @return logical
#' @noRd
.saemModeledResidHasEta <- function(ui, conds) {
  .pred <- ui$predDf
  .cols <- intersect(c("a", "b", "c", "d", "e", "lambda"), names(.pred))
  .vars <- unique(stats::na.omit(unlist(.pred[.pred$cond %in% conds, .cols, drop = FALSE])))
  .iniDf <- ui$iniDf
  .etas <- .iniDf$name[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2]
  .rhs <- list()
  for (.e in ui$lstExpr) {
    if (is.call(.e) && (identical(.e[[1]], as.name("<-")) || identical(.e[[1]], as.name("="))) && is.name(.e[[2]])) {
      .n <- as.character(.e[[2]])
      .rhs[[.n]] <- c(.rhs[[.n]], all.vars(.e[[3]]))
    }
  }
  .seen <- character(0)
  while (length(.vars) > 0L) {
    if (any(.vars %in% .etas)) {
      return(TRUE)
    }
    .seen <- c(.seen, .vars)
    .vars <- setdiff(unique(unlist(.rhs[intersect(.vars, names(.rhs))])), .seen)
  }
  FALSE
}

#' Raise an unset saem `nu` for a residual error with an eta
#'
#' @param control saem control as passed to the hook (may be `NULL`)
#' @return the updated control, or `NULL` when `nu` was set explicitly
#' @noRd
.saemAutoNu <- function(control) {
  if (is.null(control)) {
    control <- saemControl()
  }
  if (!inherits(control, "saemControl") || !isTRUE(control$mcmc$nuAuto)) {
    return(NULL)
  }
  control$mcmc$nu <- pmax(control$mcmc$nu, 4)
  control
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
#' equivalent `cp ~ add(a) + dnorm()`.  Any eta-less likelihood theta of a
#' general-likelihood fit then gets a temporary eta (`.saemAddPseudoEtas()`).
#'
#' @param ui rxode2 ui
#' @param est estimation method
#' @param data dataset (unused)
#' @param control saem control; an unset `nu` is raised when the residual
#'   error depends on an eta (`.saemAutoNu()`)
#' @return list with the rewritten ui (and control), or `NULL` when nothing changes
#' @noRd
.preProcessSaemModeledResid <- function(ui, est, data, control) {
  nlmixr2global$nlmixr2EstEnv$saemPseudoTransforms <- NULL
  if (!identical(est, "saem")) {
    return(NULL)
  }
  .orig <- ui
  .conds <- .saemModeledResidualCond(ui)
  .pred <- ui$predDf
  for (.cond in .conds) {
    .new <- .saemAddDnormToErrLine(ui$lstExpr[[.pred$line[.pred$cond == .cond]]])
    ui <- eval(bquote(rxode2::model(ui, .(.new))))
    warning(sprintf("modeled residual error for '%s'; fit as dnorm() likelihood", .cond), call. = FALSE)
  }
  # an explicit + dnorm() is the same likelihood, so it gets the same nu
  .ctl <- NULL
  .etaConds <- .saemModeledResidualCond(.orig, c("norm", "dnorm"))
  if (length(.etaConds) > 0L && .saemModeledResidHasEta(.orig, .etaConds)) {
    .ctl <- .saemAutoNu(control)
    if (!is.null(.ctl)) {
      warning("residual error depends on an eta; MCMC nu raised to ", deparse1(.ctl$mcmc$nu), call. = FALSE)
    }
  }
  .spec <- .saemPseudoEtaThetas(ui)
  if (length(.conds) == 0L && nrow(.spec) == 0L) {
    if (is.null(.ctl)) {
      return(NULL)
    }
    return(list(control = .ctl))
  }
  if (nrow(.spec) > 0L) {
    ui <- .saemAddPseudoEtas(ui, .spec)
    # a later hook can rebuild the ui (IOV) and drop these; saem puts them back
    nlmixr2global$nlmixr2EstEnv$saemPseudoTransforms <- ui$boundedTransforms
    .pre <- "temporary eta for eta-less likelihood theta(s): "
    warning(.pre, .vaeTruncList(.spec$theta, prefix = .pre), call. = FALSE)
  }
  # the reported fit shows the model as written (see .nlmixrEstUpdatesOrigModel)
  if (is.null(nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi)) {
    nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- rxode2::rxUiDecompress(.orig)
  }
  .ret <- list(ui = ui)
  if (!is.null(.ctl)) {
    .ret$control <- .ctl
  }
  .ret
}

preProcessHooksAdd(".preProcessSaemModeledResid", .preProcessSaemModeledResid)
