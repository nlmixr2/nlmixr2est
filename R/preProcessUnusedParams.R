#' Variables the objective function depends on
#'
#' Walks the normalized model (with `linCmt()` expanded) back from what feeds
#' the likelihood: the endpoint lines and any assignment to a call other than a
#' state's (`mtime()`, ...).  A state's ODE and properties (`f()`, `alag()`,
#' `x(0)`, ...) feed that state, the `linCmtA()`/`linCmtB()` arguments feed
#' the `linCmt()` states, and an assignment inside `if()`/`while()` also depends on its
#' condition.  The walk
#' ignores statement order, so it may keep a variable whose value is later
#' overwritten (never the reverse).
#'
#' @param ui rxode2 ui
#' @return character vector of variable names, or `NULL` when the model cannot
#'   be walked
#' @noRd
.objectiveDepends <- function(ui) {
  tryCatch(
    {
      .norm <- rxode2::rxGetLin(rxode2::rxNorm(ui$mv0))
      # braces let a normalized `else` start its own line
      .exprs <- str2lang(paste0("{\n", .norm, "\n}"))
      .states <- ui$mv0$state
      .linStates <- intersect(.states, c("depot", "central", "peripheral1", "peripheral2"))
      .acc <- new.env(parent = emptyenv())
      .acc$lhs <- character(0)
      .acc$rhs <- list()
      .acc$seed <- character(0)
      .linVars <- function(e) {
        if (!is.call(e)) {
          return(character(0))
        }
        if (identical(e[[1]], quote(linCmtA)) || identical(e[[1]], quote(linCmtB))) {
          return(all.vars(e))
        }
        unlist(lapply(as.list(e)[-1], .linVars))
      }
      .walk <- function(e, cond) {
        if (!is.call(e)) {
          return(invisible())
        }
        .f <- e[[1]]
        if (identical(.f, quote(`{`))) {
          for (.a in as.list(e)[-1]) {
            .walk(.a, cond)
          }
        } else if (identical(.f, quote(`if`)) || identical(.f, quote(`while`))) {
          .c <- c(cond, all.vars(e[[2]]))
          for (.a in as.list(e)[-(1:2)]) {
            .walk(.a, .c)
          }
        } else if (
          (identical(.f, quote(`=`)) || identical(.f, quote(`<-`)) || identical(.f, quote(`~`))) &&
            length(e) == 3L
        ) {
          .l <- e[[2]]
          .r <- e[[3]]
          .lin <- .linVars(.r)
          for (.s in .linStates[length(.lin) > 0L]) {
            .acc$lhs <- c(.acc$lhs, .s)
            .acc$rhs[[length(.acc$lhs)]] <- unique(c(cond, .lin))
          }
          .t <- .lhsTarget(.l, .states)
          if (is.na(.t)) {
            # mtime() and the like
            .acc$seed <- c(.acc$seed, cond, all.vars(e))
          } else {
            .acc$lhs <- c(.acc$lhs, .t)
            .acc$rhs[[length(.acc$lhs)]] <- unique(c(cond, all.vars(.r)))
          }
        }
        invisible()
      }
      .walk(.exprs, character(0))
      .lst <- ui$lstExpr
      for (.l in unique(ui$predDf$line)) {
        .acc$seed <- c(.acc$seed, all.vars(.lst[[.l]]))
      }
      .acc$seed <- c(.acc$seed, as.character(ui$predDf$var))
      .need <- unique(.acc$seed)
      repeat {
        .new <- setdiff(unique(unlist(.acc$rhs[.acc$lhs %in% .need])), .need)
        if (length(.new) == 0L) {
          break
        }
        .need <- c(.need, .new)
      }
      .need
    },
    error = function(e) NULL
  )
}

#' The variable an assignment defines
#'
#' @param l left hand side of an assignment
#' @param states model states
#' @return the assigned name, the state for `d/dt(x)`, `f(x)`, `alag(x)`,
#'   `x(0)` and the like, or `NA` for any other call
#' @noRd
.lhsTarget <- function(l, states) {
  if (is.name(l)) {
    return(as.character(l))
  }
  if (is.call(l) && identical(l[[1]], quote(`/`)) && identical(l[[2]], quote(d))) {
    l <- l[[3]]
  }
  if (is.call(l) && length(l) == 2L) {
    for (.x in list(l[[2]], l[[1]])) {
      if (is.name(.x) && as.character(.x) %in% states) {
        return(as.character(.x))
      }
    }
  }
  NA_character_
}

#' Estimated parameters that cannot change the objective function
#'
#' @param ui rxode2 ui
#' @param control estimation control; `priorMethod = "none"` ignores priors
#' @return names of non-fixed thetas and etas that the objective does not
#'   depend on
#' @noRd
.unusedEstimatedParams <- function(ui, control = NULL) {
  .need <- .objectiveDepends(ui)
  if (is.null(.need)) {
    return(character(0))
  }
  .iniDf <- ui$iniDf
  # residual parameters live on endpoint lines; off-diagonals have no symbol
  .w <- !.iniDf$fix &
    is.na(.iniDf$err) &
    (!is.na(.iniDf$ntheta) | (!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2))
  # a prior adds its own term to the objective
  .prior <- if (identical(control$priorMethod, "none")) {
    NULL
  } else {
    tryCatch(rxode2::rxUiPriors(ui)$name, error = function(e) NULL)
  }
  setdiff(.iniDf$name[.w], c(.need, .prior))
}

#' Warn about estimated parameters that do not reach the objective
#'
#' A parameter used only through an intermediate that nothing else reads (a
#' leftover from editing the model) has a zero gradient, which can stall the
#' optimizer with no hint of the cause.
#'
#' @param ui rxode2 ui model
#' @inheritParams nlmixr2
#' @return `NULL` (this hook only warns)
#' @noRd
.preProcessUnusedParams <- function(ui, est, data, control) {
  if (isTRUE(est %in% c(.preProcessHooksSimulationEst, "posthoc"))) {
    return(NULL)
  }
  # a registered likelihood contributor may read any parameter
  if (.Call("_nlmixr2est_likContribN", PACKAGE = "nlmixr2est") > 0L) {
    return(NULL)
  }
  .ui <- rxode2::rxUiDecompress(ui)
  for (.n in .unusedEstimatedParams(.ui, control)) {
    warning(paste0("'", .n, "' does not affect the objective function"), call. = FALSE)
  }
  NULL
}
preProcessHooksAdd(".preProcessUnusedParams", .preProcessUnusedParams)
