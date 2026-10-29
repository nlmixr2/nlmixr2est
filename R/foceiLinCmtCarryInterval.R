# The covariates the linCmt() carry advances an interval with.  rxode2#1483
# advances the linear compartments over (start, end] with that interval's
# covariates (locf the start record's, midpoint their mean); before it, every
# interpolation used the end record's.  The carry lines run at the end record,
# so on an rxode2 with #1483 the -5 advance, the kernel term's g and the
# lag-time kernel read each covariate at its interval value instead.

.carryIntervalEnv <- new.env(parent = emptyenv())

#' Does the loaded rxode2 advance a linCmt() interval with that interval's
#' covariates (rxode2#1483)?  Probed once per session.
#' @noRd
.rxFoceiLinCmtIntervalCovs <- function() {
  if (is.null(.carryIntervalEnv$covs)) {
    .carryIntervalEnv$covs <- tryCatch(
      {
        .m <- suppressMessages(rxode2::rxode2("cl <- rxLcProbeCov\nv <- 1\ncp <- linCmt()"))
        .ev <- data.frame(
          id = 1L,
          time = c(0, 1),
          amt = c(1, 0),
          evid = c(1L, 0L),
          cmt = 1L,
          rxLcProbeCov = c(1, 2)
        )
        .r <- rxode2::rxSolve(.m, events = .ev, covsInterpolation = "locf", returnType = "data.frame")
        # the interval (0, 1] runs at k = 1 with #1483, at the end record's k = 2 before
        isTRUE(abs(.r$cp[1] - exp(-1)) < 1e-8)
      },
      error = function(e) FALSE
    )
  }
  .carryIntervalEnv$covs
}

#' How the carry reads an interval's covariates: "end" (the end record's),
#' "start" (the previous record's) or "mean" (the mean of the two)
#' @noRd
.rxFoceiLinCmtCarryCovMode <- function(ui) {
  .interp <- rxode2::rxGetControl(ui, "rxControl", NULL)$covsInterpolation
  if (is.null(.interp)) {
    .interp <- "locf"
  } else if (is.numeric(.interp)) {
    .interp <- c("linear", "locf", "nocb", "midpoint")[as.integer(.interp) + 1L]
  }
  if (!(.interp %in% c("locf", "midpoint")) || !.rxFoceiLinCmtIntervalCovs()) {
    return("end")
  }
  if (identical(.interp, "locf")) "start" else "mean"
}

#' Tag carry pairs with the interval covariate mode and the model covariates
#' @noRd
.rxFoceiLinCmtCarryTagCovs <- function(pairs, ui) {
  if (is.null(pairs)) {
    return(NULL)
  }
  attr(pairs, "carryCovMode") <- .rxFoceiLinCmtCarryCovMode(ui)
  attr(pairs, "carryCovs") <- ui$allCovs
  pairs
}

#' The interval-covariate lines and a substitution putting them into a
#' symengine expression (the identity in "end" mode)
#' @param pairs tagged pairs
#' @param exprs symengine expressions the substitution will be applied to
#' @return list(lines, sub)
#' @noRd
.rxFoceiLinCmtCarryCovInterval <- function(pairs, exprs) {
  .mode <- attr(pairs, "carryCovMode")
  .covs <- attr(pairs, "carryCovs")
  .id <- list(lines = character(0), sub = function(e) e)
  if (is.null(.mode) || identical(.mode, "end") || length(.covs) == 0L) {
    return(.id)
  }
  .used <- unique(unlist(lapply(exprs, function(e) {
    vapply(as.list(symengine::free_symbols(e)), paste, character(1))
  })))
  .covs <- intersect(.covs, .used)
  if (length(.covs) == 0L) {
    return(.id)
  }
  .sym <- paste0("rx_lcCovI", seq_along(.covs), "_")
  .lines <- unlist(lapply(seq_along(.covs), function(i) {
    .c <- .covs[i]
    .s <- .sym[i]
    # no previous record on a subject's first row, where the advance is over dt = 0
    .l <- c(paste0(.s, "~lag(", .c, ")"), paste0(.s, "~ifelse(is.na(", .s, "),", .c, ",", .s, ")"))
    if (identical(.mode, "mean")) {
      .l <- c(.l, paste0(.s, "~(", .s, "+", .c, ")/2"))
    }
    .l
  }))
  list(
    lines = .lines,
    sub = function(e) {
      for (i in seq_along(.covs)) {
        e <- symengine::subs(e, symengine::S(.covs[i]), symengine::S(.sym[i]))
      }
      e
    }
  )
}
