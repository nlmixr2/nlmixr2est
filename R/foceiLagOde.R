# Lagged calculated variables inside ODEs (#1176): symengine binds one to a
# bare symbol, so its definition is substituted before the sensitivity ODEs
# are derived.

#' Replace rxode2's snapshots of a reassigned lagged variable by their values
#'
#' rxode2 reads an earlier value of a lagged variable assigned more than once
#' (a pruned if/else) through a snapshot symbol `rx_lagv<i>_<v>`, defined by
#' an `rx_lagv<i>_<v>=<v>` lhs.  symengine cannot differentiate through that
#' symbol, so each snapshot is replaced by the value it holds.
#' @param s symengine environment from `rxode2::rxS()`
#' @return `s`, invisibly
#' @noRd
.foceiLagSnapResolve <- function(s) {
  .lhs <- s$..lhs
  .isSnap <- grepl("^rx_lagv[0-9]+_[^=]+=", .lhs)
  if (!any(.isSnap)) {
    return(invisible(s))
  }
  .subTxt <- function(x, nm, val) {
    for (j in seq_along(nm)) {
      x <- gsub(paste0("(?<![A-Za-z0-9_.])", nm[j], "(?![A-Za-z0-9_.])"), val[j], x, perl = TRUE)
    }
    x
  }
  .nm <- .val <- character(0)
  for (i in which(.isSnap)) {
    .v <- sub("^[^=]*=", "", .lhs[i])
    .prev <- which(startsWith(.lhs[seq_len(i - 1L)], paste0(.v, "=")))
    .nm <- c(.nm, sub("=.*$", "", .lhs[i]))
    .val <- c(.val, paste0("(", .subTxt(sub("^[^=]*=", "", .lhs[max(.prev)]), .nm, .val), ")"))
  }
  .se <- lapply(.val, function(x) symengine::S(rxode2::rxToSE(x)))
  rm(list = intersect(.nm, ls(s, all.names = TRUE)), envir = s)
  for (n in ls(s, all.names = TRUE)) {
    .e <- get(n, envir = s)
    if (inherits(.e, "Basic")) {
      for (j in seq_along(.nm)) {
        .e <- symengine::subs(.e, symengine::S(.nm[j]), .se[[j]])
      }
      assign(n, .e, envir = s)
    } else if (startsWith(n, "..") && is.character(.e) && any(grepl("rx_lagv", .e, fixed = TRUE))) {
      assign(n, .subTxt(.e, .nm, .val), envir = s)
    }
  }
  s$..lhs <- .subTxt(.lhs[!.isSnap], .nm, .val)
  invisible(s)
}

#' Put the definitions of lagged variables into the ODEs that use them
#'
#' Run before the sensitivity ODEs are derived, so they do not treat a lagged
#' variable as a constant.  A variable defined through a history call is left
#' alone (see `.foceiLagInOde()`).
#' @param s symengine environment, before `.sensEtaOrTheta()`
#' @return `s`, invisibly
#' @noRd
.foceiLagIntoOde <- function(s) {
  .defs <- .foceiLagDefs(s)
  if (length(.defs) == 0L) {
    return(invisible(s))
  }
  .var <- sub("=.*$", "", .defs)
  .rhs <- sub("^[^=]*=", "", .defs)
  .ddt <- ls(s, pattern = "^rx__d_dt_.*__$", all.names = TRUE)
  .ddt <- .ddt[
    !vapply(
      .ddt,
      function(d) {
        # rxFromSE() is NSE: hand it a plain variable
        .e <- get(d, envir = s)
        .foceiLagRefs(rxode2::rxFromSE(.e), .var, hist = TRUE)
      },
      logical(1)
    )
  ]
  # last definition first, so one that uses an earlier one is fully expanded
  for (v in rev(unique(.var))) {
    .w <- which(.var == v)
    if (.foceiLagRefs(.rhs[.w], .var, hist = TRUE)) {
      next
    }
    # a pruned if/else defines it more than once, each from the one before
    .def <- NULL
    for (.txt in .rhs[.w]) {
      # rxToSE() is NSE: hand it a plain variable
      .new <- symengine::S(rxode2::rxToSE(.txt))
      .def <- if (is.null(.def)) .new else symengine::subs(.new, symengine::S(v), .def)
    }
    for (d in .ddt) {
      assign(d, symengine::subs(get(d, envir = s), symengine::S(v), .def), envir = s)
    }
  }
  # the ODE text too: the lagged variable is defined after the ODEs
  for (d in .ddt) {
    .pre <- paste0("d/dt(", sub("^rx__d_dt_(.*)__$", "\\1", d), ")=")
    .w <- which(startsWith(s$..ddt, .pre))
    if (length(.w) == 1L && .foceiLagRefs(s$..ddt[.w], .var)) {
      .e <- get(d, envir = s)
      s$..ddt[.w] <- paste0(.pre, rxode2::rxFromSE(.e))
    }
  }
  invisible(s)
}

#' Whether an ODE still uses a lagged calculated variable
#'
#' @param s symengine environment, after `.foceiLagIntoOde()`
#' @return `TRUE` when the sensitivity ODEs would treat one as a constant
#' @noRd
.foceiLagInOde <- function(s) {
  .defs <- .foceiLagDefs(s)
  .ddt <- ls(s, pattern = "^rx__d_dt_.*__$", all.names = TRUE)
  if (length(.defs) == 0L || length(.ddt) == 0L) {
    return(FALSE)
  }
  .txt <- vapply(
    .ddt,
    function(d) {
      .e <- get(d, envir = s)
      rxode2::rxFromSE(.e)
    },
    character(1)
  )
  .foceiLagRefs(.txt, unique(sub("=.*$", "", .defs)))
}
