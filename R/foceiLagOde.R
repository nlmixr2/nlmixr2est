# Lagged calculated variables inside ODEs (#1176): symengine binds one to a
# bare symbol, so its definition is substituted before the sensitivity ODEs
# are derived.

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
    # a snapshot is expanded with the variable it holds a value of
    if (grepl("^rx_lagv[0-9]+_", v) || .foceiLagRefs(.rhs[.w], .var, hist = TRUE)) {
      next
    }
    # a pruned if/else defines it more than once, each from the one before;
    # newer rxode2 reads value i through the snapshot rx_lagv<i>_<v>
    .def <- NULL
    .snap <- list()
    for (.txt in .rhs[.w]) {
      # rxToSE() is NSE: hand it a plain variable
      .new <- symengine::S(rxode2::rxToSE(.txt))
      if (!is.null(.def)) {
        .new <- .foceiLagSubsSnap(symengine::subs(.new, symengine::S(v), .def), v, .snap)
      }
      .def <- .new
      .snap[[length(.snap) + 1L]] <- .def
    }
    for (d in .ddt) {
      .e <- symengine::subs(get(d, envir = s), symengine::S(v), .def)
      assign(d, .foceiLagSubsSnap(.e, v, .snap), envir = s)
    }
  }
  # the ODE text too: the lagged variable is defined after the ODEs
  for (d in .ddt) {
    .pre <- paste0("d/dt(", sub("^rx__d_dt_(.*)__$", "\\1", d), ")=")
    .w <- which(startsWith(s$..ddt, .pre))
    if (length(.w) == 1L && .foceiLagRefs(s$..ddt[.w], c(.var, .foceiLagSnapNames(s$..ddt[.w])))) {
      .e <- get(d, envir = s)
      s$..ddt[.w] <- paste0(.pre, rxode2::rxFromSE(.e))
    }
  }
  invisible(s)
}

#' Substitute the snapshots of a lagged variable's earlier values
#'
#' @param e symengine expression
#' @param v lagged variable name
#' @param snap list of `v`'s values, one per assignment so far
#' @return `e` with each `rx_lagv<i>_<v>` replaced by `snap[[i]]`
#' @noRd
.foceiLagSubsSnap <- function(e, v, snap) {
  for (.i in seq_along(snap)) {
    e <- symengine::subs(e, symengine::S(paste0("rx_lagv", .i, "_", v)), snap[[.i]])
  }
  e
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
  .var <- unique(sub("=.*$", "", .defs))
  .foceiLagRefs(.txt, c(.var, .foceiLagSnapNames(.txt)))
}

#' Snapshot names (`rx_lagv<i>_<v>`) used in rxode2 text
#'
#' @param txt rxode2 text
#' @return character vector of the snapshot names
#' @noRd
.foceiLagSnapNames <- function(txt) {
  unique(unlist(regmatches(txt, gregexpr("\\brx_lagv[0-9]+_[A-Za-z0-9_.]+", txt, perl = TRUE))))
}
