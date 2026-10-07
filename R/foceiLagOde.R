# Lagged calculated variables inside ODEs (#1176): symengine binds one to a
# bare symbol, so its definition is substituted before the sensitivity ODEs
# are derived.

#' Put the definitions of lagged variables into the ODEs that use them
#'
#' Run before the sensitivity ODEs are derived, so they do not treat a lagged
#' variable as a constant.  A variable defined more than once, or through a
#' history call, is left alone (see `.foceiLagInOde()`).
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
    if (length(.w) != 1L || .foceiLagRefs(.rhs[.w], .var, hist = TRUE)) {
      next
    }
    # rxToSE() is NSE: hand it a plain variable
    .txt <- .rhs[.w]
    .def <- symengine::S(rxode2::rxToSE(.txt))
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
