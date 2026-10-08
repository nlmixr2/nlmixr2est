# Lagged calculated variables (#1176): which ones a model has, the lhs lines
# that define them, and where they are used.

.foceiHistFn <- c("lag", "lead", "diff", "first", "last", "lag0", "lead0", "diff0")

#' The lhs lines that define lagged calculated variables
#'
#' rxode2 reads an earlier value of a reassigned lagged variable through a
#' `rx_lagv<i>_<var>` snapshot, so those lines are kept too, in model order.
#' @param s symengine environment
#' @param ar when `FALSE`, drop the AR(1) residual's own lagged variables
#'   (`rx_ar*`); the AR(1) gradient correction handles those
#' @return the `var=expr` lines of `s$..lhs` that define them
#' @noRd
.foceiLagLhs <- function(s, ar = TRUE) {
  .v <- s$..laggedVars
  if (!ar) {
    .v <- .v[!grepl("^rx_ar", .v)]
  }
  if (length(.v) == 0L || is.null(s$..lhs)) {
    return(character(0))
  }
  .n <- sub("=.*$", "", s$..lhs)
  .snap <- grepl("^rx_lagv[0-9]+_", .n)
  s$..lhs[.n %in% .v | (.snap & sub("^rx_lagv[0-9]+_", "", .n) %in% .v)]
}

#' The definitions of the lagged calculated variables
#'
#' @param s symengine environment
#' @return `.foceiLagLhs()` without the AR(1) residual's variables
#' @noRd
.foceiLagDefs <- function(s) {
  .foceiLagLhs(s, ar = FALSE)
}

#' Whether rxode2 text uses any of the given variables
#'
#' @param txt rxode2 text
#' @param vars variable names
#' @param hist when `TRUE`, only count a use as a history function argument
#' @return logical
#' @noRd
.foceiLagRefs <- function(txt, vars, hist = FALSE) {
  .pre <- if (hist) paste0("\\b(", paste(.foceiHistFn, collapse = "|"), ")\\(\\s*") else "(?<![A-Za-z0-9_.])"
  .post <- if (hist) "\\s*[,)]" else "(?![A-Za-z0-9_.])"
  any(vapply(
    vars,
    function(v) any(grepl(paste0(.pre, gsub(".", "\\.", v, fixed = TRUE), .post), txt, perl = TRUE)),
    logical(1)
  ))
}

#' Whether a model uses a history function of a calculated variable
#'
#' A dosing `lag(cmt) <-` is not a history function, and a covariate's
#' history needs no sensitivity.
#' @param ui rxode2 UI
#' @return logical
#' @noRd
.foceiUsesLagVar <- function(ui) {
  .lhs <- tryCatch(rxode2::rxModelVars(ui)$lhs, error = function(e) character(0))
  .lhs <- .lhs[!grepl("^rx_ar", .lhs)]
  if (length(.lhs) == 0L) {
    return(FALSE)
  }
  .found <- FALSE
  .walk <- function(e) {
    if (.found || !is.call(e)) {
      return(invisible())
    }
    .f <- e[[1]]
    if (is.name(.f) && as.character(.f) %in% c("<-", "=", "~")) {
      .walk(e[[3]])
      return(invisible())
    }
    if (
      is.name(.f) &&
        as.character(.f) %in% .foceiHistFn &&
        length(e) >= 2L &&
        is.name(e[[2]]) &&
        as.character(e[[2]]) %in% .lhs
    ) {
      .found <<- TRUE
      return(invisible())
    }
    lapply(as.list(e)[-1], .walk)
    invisible()
  }
  .lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  lapply(.lst, .walk)
  .found
}
