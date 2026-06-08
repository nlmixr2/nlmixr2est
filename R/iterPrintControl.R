#' Control iteration-time print formatting
#'
#' Bundles every option that controls the formatting of the iteration
#' progress output emitted by `nlmixr2` estimators (focei, saem,
#' bobyqa, nlm, optim, nls, nlminb, lbfgsb3c, n1qn1, newuoa, uobyqa).
#' The returned list is consumed by the shared C++ helper
#' (`scaleApplyIterPrintControl` in `src/scale.h`) so every estimator
#' formats its iteration trace through one code path.
#'
#' Pass this as the `print` argument to any of the `*Control()`
#' functions, e.g. `foceiControl(print = iterPrintControl(every = 5,
#' headerEvery = 20))`.  Equivalently, the existing scalar form
#' `foceiControl(print = 5, printNcol = 8)` continues to work — the
#' outer `*Control()` wraps your individual arguments into an
#' `iterPrintControl()` call internally.
#'
#' @param every Integer.  Print one iteration row every `every`
#'   parameter evaluations.  `0` suppresses iteration output entirely.
#'   Defaults to `1L`.
#' @param ncol Integer or `NULL`.  Number of parameter columns to emit
#'   per row before wrapping to a continuation row.  `NULL` (default)
#'   uses `floor((getOption("width") - 23) / 12)`, which fits an
#'   80-column terminal.
#' @param headerEvery Integer or `NULL`.  Re-emit the column header
#'   every `headerEvery` parameter-print events (not raw iterations).
#'   With `every = 5` and `headerEvery = 10`, the header re-prints
#'   every 50 iterations.  `0` prints the header once at fit start
#'   only.  `NULL` (default) uses `10L`.
#' @param useColor Logical or `NULL`.  Whether to emit ANSI bold/color
#'   escapes.  `NULL` (default) uses [crayon::has_color()].
#' @return A list with the validated, defaulted iteration-print
#'   options.  Has class `"iterPrintControl"` so the outer `*Control()`
#'   functions can distinguish a pre-built object from a scalar
#'   `print = N`.
#' @author Bill Denney, Matthew L. Fidler
#' @export
#' @examples
#' iterPrintControl()
#' iterPrintControl(every = 5, headerEvery = 0)
iterPrintControl <- function(every = 1L,
                             ncol = NULL,
                             headerEvery = NULL,
                             useColor = NULL) {
  if (is.null(ncol))        ncol        <- floor((getOption("width") - 23) / 12)
  if (is.null(headerEvery)) headerEvery <- 10L
  if (is.null(useColor))    useColor    <- crayon::has_color()
  checkmate::assertIntegerish(every,       len = 1, lower = 0, any.missing = FALSE)
  checkmate::assertIntegerish(ncol,        len = 1, lower = 1, any.missing = FALSE)
  checkmate::assertIntegerish(headerEvery, len = 1, lower = 0, any.missing = FALSE)
  checkmate::assertLogical(useColor,       len = 1, any.missing = FALSE)
  .ret <- list(
    every       = as.integer(every),
    ncol        = as.integer(ncol),
    headerEvery = as.integer(headerEvery),
    useColor    = as.logical(useColor)
  )
  class(.ret) <- c("iterPrintControl", "list")
  .ret
}

#' Wrap scalar or list arguments into an iterPrintControl object
#'
#' Internal helper used by every `*Control()` function to absorb
#' historical scalar arguments (`print`, `printNcol`, `useColor`) into
#' a single [iterPrintControl()] sub-list.  If the user already passed
#' a pre-built `iterPrintControl()` object via the `print` argument,
#' return it directly (and warn if any of the historical scalar args
#' were also set, since they would be silently ignored).
#'
#' @param print Either an integer print-frequency or an
#'   `iterPrintControl` object.
#' @param printNcol,useColor Scalar arguments that historically lived
#'   on each `*Control()` function.  Forwarded to [iterPrintControl()]
#'   only when `print` is a scalar.
#' @param printHeader Deprecated.  Was briefly exposed as a top-level
#'   `*Control()` argument; now lives only in [iterPrintControl()]
#'   (the `headerEvery` argument there).  Accepted here only to issue
#'   a clear error pointing users at the new location.
#' @return An `iterPrintControl` list.
#' @noRd
.absorbIterPrintControl <- function(print = 1L,
                                    printNcol = NULL,
                                    useColor = NULL,
                                    printHeader = NULL,
                                    printUseColor = NULL,
                                    iterPrintControl = NULL) {
  if (!is.null(printHeader)) {
    stop("`printHeader` is no longer a top-level *Control() argument; pass ",
         "`print = iterPrintControl(headerEvery = ", printHeader, ")` instead.",
         call. = FALSE)
  }
  if (!is.null(printUseColor)) {
    stop("`printUseColor` is no longer a top-level *Control() argument; pass ",
         "`print = iterPrintControl(useColor = ", printUseColor, ")` instead.",
         call. = FALSE)
  }
  # If a fully-built iterPrintControl object came in via the dedicated
  # `iterPrintControl` argument (typical when a *Control() return value
  # round-trips back through do.call(*Control, .ctl)), use it directly.
  if (!is.null(iterPrintControl)) {
    if (!inherits(iterPrintControl, "iterPrintControl")) {
      stop("`iterPrintControl` must be the result of iterPrintControl()",
           call. = FALSE)
    }
    return(iterPrintControl)
  }
  if (inherits(print, "iterPrintControl")) {
    .conflicts <- character(0)
    if (!is.null(printNcol)) .conflicts <- c(.conflicts, "printNcol")
    if (!is.null(useColor))  .conflicts <- c(.conflicts, "useColor")
    if (length(.conflicts)) {
      warning("ignoring `", paste(.conflicts, collapse = "`, `"),
              "` because `print` was passed as an iterPrintControl object",
              call. = FALSE)
    }
    return(print)
  }
  iterPrintControl(every = print, ncol = printNcol, useColor = useColor)
}
