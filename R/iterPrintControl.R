#' Iteration-print configuration parameters (documentation stub)
#'
#' Documentation-only stub.  This is the single canonical location
#' for every argument that controls iteration-time print formatting
#' in `nlmixr2est`.  Two parallel sets of names live here:
#'
#' - The scalar arguments accepted by every `*Control()` function
#'   (`print`, `printNcol`, `useColor`), and
#' - The arguments of [iterPrintControl()] itself (`every`, `ncol`,
#'   `headerEvery`, `useColor`, `simple`).
#'
#' Both sets are documented here so a reader looking up any one
#' argument finds all of them — and downstream functions pull just
#' the entries they need via `@inheritParams iterPrintParams`.
#'
#' @param print Either a scalar print-frequency (`0` = suppress
#'   iteration output; `1` (default) = print every parameter
#'   evaluation; `N` = print every Nth evaluation), OR a pre-built
#'   [iterPrintControl()] object bundling all iteration-print options
#'   (column wrap, header cadence, color, simple/three-row mode).
#'   The scalar form is equivalent to
#'   `iterPrintControl(every = print, ncol = printNcol,
#'   useColor = useColor)`.
#' @param printNcol Integer (or `NULL`) — number of parameter columns
#'   emitted per row before wrapping to a continuation row.  `NULL`
#'   (the default) defers to [iterPrintControl()]'s default
#'   (`floor((getOption("width") - 23) / 12)`, which fits an 80-column
#'   terminal).  Equivalent to `print = iterPrintControl(ncol = ...)`.
#' @param every Integer.  Print one iteration row every `every`
#'   parameter evaluations.  `0` suppresses iteration output entirely.
#'   Defaults to `1L`.  Inside [iterPrintControl()] this is the
#'   canonical name for what the outer `*Control()` functions call
#'   `print`.
#' @param ncol Integer or `NULL`.  Number of parameter columns to
#'   emit per row before wrapping to a continuation row.  `NULL`
#'   (default) uses `floor((getOption("width") - 23) / 12)`, which
#'   fits an 80-column terminal.  Inside [iterPrintControl()] this
#'   is the canonical name for what the outer `*Control()` functions
#'   call `printNcol`.
#' @param headerEvery Integer or `NULL`.  Re-emit the column header
#'   every `headerEvery` parameter-print events (not raw iterations).
#'   With `every = 5` and `headerEvery = 10`, the header re-prints
#'   every 50 iterations.  `0` prints the header once at fit start
#'   only.  `NULL` (default) uses `10L`.
#' @param useColor Logical (or `NULL`) — whether to emit ANSI
#'   bold/color escapes in the iteration print.  `NULL` (the default)
#'   defers to [iterPrintControl()]'s default ([crayon::has_color()]).
#'   Equivalent to `print = iterPrintControl(useColor = ...)`.
#' @param simple Logical.  When `TRUE`, the printer emits a single
#'   row per iteration (just the optimizer-scale parameters) and
#'   suppresses the unscaled (`U`) / back-transformed (`X`) follow-up
#'   rows.  Used by estimators (like saem) that have no internal
#'   optimizer scaling, where U and X would be degenerate copies of
#'   the first row.  Defaults to `FALSE` (full three-row output).
#' @return Nothing; this is a documentation-only helper.
#' @keywords internal
#' @name iterPrintParams
NULL

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
#' All argument descriptions live in [iterPrintParams] so the same
#' text is shared with every `*Control()` function.
#'
#' @inheritParams iterPrintParams
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
                             useColor = NULL,
                             simple = FALSE) {
  if (is.null(ncol))        ncol        <- floor((getOption("width") - 23) / 12)
  if (is.null(headerEvery)) headerEvery <- 10L
  if (is.null(useColor))    useColor    <- crayon::has_color()
  checkmate::assertIntegerish(every,       len = 1, lower = 0, any.missing = FALSE)
  checkmate::assertIntegerish(ncol,        len = 1, lower = 1, any.missing = FALSE)
  checkmate::assertIntegerish(headerEvery, len = 1, lower = 0, any.missing = FALSE)
  checkmate::assertLogical(useColor,       len = 1, any.missing = FALSE)
  checkmate::assertLogical(simple,         len = 1, any.missing = FALSE)
  .ret <- list(
    every       = as.integer(every),
    ncol        = as.integer(ncol),
    headerEvery = as.integer(headerEvery),
    useColor    = as.logical(useColor),
    simple      = as.logical(simple)
  )
  class(.ret) <- c("iterPrintControl", "list")
  .ret
}

#' Derive every iteration-print transform vector from a ui object
#'
#' Pure inspection helper.  Walks `ui$muRefCurEval` against
#' `ui$iniDf` and emits, in one pass, every transform vector any
#' estimator's iteration printer or C-side setup needs.
#'
#' Print-vector-ordered (length == length(printNames), or all
#' thetas in `ntheta` order when `printNames` is NULL) — index `i`
#' describes the i-th printed parameter:
#'
#'   `xPar`              integer log/logit code per printed param:
#'                       `1`  = log-transformed; X shows `exp(value)`
#'                       `-m` = m-th logit-transformed parameter
#'                              (1-based); X shows
#'                              `expit(value, logitThetaLow[m-1],
#'                                    logitThetaHi[m-1])`
#'                       `0`  = no log/logit transform.
#'   `probitIdx`         integer probit index per printed param:
#'                       `k`  = k-th probit-transformed parameter
#'                              (1-based); X shows
#'                              `probitInv(value, probitThetaLow[k-1],
#'                                        probitThetaHi[k-1])`
#'                       `0`  = no probit transform.
#'                       Probit is carried as a parallel index rather
#'                       than folded into `xPar` so omega xPar codes
#'                       (2-5, used by `scaleGetScaleC`) cannot collide.
#'   `logitThetaLow`     lower bounds (one per logit entry, in the
#'                       occurrence order encoded in `xPar`).
#'   `logitThetaHi`      upper bounds, same ordering as logitThetaLow.
#'   `probitThetaLow`    lower bounds for probit, in occurrence order.
#'   `probitThetaHi`     upper bounds for probit, same ordering.
#'
#' ntheta-indexed (always populated from iniDf$ntheta, independent
#' of `printNames`):
#'
#'   `logNthetas`        ntheta values for log-transformed thetas.
#'   `logitNthetas`      ntheta values for logit-transformed thetas.
#'   `logitNthetasLow`   matching lower bounds.
#'   `logitNthetasHi`    matching upper bounds.
#'   `probitNthetas`     ntheta values for probit-transformed thetas.
#'   `probitNthetasLow`  matching lower bounds.
#'   `probitNthetasHi`   matching upper bounds.
#'
#' Function is pure — never mutates the ui or any environment.
#' Callers do their own assignments.
#'
#' Names in `printNames` not present in `ui$muRefCurEval`
#' (e.g. saem's `V(eta.*)` omega-variance names or residual-error
#' names) silently get `xPar = 0` and `probitIdx = 0`.
#'
#' @param ui rxode2 ui object.
#' @param printNames Character vector of parameter names in the same
#'   order as the printed parameter vector.  When `NULL` (default)
#'   uses all thetas (fixed + unfixed) in `ntheta` order — the form
#'   focei's C-side consumes for its theta-indexed back-transform.
#' @return Named list of integer / numeric vectors as documented above.
#' @noRd
.iterPrintXParFromUi <- function(ui, printNames = NULL) {
  iniThetas <- ui$iniDf[!is.na(ui$iniDf$ntheta), c("ntheta", "name")]
  iniThetas <- iniThetas[order(iniThetas$ntheta), ]
  if (is.null(printNames)) printNames <- iniThetas$name
  printNames <- as.character(printNames)
  muRef <- ui$muRefCurEval
  xPar <- integer(length(printNames))
  probitIdx <- integer(length(printNames))
  logitThetaLow <- numeric(0)
  logitThetaHi <- numeric(0)
  probitThetaLow <- numeric(0)
  probitThetaHi <- numeric(0)
  empty <- function() list(
    xPar             = xPar,
    probitIdx        = probitIdx,
    logitThetaLow    = logitThetaLow,
    logitThetaHi     = logitThetaHi,
    probitThetaLow   = probitThetaLow,
    probitThetaHi    = probitThetaHi,
    logNthetas       = integer(0),
    logitNthetas     = integer(0),
    logitNthetasLow  = numeric(0),
    logitNthetasHi   = numeric(0),
    probitNthetas    = integer(0),
    probitNthetasLow = numeric(0),
    probitNthetasHi  = numeric(0)
  )
  if (is.null(muRef) || nrow(muRef) == 0L) return(empty())
  if (!is.null(ui$boundedTransforms)) {
    for (.tr in ui$boundedTransforms) {
      .w <- which(muRef$parameter == .tr$internalName)
      if (length(.w) > 0L) {
        muRef$low[.w] <- .tr$lower
        muRef$hi[.w] <- .tr$upper
      }
    }
  }
  # Per-printed-name xPar / probitIdx / bounds, in printNames order.
  for (i in seq_along(printNames)) {
    nm <- printNames[i]
    idx <- which(muRef$parameter == nm)
    if (length(idx) == 0L) next
    idx <- idx[1L]
    ce <- muRef$curEval[idx]
    if (isTRUE(ce == "exp")) {
      xPar[i] <- 1L
    } else if (isTRUE(ce == "expit")) {
      logitThetaLow <- c(logitThetaLow, muRef$low[idx])
      logitThetaHi  <- c(logitThetaHi,  muRef$hi[idx])
      xPar[i] <- -as.integer(length(logitThetaLow))
    } else if (isTRUE(ce == "probitInv")) {
      probitThetaLow <- c(probitThetaLow, muRef$low[idx])
      probitThetaHi  <- c(probitThetaHi,  muRef$hi[idx])
      probitIdx[i] <- as.integer(length(probitThetaLow))
    }
  }
  # ntheta-indexed views over the unfixed-theta order.  These are
  # always emitted regardless of `printNames` because
  # .postEstimationBoundedTransform consumes them on the env under
  # logThetasF/logitThetasF/probitThetasF names.
  tr <- merge(iniThetas, muRef, by.x = "name", by.y = "parameter")
  tr <- tr[order(tr$ntheta), ]
  list(
    xPar             = xPar,
    probitIdx        = probitIdx,
    logitThetaLow    = logitThetaLow,
    logitThetaHi     = logitThetaHi,
    probitThetaLow   = probitThetaLow,
    probitThetaHi    = probitThetaHi,
    logNthetas       = as.integer(tr[which(tr$curEval == "exp"),       "ntheta"]),
    logitNthetas     = as.integer(tr[which(tr$curEval == "expit"),     "ntheta"]),
    logitNthetasLow  = as.double( tr[which(tr$curEval == "expit"),     "low"]),
    logitNthetasHi   = as.double( tr[which(tr$curEval == "expit"),     "hi"]),
    probitNthetas    = as.integer(tr[which(tr$curEval == "probitInv"), "ntheta"]),
    probitNthetasLow = as.double( tr[which(tr$curEval == "probitInv"), "low"]),
    probitNthetasHi  = as.double( tr[which(tr$curEval == "probitInv"), "hi"])
  )
}

#' Wrap scalar or list arguments into an iterPrintControl object
#'
#' Package-set internal helper used by every `*Control()` function in
#' `nlmixr2est` and downstream packages (e.g. `babelmixr2`) to absorb
#' the scalar `print` / `printNcol` / `useColor` arguments into a
#' single [iterPrintControl()] sub-list.  If the user already passed
#' a pre-built `iterPrintControl()` object via the `print` argument
#' (or, on round-trip, via an `iterPrintControl =` slot in `...`),
#' return it directly.
#'
#' Exported under a leading-dot name to mark it as a package-set
#' internal — callable from sibling packages in the nlmixr2 family
#' but not meant for end users.
#'
#' @param print Either an integer print-frequency or an
#'   `iterPrintControl` object.
#' @param printNcol,useColor Scalar `*Control()` arguments forwarded
#'   to [iterPrintControl()] only when `print` is a scalar.
#' @param iterPrintControl Optional pre-built [iterPrintControl()]
#'   object.  Wins over `print` and the other scalars when supplied —
#'   used by the round-trip case where a returned control list is
#'   passed back through `do.call(*Control, .ctl)`.
#' @return An `iterPrintControl` list.
#' @keywords internal
#' @export
.absorbIterPrintControl <- function(print = 1L,
                                    printNcol = NULL,
                                    useColor = NULL,
                                    iterPrintControl = NULL) {
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
