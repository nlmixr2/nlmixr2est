#' Iteration-print configuration parameters (documentation stub)
#'
#' Documentation-only stub holding every argument that controls
#' iteration-time print formatting in `nlmixr2est`: both the scalar
#' `*Control()` arguments (`print`, `printNcol`, `useColor`) and the
#' [iterPrintControl()] arguments (`every`, `ncol`, `headerEvery`,
#' `useColor`, `simple`); downstream functions pull entries via
#' `@inheritParams iterPrintParams`.
#'
#' @param print Either a scalar print-frequency (`0` = suppress, `1`
#'   (default) = every evaluation, `N` = every Nth), OR a pre-built
#'   [iterPrintControl()] object. The scalar form is equivalent to
#'   `iterPrintControl(every = print, ncol = printNcol, useColor = useColor)`.
#' @param printNcol Integer (or `NULL`) — parameter columns per row before
#'   wrapping. `NULL` (default) uses `floor((getOption("width") - 23) / 12)`
#'   (fits an 80-column terminal). Equivalent to `print =
#'   iterPrintControl(ncol = ...)`.
#' @param every Integer. Print one iteration row every `every` parameter
#'   evaluations; `0` suppresses output. Defaults to `1L`. This is
#'   [iterPrintControl()]'s name for what `*Control()` calls `print`.
#' @param ncol Integer or `NULL`. Parameter columns per row before wrapping.
#'   `NULL` (default) uses `floor((getOption("width") - 23) / 12)`. This is
#'   [iterPrintControl()]'s name for what `*Control()` calls `printNcol`.
#' @param headerEvery Integer or `NULL`. Re-emit the column header every
#'   `headerEvery` parameter-print events (not raw iterations); `0` prints
#'   it once at fit start. `NULL` (default) uses `10L`.
#' @param useColor Logical (or `NULL`) — emit ANSI bold/color escapes in the
#'   iteration print. `NULL` (default) defers to [crayon::has_color()].
#'   Equivalent to `print = iterPrintControl(useColor = ...)`.
#' @param simple Logical. When `TRUE`, print a single row per iteration
#'   (optimizer-scale parameters only), suppressing the unscaled (`U`) /
#'   back-transformed (`X`) rows — used by estimators like saem with no
#'   internal optimizer scaling. Defaults to `FALSE`.
#' @return Nothing; this is a documentation-only helper.
#' @keywords internal
#' @name iterPrintParams
NULL

#' Control iteration-time print formatting
#'
#' Bundles every option controlling the iteration progress output emitted by
#' `nlmixr2` estimators; consumed by the shared C++ helper
#' (`scaleApplyIterPrintControl` in `src/scale.h`) so every estimator formats
#' its trace through one code path. Pass as the `print` argument to any
#' `*Control()` function, e.g. `foceiControl(print = iterPrintControl(every =
#' 5, headerEvery = 20))`; the scalar form `foceiControl(print = 5, printNcol
#' = 8)` still works and is wrapped into an `iterPrintControl()` internally.
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
#' `ui$iniDf` and emits, in one pass, every per-printed-parameter
#' transform vector any estimator's iteration printer or C-side
#' setup needs.  All output vectors are aligned to `printNames` —
#' index `i` describes the i-th printed parameter.
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
    xPar           = xPar,
    probitIdx      = probitIdx,
    logitThetaLow  = logitThetaLow,
    logitThetaHi   = logitThetaHi,
    probitThetaLow = probitThetaLow,
    probitThetaHi  = probitThetaHi
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
  list(
    xPar           = xPar,
    probitIdx      = probitIdx,
    logitThetaLow  = logitThetaLow,
    logitThetaHi   = logitThetaHi,
    probitThetaLow = probitThetaLow,
    probitThetaHi  = probitThetaHi
  )
}

#' Wrap scalar or list arguments into an iterPrintControl object
#'
#' Internal helper used by every `*Control()` function (and downstream
#' packages like `babelmixr2`) to absorb the scalar `print`/`printNcol`/
#' `useColor` arguments into a single [iterPrintControl()] object, or pass
#' through an already pre-built one. Exported under a leading-dot name to
#' mark it as package-internal, not for end users.
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
