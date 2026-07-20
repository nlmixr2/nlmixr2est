#' Nudge exactly-zero theta initial estimates off zero (FOCEi family)
#'
#' FOCEi scales a linear parameter by its native magnitude `|init|`, which
#' is `0` (no scale) when a population parameter is initialized at exactly
#' `0` (the parameter then freezes).  This hook moves every estimated,
#' non-fixed `theta` whose initial estimate is exactly `0` to `+zeroTheta`
#' when that is within its bounds, otherwise `-zeroTheta`; if neither is
#' within the bounds it errors.  Only runs for the FOCEi family (a
#' `foceiControl`), and runs before `.preProcessBoundedTransform` so the
#' nudged value is what gets transformed.
#'
#' @param ui rxode2 ui model
#' @inheritParams nlmixr2
#' @return list with the ui (possibly modified)
#' @export
#' @author Matthew L. Fidler
.preProcessZeroTheta <- function(ui, est, data, control) {
  if (!inherits(control, "foceiControl")) return(NULL)
  .mag <- control$zeroTheta
  if (is.null(.mag) || !is.numeric(.mag) || length(.mag) != 1L ||
        !is.finite(.mag) || .mag <= 0) {
    return(NULL)
  }
  .ui <- rxode2::rxUiDecompress(ui)
  .iniDf <- .ui$iniDf
  .w <- which(!is.na(.iniDf$ntheta) & .iniDf$est == 0 & !.iniDf$fix)
  if (length(.w) == 0L) return(NULL)
  .changed <- character(0)
  for (.i in .w) {
    .lower <- .iniDf$lower[.i]
    .upper <- .iniDf$upper[.i]
    if (.mag > .lower && .mag < .upper) {
      .new <- .mag
    } else if (-.mag > .lower && -.mag < .upper) {
      .new <- -.mag
    } else {
      stop("cannot move the zero initial estimate of '", .iniDf$name[.i],
           "' off 0: neither ", .mag, " nor ", -.mag,
           " is within its bounds (", .lower, ", ", .upper, ")",
           call.=FALSE)
    }
    .iniDf$est[.i] <- .new
    .changed <- c(.changed, .iniDf$name[.i])
  }
  .ui$iniDf <- .iniDf
  .minfo(paste0("moved zero initial estimate(s) off 0 (foceiControl(zeroTheta)): ",
                paste(.changed, collapse=", ")))
  list(ui=.ui)
}
preProcessHooksAdd(".preProcessZeroTheta", .preProcessZeroTheta)
