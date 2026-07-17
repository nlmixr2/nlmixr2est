#' Warn about dosing properties on the wrong linCmt() compartment
#'
#' In a first-order absorption `linCmt()` model the dose enters the implied
#' `depot` compartment, so `f()`/`alag()`/`rate()`/`dur()` on `central` never
#' reach the solver and are silently ignored (nlmixr2est issue #358).  Warn so
#' the user moves the property to `depot`.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return `NULL` (only emits a warning)
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.preProcessLinCmtDosing <- function(ui, est, data, control) {
  .props <- try(ui$props, silent=TRUE)
  if (inherits(.props, "try-error") || !isTRUE(.props$linCmt)) {
    return(NULL)
  }
  # Only first-order absorption models have an implied 'depot' dosing
  # compartment distinct from 'central'; a pure IV linCmt doses 'central'
  # directly, where these properties are legitimate.
  if (!("depot" %in% .props$cmt)) {
    return(NULL)
  }
  .cmtProp <- .props$cmtProp
  if (!is.data.frame(.cmtProp) || nrow(.cmtProp) == 0L) {
    return(NULL)
  }
  .dosingProps <- c("f", "alag", "rate", "dur")
  .bad <- .cmtProp[.cmtProp$Compartment == "central" &
                     tolower(.cmtProp$Property) %in% .dosingProps, , drop=FALSE]
  if (nrow(.bad) == 0L) {
    return(NULL)
  }
  .pretty <- c(f="bioavailability 'f(central)'", alag="lag time 'alag(central)'",
               rate="'rate(central)'", dur="'dur(central)'")
  .msg <- unname(.pretty[tolower(.bad$Property)])
  warning("in a first-order absorption linCmt() model the dose enters the 'depot' compartment; ",
          paste(.msg, collapse=", "),
          " is ignored while solving -- apply it to 'depot' instead (e.g. 'f(depot)')",
          call.=FALSE)
  NULL
}

preProcessHooksAdd(".preProcessLinCmtDosing", .preProcessLinCmtDosing)
