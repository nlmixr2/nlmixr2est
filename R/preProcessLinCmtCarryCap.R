#' Carry columns a `linCmt()` model's eta-gradient carry needs
#'
#' Zero when the model cannot have more than `cap` (no `linCmt()`, no
#' covariate, or too few etas to fill the columns), so only a model that may
#' be past the cap pays for the symbolic eligibility check.
#'
#' @param ui rxode2 ui
#' @param cap the carry columns the loaded rxode2 has
#' @return integer, the columns needed
#' @noRd
#' @author Matthew L. Fidler
.rxFoceiLinCmtCarryColumnsNeeded <- function(ui, cap) {
  if (rxode2::.rxLinNcmt(ui)["numLin"] <= 0L || length(ui$allCovs) == 0L) {
    return(0L)
  }
  .iniDf <- ui$iniDf
  .nEta <- sum(!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2)
  # one column per pair, one tracker per slot (7 slots), one amounts and one
  # lag column at most
  if (.nEta + min(.nEta, 7L) + 2L <= cap) {
    return(0L)
  }
  .pairs <- tryCatch(.foceiLinCmtCarryPairs(ui), error = function(e) NULL)
  .pairs <- .rxFoceiLinCmtCarryJumpFilter(.pairs)
  if (is.null(.pairs)) {
    return(0L)
  }
  as.integer(.rxFoceiLinCmtCarryColumns(.pairs))
}

#' Solve a `linCmt()` model past the carry cap as ODEs
#'
#' The eta-gradient carry has `.rxFoceiLinCmtCarryMaxPairs()` columns.  The
#' etas past them would keep the standard gradient, which is wrong when a
#' covariate changes a `linCmt()` parameter within a subject and makes the
#' objective depend on the order the etas are declared in (#1318).  The ODE
#' translation has exact sensitivities for every eta.  Named to run after
#' `.uiApplyIov()`, whose occasion etas count toward the cap.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the translated ui, or `NULL` when nothing to do
#' @export
#' @author Matthew L. Fidler
.uiApplyLinCmtCarryCap <- function(ui, est, data, control) {
  if (!(est %in% .linCmtFoceiEstFamily)) {
    return(NULL)
  }
  # .preProcessLinCmtOde() translates these
  if (.uiIsMixedLinCmtOde(ui)) {
    return(NULL)
  }
  .ctl <- if (is.list(control)) control else list()
  .carry <- .ctl$linCmtSensCarry
  if (is.null(.carry)) {
    .carry <- "auto"
  }
  if (!.rxFoceiLinCmtCarryOn(.carry, .ctl$rxControl$covsInterpolation)) {
    return(NULL)
  }
  .cap <- .rxFoceiLinCmtCarryMaxPairs()
  .need <- .rxFoceiLinCmtCarryColumnsNeeded(ui, .cap)
  if (.need <= .cap) {
    return(NULL)
  }
  .ui <- try(rxode2::linToOde(ui), silent = TRUE)
  if (inherits(.ui, "try-error")) {
    return(NULL)
  }
  warning(
    sprintf("linCmt() solved as ODEs: etas need %d of %d carry columns", .need, .cap),
    call. = FALSE
  )
  list(ui = .ui)
}

preProcessHooksAdd(".uiApplyLinCmtCarryCap", .uiApplyLinCmtCarryCap)
