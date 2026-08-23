# Data-dependent confirmation for the linCmt() sensitivity carry, and the
# UI-facing entry point.  "Time-varying" is a data property: with no data a
# covariate in ui$allCovs is only a CANDIDATE (varying = NA); with data,
# actual within-subject variation is confirmed (TRUE/FALSE).  The symbolic
# detection itself lives in foceiLinCmtCarryEligible.R.

#' Does a covariate vary within any subject in the data?
#'
#' @return `NA` when the data has no usable ID/covariate column, otherwise
#'   `TRUE`/`FALSE`
#' @noRd
.rxFoceiCarryCovVaries <- function(data, cov) {
  .n <- names(data)
  .idCol <- .n[tolower(.n) == "id"]
  if (length(.idCol) != 1L || !(cov %in% .n)) return(NA)
  any(vapply(split(data[[cov]], data[[.idCol]]),
             function(v) {
               v <- v[!is.na(v)]
               length(unique(v)) > 1L
             }, logical(1)))
}

#' NA without data, else whether any of the covariates varies within a subject
#' @noRd
.rxFoceiCarryVarying <- function(covs, data) {
  if (is.null(data)) return(NA)
  .v <- vapply(covs, function(cv) .rxFoceiCarryCovVaries(data, cv), logical(1))
  if (all(is.na(.v))) NA else isTRUE(any(.v, na.rm = TRUE))
}

#' Empty carry-eligibility result (fixed column layout)
#' @noRd
.rxFoceiCarryEmpty <- function() {
  data.frame(slot = integer(0), slotName = character(0),
             eta = character(0), etaName = character(0),
             covs = character(0), shape = character(0),
             formula = character(0), dEtaFormula = character(0),
             varying = logical(0), stringsAsFactors = FALSE)
}

#' Carry-eligible pairs straight from a model UI (test/entry convenience)
#'
#' @param ui rxode2 UI model
#' @inheritParams .rxFoceiLinCmtCarryEligible
#' @return see `.rxFoceiLinCmtCarryEligible`
#' @noRd
.foceiLinCmtCarryPairs <- function(ui, data = NULL,
                                   interpolation = c("locf", "nocb", "midpoint", "linear")) {
  .ui <- rxode2::assertRxUi(ui)
  .s <- .ui$foceiEtaS
  .etaVars <- paste0("ETA_", seq_len(.s$..maxEta), "_")
  .rxFoceiLinCmtCarryEligible(list(.ui), .s, .etaVars, data = data, # nolint: object_usage_linter.
                              interpolation = match.arg(interpolation))
}
