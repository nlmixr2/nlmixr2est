#' Does the installed rxode2 estimate tte() endpoints?
#'
#' Detected by feature, not version: dev rxode2 5.1.8.9000 builds exist both
#' with and without `tte()` endpoints.
#'
#' @return logical
#' @noRd
.nlmixr2estRxode2HasTte <- function() {
  exists(".rxTteOptions", envir = asNamespace("rxode2"), inherits = FALSE)
}

#' Does this ui have a tte() endpoint?
#'
#' @param ui rxode2 ui
#' @return logical
#' @noRd
.nlmixr2estUiHasTte <- function(ui) {
  if (!.nlmixr2estRxode2HasTte()) {
    return(FALSE)
  }
  .predDf <- tryCatch(ui$predDf, error = function(e) NULL)
  !is.null(.predDf) && any(.predDf$distribution == "tte")
}

#' Check the data of tte() endpoints before estimation
#'
#' A tte endpoint's records have `DV` 1 at an event (or the end of an interval
#' holding one) and 0 otherwise.  Unless `tte(origin = TRUE)`, the first record
#' of each subject is the start of observation, so an event there is not
#' scored; that is warned about.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the ui (unmodified; this hook only checks)
#' @noRd
.preProcessTte <- function(ui, est, data, control) {
  if (!is.data.frame(data)) {
    return(list(ui = ui))
  }
  if (!.nlmixr2estUiHasTte(ui)) {
    return(list(ui = ui))
  }
  .predDf <- ui$predDf
  .tteOptions <- utils::getFromNamespace(".rxTteOptions", "rxode2")
  .id <- .nlmixrDataCol(data, "id")
  if (is.null(.id)) {
    .id <- rep(1L, nrow(data))
  }
  .dv <- .nlmixrDataCol(data, "dv")
  .time <- .nlmixrDataCol(data, "time")
  if (nrow(.predDf) > 1L && is.null(.nlmixrDataCol(data, "cmt")) && is.null(.nlmixrDataCol(data, "dvid"))) {
    stop(
      "a tte() endpoint beside other endpoints requires a 'cmt' or 'dvid' column",
      call. = FALSE
    )
  }
  .obs <- .nlmixrDataObsRows(data)
  for (.i in which(.predDf$distribution == "tte")) {
    .pred1 <- .predDf[.i, ]
    .is <- .nlmixrDataEndpointRows(data, .predDf, .i, .obs)
    .v <- .dv[.is]
    if (any(is.na(.v) | !(.v %in% c(0, 1)))) {
      stop(
        "the 'DV' of tte() endpoint '",
        .pred1$cond,
        "' must be 0 or 1",
        call. = FALSE
      )
    }
    if (isTRUE(.tteOptions(ui, .pred1)$origin)) {
      next
    }
    .o <- order(.id[.is], .time[.is])
    .first <- !duplicated(.id[.is][.o])
    if (any(.v[.o][.first] == 1)) {
      warning(
        "the first record of tte() endpoint '",
        .pred1$cond,
        "' is an event for some subjects; it is the start of observation and is not scored ",
        "(add a DV=0 start record, or use tte(origin = TRUE))",
        call. = FALSE
      )
    }
  }
  list(ui = ui)
}
preProcessHooksAdd(".preProcessTte", .preProcessTte)
