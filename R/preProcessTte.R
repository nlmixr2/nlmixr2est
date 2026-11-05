#' Does the installed rxode2 estimate tte() endpoints?
#'
#' rxode2 after 5.1.8 adds `tte()` endpoints; with an older rxode2 no tte
#' endpoint can be parsed, so every tte path stays inactive.
#'
#' @return logical
#' @noRd
.nlmixr2estRxode2HasTte <- function() {
  utils::packageVersion("rxode2") > "5.1.8"
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
  .nms <- tolower(names(data))
  .col <- function(.n) {
    .w <- which(.nms == .n)
    if (length(.w) == 0L) NULL else data[[.w[1]]]
  }
  .id <- .col("id")
  if (is.null(.id)) .id <- rep(1L, nrow(data))
  .dv <- .col("dv")
  .time <- .col("time")
  .evid <- .col("evid")
  .obs <- if (is.null(.evid)) rep(TRUE, nrow(data)) else .evid == 0
  .cmt <- .col("cmt")
  for (.i in which(.predDf$distribution == "tte")) {
    .pred1 <- .predDf[.i, ]
    .is <- .obs
    if (!is.null(.cmt)) {
      .is <- .is & !is.na(.cmt) & if (is.character(.cmt) || is.factor(.cmt)) {
        as.character(.cmt) == .pred1$cond
      } else {
        .cmt == .pred1$cmt
      }
    } else if (nrow(.predDf) > 1L) {
      stop("estimating a tte() endpoint beside other endpoints requires a 'cmt' column", call. = FALSE)
    }
    .v <- .dv[.is]
    if (any(is.na(.v) | !(.v %in% c(0, 1)))) {
      stop("the 'DV' of tte() endpoint '", .pred1$cond, "' must be 0 or 1", call. = FALSE)
    }
    if (isTRUE(.tteOptions(ui, .pred1)$origin)) next
    .o <- order(.id[.is], .time[.is])
    .first <- !duplicated(.id[.is][.o])
    if (any(.v[.o][.first] == 1)) {
      warning("the first record of tte() endpoint '", .pred1$cond,
              "' is an event for some subjects; it is the start of observation and is not scored ",
              "(add a DV=0 start record, or use tte(origin = TRUE))", call. = FALSE)
    }
  }
  list(ui = ui)
}
preProcessHooksAdd(".preProcessTte", .preProcessTte)
