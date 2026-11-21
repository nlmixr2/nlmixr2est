#' Preprocess the zero omegas
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the ui (possibly modified)
#' @export
#' @author Matthew L. Fidler
.preProcessDataUi <- function(ui, est, data, control) {
  .ui <- ui
  if (isTRUE(.ui$uiUseData)) {
    rxode2::rxUdfUiReset()
    on.exit({
      rxode2::rxUdfUiReset()
    })
    if (rxode2::rxIs(data, "event.data.frame")) {
      rxode2::rxUdfUiData(data)
    }
    rxode2::rxUdfUiEst(est)
    # Add for next version of rxode2
    ## rxode2::rxUdfUiControl(control)
    # user functions may register likelihood penalties while expanding
    nlmixr2global$likPenaltyUdf <- list()
    on.exit(nlmixr2global$likPenaltyUdf <- NULL, add = TRUE)
    .pen <- .likPenaltyAsList(.ui$nlmixr2LikPenalty)
    rxode2::rxUdfUiMv(rxode2::rxModelVars(.ui))
    # an rxUi is returned as-is by rxode2(); re-parse from its function
    .ui <- suppressMessages(rxode2::rxode2(as.function(.ui)))
    if (length(.pen) > 0L) {
      .ui <- do.call(nlmixr2LikPenaltyAdd, c(list(.ui), .pen, list(replace = TRUE)))
    }
    .ui <- .likPenaltyUdfAttach(.ui, nlmixr2global$likPenaltyUdf)
  }
  list(ui = .ui)
}
preProcessHooksAdd(".preProcessDataUi", .preProcessDataUi)
