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
    .ui <- rxode2::rxode2(.ui)
  }
  list(ui=.ui)
}
preProcessHooksAdd(".preProcessDataUi", .preProcessDataUi)
