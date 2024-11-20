

#' This preprocesses the UI with any needed modifications
#'
#'
#' This is used to apply the mu referencing bug fix for `rxode2`
#'
#' @param ui rxode2 UI
#' @param control -- control list for the fit
#' @param est -- the fit estimation method
#' @return correct ui for say a fit (possibly a simulation)
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
.nlmixrPreprocessUi <- function(ui, control, est) {
  ui <- rxode2::assertRxUi(ui)
  .ret <- rxode2::rxUiDecompress(ui)
  .checkLiteralFix <- TRUE
  if (is.null(control)) {
  } else if (checkmate::testLogical(control$literalFix, any.missing=FALSE, len=1, null.ok=FALSE)) {
    .checkLiteralFix <- control$literalFix
  }
  if (.checkLiteralFix) {
    .ui <- try(rxode2::rxFixPop(ui, returnNull=TRUE))
    if (inherits(.ui, "try-error")) .ui <- NULL
    if (!is.null(.ui)) {
      .ret <- rxode2::rxUiDecompress(.ui)
      nlmixr2global$nlmixr2EstEnv$uiUnfix <- ui
    }
  }
  .zeroEtas <- .getZeroEtasFromModel(.ret)
  if (length(.zeroEtas) > 0) {
    nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- rxode2::rxUiDecompress(.ret)
    .minfo(paste0("the following etas are removed from the model since their initial estimates are zero: ",
           paste(.zeroEtas, collapse=", ")))
    .ret <- .downgradeEtas(ui, zeroEtas=.zeroEtas)
  }
  .ret
}
