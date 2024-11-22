#' This literally fixes parameters in the model
#'
#' Whenever there is a fixed parameter in the model, the parameter is
#' replaced with the literal value inside of the model and dropped
#' from the `ini` block.  This only occurs when the
#' `control$literalFix=TRUE`.
#'
#' @param ui model function/object
#' @inheritParams nlmixr2
#' @return list with possibly updated ui
#' @export
#' @author Matthew L. Fidler
.nlmixrPreprocessLiteralFix <- function(ui, est, data, control) {
  .checkLiteralFix <- TRUE
  .ui <- ui
  if (is.null(control)) {
  } else if (checkmate::testLogical(control$literalFix,
                                    any.missing=FALSE, len=1,
                                    null.ok=FALSE)) {
    .checkLiteralFix <- control$literalFix
  }
  if (.checkLiteralFix) {
    .ui <- try(rxode2::rxFixPop(ui, returnNull=TRUE))
    if (inherits(.ui, "try-error")) .ui <- NULL
    if (!is.null(.ui)) {
      .ui <- rxode2::rxUiDecompress(.ui)
      nlmixr2global$nlmixr2EstEnv$uiUnfix <- ui
    } else {
      .ui <- ui
    }
  }
  list(ui=.ui)
}
preProcessHooksAdd(".nlmixrPreprocessLiteralFix", .nlmixrPreprocessLiteralFix)
