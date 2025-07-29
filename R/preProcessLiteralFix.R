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
  .checkLiteralFixRes <- FALSE
  .ui <- ui
  if (is.null(control)) {
  } else {
    if (checkmate::testLogical(control$literalFix,
                               any.missing=FALSE, len=1,
                               null.ok=FALSE)) {
      .checkLiteralFix <- control$literalFix
    } else if (is.null(control$literalFix)) {
      .checkLiteralFix <- FALSE
    }
    if (checkmate::testLogical(control$literalFixRes,
                               any.missing=FALSE, len=1,
                               null.ok=FALSE)) {
      .checkLiteralFixRes <- control$literalFixRes
    } else if (is.null(control$literalFixRes)) {
      .checkLiteralFixRes <- FALSE
    }
  }
  .assignUnfix <- FALSE
  if (.checkLiteralFix || .checkLiteralFixRes) {
    if (.checkLiteralFix) {
      .ui <- try(rxode2::rxFixPop(ui, returnNull=TRUE))
      if (inherits(.ui, "try-error")) .ui <- NULL
      if (!is.null(.ui)) {
        .ui <- rxode2::rxUiDecompress(.ui)
        .assignUnfix <- TRUE
      } else {
        .ui <- ui
      }
    }
    if (.checkLiteralFixRes && utils::packageVersion("rxode2") >= "4.0.0") {
      .ui2 <- .Call(`_rxode2rxFixRes`, .ui, TRUE)
      if (inherits(.ui2, "try-error")) .ui2 <- NULL
      if (!is.null(.ui2)) {
        .ui2 <- rxode2::rxUiDecompress(.ui2)
        .assignUnfix <- TRUE
      } else {
        .ui2 <- .ui
      }
      .ui <- .ui2
    }
  }
  if (.assignUnfix) {
    nlmixr2global$nlmixr2EstEnv$uiUnfix <- ui
  }
  list(ui=.ui)
}
preProcessHooksAdd(".nlmixrPreprocessLiteralFix", .nlmixrPreprocessLiteralFix)
