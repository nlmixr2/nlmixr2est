#' @export
model.nlmixr2FitCore <- function(x, ..., envir = parent.frame()) {
  .modelLines <- rxode2::.quoteCallInfoLines(match.call(expand.dots = TRUE)[-(1:2)],
                                             envir = envir)
  .ret <- rxode2::.copyUi(x$ui)
  rxode2::.modelHandleModelLines(.modelLines, .ret, modifyIni = FALSE,
                                 envir)
}

#' @export
model.nlmixr2FitData <- model.nlmixr2FitCore
