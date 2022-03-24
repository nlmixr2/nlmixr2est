#' @export
update.nlmixr2FitCore <- function(object, ...) {
  .modelLines <- rxode2::.quoteCallInfoLines(match.call(expand.dots = TRUE)[-(1:2)],
                                             envir = envir)
  x <- object$ui
  .ret <- rxode2::.copyUi(x)
  rxode2::.modelHandleModelLines(.modelLines, .ret, modifyIni = TRUE,
                                 envir = envir)
}

#' @export
update.nlmixr2FitData <- update.nlmixr2FitCore

#' @export
update.nlmixr2FitCoreSilent <- update.nlmixr2FitCore
