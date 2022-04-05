#' @export
update.nlmixr2FitCore <- function(object, ...) {
  .nlmixr2savePipe(object)

  .modelLines <- rxode2::.quoteCallInfoLines(match.call(expand.dots = TRUE)[-(1:2)],
                                             envir = parent.frame(2))
  x <- object$ui
  .ret <- rxode2::.copyUi(x)
  rxode2::.modelHandleModelLines(.modelLines, .ret, modifyIni = TRUE,
                                 envir = parent.frame(2))
}

#' @export
update.nlmixr2FitData <- update.nlmixr2FitCore

#' @export
update.nlmixr2FitCoreSilent <- update.nlmixr2FitCore
