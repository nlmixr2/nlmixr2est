#' @export
model.nlmixr2FitCore <- function(x, ..., append = FALSE,
                                 auto = getOption("rxode2.autoVarPiping", TRUE),
                                 envir = parent.frame()) {
  .nlmixr2savePipe(x)
  .modelLines <- rxode2::.quoteCallInfoLines(match.call(expand.dots = TRUE)[-(1:2)],
                                             envir = envir)
  .ret <- rxode2::.copyUi(x$ui)
  rxode2::.modelHandleModelLines(
    modelLines = .modelLines,
    rxui = .ret,
    append = append,
    auto = auto,
    modifyIni = FALSE,
    envir=envir
  )
}

#' @export
model.nlmixr2FitData <- model.nlmixr2FitCore

#' @export
model.nlmixr2FitCoreSilent <- model.nlmixr2FitCore
