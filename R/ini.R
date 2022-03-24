#' @export
ini.nlmixr2FitCore <- function(x, ..., envir = parent.frame()) {
  .ret <- rxode2::.copyUi(x$ui)
  .iniLines <- rxode2::.quoteCallInfoLines(match.call(expand.dots = TRUE)[-(1:2)],
                                           envir = envir)
  lapply(.iniLines, function(line) {
    rxode2::.iniHandleFixOrUnfix(line, .ret, envir = envir)
  })
  .ret
}

#' @export
ini.nlmixr2FitData <- ini.nlmixr2FitCore

#' @export
ini.nlmixr2FitCoreSilent <- ini.nlmixr2FitCore
