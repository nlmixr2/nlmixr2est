#' @export
as.function.nlmixr2FitData <- function(x, ...) {
  x$ui$fun
}

#' @export
as.function.nlmixr2FitCore <- as.function.nlmixr2FitData
