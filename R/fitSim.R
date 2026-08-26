#' Method for getting simulation rxode2 classic models based on fits
#'
#' @param x list where first element is the fit.  The class represents the estimation method.
#' @return  model for fit$simulationModel
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
getBaseSimModelFit <- function(x) {
  UseMethod("getBaseSimModelFit")
}

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.default <- function(x) {
  .obj <- x[[1]]
  .ui <- .obj$ui
  rxode2::getBaseSimModel(.ui)
}

## The focei-family methods used to build a `predOnly`-based simulation model
## and then throw it away, and to call the default method twice -- once with
## its result discarded.  They were therefore already identical to the default
## in behavior, at about three times the cost (two extra lowerings of the fit
## per call, one of them a `rxNorm()` of the focei `predOnly` model).  They are
## kept as aliases so dispatch, and the exported methods, stay as they were.
#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.focei <- getBaseSimModelFit.default

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.foce <- getBaseSimModelFit.default

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.focep <- getBaseSimModelFit.default

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.fo <- getBaseSimModelFit.default

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.foi <- getBaseSimModelFit.default

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.posthoc <- getBaseSimModelFit.default

getBaseSimModel.nlmixr2FitCoreSilent <- function(obj) {
  .est <- obj$est
  .ret <- list(obj)
  class(.ret) <- c(.est, "getBaseSimModelFit")
  return(getBaseSimModelFit(.ret))
}

getBaseSimModel.nlmixr2FitData <- getBaseSimModel.nlmixr2FitCoreSilent
getBaseSimModel.nlmixr2FitCore <- getBaseSimModel.nlmixr2FitCoreSilent
