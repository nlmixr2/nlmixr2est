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
getBaseSimModelFit.focei <- function(x) {
  obj <- x[[1]]
  .base <-   rxode2::rxCombineErrorLines(obj$ui)
  .base[[1]] <- quote(`rxModelVars`)
  .base <- eval(.base)
  # focei models are always pruned
  .malert("pruning simulation model...")
  .base <- rxode2::rxPrune(.base)
  .msuccess("done")
  .ctl <- obj$foceiControl
  if (.ctl$sumProd) {
    .malert("stabilizing round off errors in simulation problem...")
    .base <- rxode2::rxSumProdModel(.base)
    .msuccess("done")
  }
  if (.ctl$optExpression) {
    .base <- rxode2::rxOptExpr(.base, ifelse(.getRxPredLlikOption(),
                                             "simulation llik model",
                                             "simulation model"))
    .msuccess("done")
  }
  eval(parse(text=paste0("quote(rxode2({\n", .base, "\n}))")))
}

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.foce <- getBaseSimModelFit.focei

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.fo <- getBaseSimModelFit.focei
#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.foi <- getBaseSimModelFit.focei

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.posthoc <- getBaseSimModelFit.focei

#' @rdname getBaseSimModelFit
#' @export
getBaseSimModelFit.default <- function(x) {
  .obj <- x[[1]]
  .ui <- .obj$ui
  getBaseSimModel(.ui)
}

getBaseSimModel.nlmixr2FitCoreSilent <- function(obj) {
  .est <- obj$est
  .ret <- list(obj)
  class(.ret) <- c(.est, "getBaseSimModelFit")
  return(getBaseSimModelFit(.ret))
}

getBaseSimModel.nlmixr2FitData <- getBaseSimModel.nlmixr2FitCoreSilent
getBaseSimModel.nlmixr2FitCore <- getBaseSimModel.nlmixr2FitCoreSilent
