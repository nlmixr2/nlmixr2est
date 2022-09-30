#' Add CWRES
#'
#' This returns a new fit object with CWRES attached
#'
#' @param fit nlmixr2 fit without WRES/CWRES
#' @param focei Boolean indicating if the focei objective function is
#'   added.  If not the foce objective function is added.
#' @param updateObject Boolean indicating if the original fit object
#'   should be updated. By default this is true.
#' @param envir Environment that should be checked for object to
#'   update.  By default this is the global environment.
#' @return fit with CWRES
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'   ini({
#'     ## You may label each parameter with a comment
#'     tka <- 0.45 # Log Ka
#'     tcl <- log(c(0, 2.7, 100)) # Log Cl
#'     ## This works with interactive models
#'     ## You may also label the preceding line with label("label text")
#'     tv <- 3.45; label("log V")
#'     ## the label("Label name") works with all models
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' f <- try(nlmixr2(one.cmt, theo_sd, "saem"))
#'
#' print(f)
#'
#' # even though you may have forgotten to add the cwres, you can add it to the data.frame:
#'
#' if (!inherits(f, "try-error")) {
#'   f <- try(addCwres(f))
#'   print(f)
#' }
#'
#' # Note this also adds the FOCEi objective function
#' }
#' @author Matthew L. Fidler
#' @export
addCwres <- function(fit, focei=TRUE, updateObject = TRUE, envir = parent.frame(1)) {
  assignInMyNamespace(".finalUiCompressed", FALSE)
  on.exit(assignInMyNamespace(".finalUiCompressed", TRUE))
  assertNlmixrFitData(fit)
  checkmate::assertLogical(updateObject, len=1, any.missing=FALSE)
  checkmate::assertLogical(focei, len=1, any.missing=FALSE)
  if (any(names(fit) == "CWRES")) {
    return(fit)
  }
  nlmixrWithTiming("CWRES", {
    .objName <- as.character(substitute(fit))
    .foceiControl <- fit$foceiControl
    .foceiControl$maxOuterIterations <- 0L
    .foceiControl$maxInnerIterations <- 0L
    .foceiControl$etaMat = as.matrix(fit$eta[, -1])
    .foceiControl$compress <- FALSE
    .foceiControl$covMethod <- 0L
    .foceiControl$interaction <- focei
    .newFit <- nlmixr2(fit, data=nlme::getData(fit), est="focei",
                       control = .foceiControl)
    .extra <- setdiff(names(.newFit), names(fit))
    .extra <- as.data.frame(.newFit)[, .extra]
    .origFitEnv <- fit$env
    .fit <- nlmixrClone(fit)
    .new <- nlmixrCbind(.fit, .extra)
    .objDf <- .newFit$objDf
    .type <- rownames(.objDf)
    nlmixrAddObjectiveFunctionDataFrame(.new, .objDf, .type)
    if (updateObject) {
      nlmixrUpdateObject(.new, .objName, envir, .origFitEnv)
    }
    invisible(.new)
  },
  envir=fit)
}
#' @rdname nmObjGetData
#' @export
nmObjGetData.addCwres <- function(x, ...) {
  addCwres(x[[1]], updateObject = FALSE, envir=parent.frame(2))
}
attr(nmObjGetData.addCwres, "desc") <- "Add CWRES to object if needed"
