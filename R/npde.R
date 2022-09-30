#' NPDE calculation for nlmixr2
#'
#' @param object nlmixr2 fit object
#' @param updateObject Boolean indicating if original object should be updated.  By default this is TRUE.
#' @param table `tableControl()` list of options
#' @inheritParams foceiControl
#' @inheritParams addCwres
#' @param ... Other ignored parameters.
#'
#'
#' @return New nlmixr2 fit object
#' @author Matthew L. Fidler
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
#' f <- nlmixr2(one.cmt, theo_sd, "saem")
#'
#' # even though you may have forgotten to add the NPDE, you can add it to the data.frame:
#'
#' f <- addNpde(f)
#'
#' }
#' @export
addNpde <- function(object, updateObject = TRUE,
                    table = tableControl(), ...,
                    envir=parent.frame(1)) {
  assignInMyNamespace(".finalUiCompressed", FALSE)
  on.exit(assignInMyNamespace(".finalUiCompressed", TRUE))
  assertNlmixrFitData(object)
  if (any(names(object) == "NPDE")) {
    warning("already contains NPDE", call.=FALSE)
    return(object)
  }
  checkmate::assertLogical(updateObject, len=1, any.missing=FALSE)
  nlmixrWithTiming("NPDE", {
    .objName <- as.character(substitute(object))
    if (missing(table)) table <- object$table
    .malert("Add NPDE")
    if(missing(table)) {
      table <- object$table
    }
    table$npde <- TRUE
    .fitEnv <- object$env
    .npde <- .calcNpde(object, dv=object$DV, table=table)
    .fit <- nlmixrClone(object)
    .new <- nlmixrCbind(.fit, .npde[[2]])
    if (updateObject) {
      nlmixrUpdateObject(.new, .objName, envir, .fitEnv)
    }
    .msuccess("done")
    .new
  }, object$env)
}
