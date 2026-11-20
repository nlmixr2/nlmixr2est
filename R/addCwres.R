#' Integrates items like phiC/phiH into focei environment
#'
#' @param env environment where focei information needs to be integrated
#' @param newFit new fit object with focei information in it.
#' @return Nothing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.addFoceiInfoToFit <- function(env, newFit) {
  for (.v in c("phiC", "phiH", "llikObs")) {
    if (exists(.v, envir = newFit$env, inherits = FALSE)) {
      assign(.v, get(.v, envir = newFit$env, inherits = FALSE), envir = env)
    }
  }
}

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
addCwres <- function(fit, focei = TRUE, updateObject = TRUE, envir = parent.frame(1)) {
  nlmixr2global$finalUiCompressed <- FALSE
  on.exit(nlmixr2global$finalUiCompressed <- TRUE)
  assertNlmixrFitData(fit)
  checkmate::assertLogical(updateObject, len = 1, any.missing = FALSE)
  checkmate::assertLogical(focei, len = 1, any.missing = FALSE)
  ## event bus: `fit` is forced above, outside the scope; the internal table
  ## and objective below are silent and one fitUpdate is emitted on exit (see
  ## rxEvents.R)
  .nlmixr2EventEnter()
  .evOrig <- fit
  .evName <- .nlmixr2EventObjName(substitute(fit))
  .evRebound <- FALSE
  on.exit(.nlmixr2EventExitUpdate(returnValue(), .evOrig, .evName, "cwres", .evRebound), add = TRUE)
  if (is.null(fit$eta)) {
    stop("cannot add CWRES to a model without etas", call. = FALSE)
  } else if (any(names(fit) == "CWRES")) {
    return(fit)
  }
  nlmixrWithTiming(
    "CWRES",
    {
      .objName <- as.character(substitute(fit))
      # the fit's own table with CWRES, then the objective setOfv() adds
      .table <- fit$table
      if (!inherits(.table, "tableControl")) {
        .table <- do.call(tableControl, as.list(.table))
      }
      .table$cwres <- TRUE
      # addTable() updates the environment of the fit it is given
      .new <- nlmixrClone(fit)
      .tab <- .nlmixr2EventScope(addTable(.new, updateObject = FALSE, table = .table))
      .extra <- setdiff(names(.tab), names(fit))
      .extra <- as.data.frame(.tab)[, .extra, drop = FALSE]
      .origFitEnv <- fit$env
      .new <- nlmixrCbind(.new, .extra)
      .setOfvFo(.new, if (focei) "focei" else "foce")
      if (updateObject) {
        .evRebound <- nlmixrUpdateObject(.new, .objName, envir, .origFitEnv)
      }
      invisible(.new)
    },
    envir = fit
  )
}
#' @rdname nmObjGetData
#' @export
nmObjGetData.addCwres <- function(x, ...) {
  ## an accessor, not a user update: silent on the event bus
  .nlmixr2EventScope(addCwres(x[[1]], updateObject = FALSE, envir = parent.frame(2)))
}
attr(nmObjGetData.addCwres, "desc") <- "Add CWRES to object if needed"
