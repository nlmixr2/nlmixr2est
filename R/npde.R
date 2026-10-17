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
addNpde <- function(object, updateObject = TRUE, table = tableControl(), ..., envir = parent.frame(1)) {
  nlmixr2global$finalUiCompressed <- FALSE
  on.exit(nlmixr2global$finalUiCompressed <- TRUE)
  assertNlmixrFitData(object)
  ## event bus: one fitUpdate on exit (see rxEvents.R)
  .nlmixr2EventEnter()
  .evOrig <- object
  .evName <- .nlmixr2EventObjName(substitute(object))
  .evRebound <- FALSE
  on.exit(.nlmixr2EventExitUpdate(returnValue(), .evOrig, .evName, "npde", .evRebound), add = TRUE)
  if (any(names(object) == "NPDE")) {
    warning("already contains NPDE", call. = FALSE)
    return(object)
  }
  checkmate::assertLogical(updateObject, len = 1, any.missing = FALSE)
  nlmixrWithTiming(
    "NPDE",
    {
      .objName <- as.character(substitute(object))
      if (missing(table)) {
        table <- object$table
      }
      .malert("Add NPDE")
      if (missing(table)) {
        table <- object$table
      }
      table$npde <- TRUE
      .fitEnv <- object$env
      if (.uiHasIov(object$ui)) {
        .npde <- .addNpdeIovRefit(object, table)
      } else {
        .npde <- .calcNpde(object, dv = object$DV, table = table)[[2]]
      }
      .fit <- nlmixrClone(object)
      .new <- nlmixrCbind(.fit, .npde)
      if (updateObject) {
        .evRebound <- nlmixrUpdateObject(.new, .objName, envir, .fitEnv)
      }
      .msuccess("done")
      .new
    },
    object$env
  )
}

#' Whether a model has an occasion-level (IOV) eta
#'
#' An eta whose base condition is not `"id"`, the rows `.uiApplyIov()`
#' expands; a repeated subject-level block (`"id:same:<eta>"`) is not one.
#'
#' @param ui rxode2 user interface model
#' @return `TRUE` when the model has an occasion eta
#' @noRd
#' @author Matthew L. Fidler
.uiHasIov <- function(ui) {
  .iniDf <- ui$iniDf
  .cnd <- .iniDf$condition
  any(!is.na(.cnd) & is.na(.iniDf$err) & lotri::lotriBaseCondition(.cnd) != "id")
}

#' A focei refit of a fit with IOV, for its tables
#'
#' A finished IOV fit keeps the IOV-expanded model, but its thetas, etas and
#' data are the user's model's: the occasion etas, their magnitude thetas and
#' the occasion column are gone, so `.calcTables()` and `.calcNpde()` cannot
#' solve it.  A focei refit with no iterations at the fit's etas, as
#' `addCwres()` refits, rewrites the model again, so its table step sees the
#' expanded model.
#'
#' @param fit nlmixr2 fit with IOV
#' @param table `tableControl()` for the refit's tables
#' @return the refit
#' @noRd
#' @author Matthew L. Fidler
.iovRefit <- function(fit, table) {
  .ctl <- fit$foceiControl
  .ctl$maxOuterIterations <- 0L
  .ctl$maxInnerIterations <- 0L
  .ctl$etaMat <- fit$etaMat
  .ctl$compress <- FALSE
  .ctl$covMethod <- 0L
  .ctl$nAGQ <- 0L
  # the refit exists for its tables, whether or not the fit computed its own
  .ctl$calcTables <- TRUE
  .nlmixr2PriorGateBypass(
    nlmixr2(fit, data = nlme::getData(fit), est = "focei", control = .ctl, table = table)
  )
}

#' The NPDE columns of a fit with IOV, from `.iovRefit()`
#'
#' @param fit nlmixr2 fit with IOV
#' @param table `tableControl()` with `npde = TRUE`
#' @return data frame of the refit's NPDE columns
#' @noRd
#' @author Matthew L. Fidler
.addNpdeIovRefit <- function(fit, table) {
  .newFit <- .iovRefit(fit, table)
  # the columns _nlmixr2est_npdeCalc() (src/npde.cpp) returns; the refit's
  # other new columns (a focei refit's CWRES) are not NPDE's
  .cols <- c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD")
  .miss <- setdiff(.cols, names(.newFit))
  if (length(.miss) > 0L) {
    stop("the IOV refit for NPDE did not give: ", paste(.miss, collapse = ", "), call. = FALSE)
  }
  as.data.frame(.newFit)[, .cols, drop = FALSE]
}
