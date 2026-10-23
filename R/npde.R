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
      .npde <- .calcNpde(object, dv = object$DV, table = table)
      .fit <- nlmixrClone(object)
      .new <- nlmixrCbind(.fit, .npdeFitRows(object, .npde[[2]]))
      if (updateObject) {
        .evRebound <- nlmixrUpdateObject(.new, .objName, envir, .fitEnv)
      }
      .msuccess("done")
      .new
    },
    object$env
  )
}

#' NPDE columns with one row per row of the fit
#'
#' The NPDE is computed for the fit's observations.  The other rows of its
#' table get NA, as in a fit run with `tableControl(npde = TRUE)`: those of a
#' subject with no observation, which the table reinserts from its data, and
#' the dosing rows of a table with them (`addDosing`, which adds `EVID`).  The
#' fit's table and the NPDE follow the same order for the observations.
#'
#' @param fit nlmixr2 fit
#' @param npde NPDE columns for the fit's observations
#' @return `npde`, with NA rows for the fit's other rows
#' @noRd
#' @author Matthew L. Fidler
.npdeFitRows <- function(fit, npde) {
  if (nrow(npde) == nrow(fit)) {
    return(npde)
  }
  .in <- as.character(fit$ID) %in% as.character(fit$idLvl)
  if ("EVID" %in% names(fit)) {
    .in <- .in & fit$EVID == 0
  }
  if (sum(.in) != nrow(npde)) {
    # without EVID, a table's dosing rows cannot be told from its observations
    stop(
      "NPDE has ",
      nrow(npde),
      " rows for the fit's ",
      sum(.in),
      " observations",
      if (!("EVID" %in% names(fit))) "; a table with dosing rows needs its 'EVID' column",
      call. = FALSE
    )
  }
  .ret <- npde[rep(NA_integer_, nrow(fit)), , drop = FALSE]
  .ret[.in, ] <- npde
  rownames(.ret) <- NULL
  .ret
}
