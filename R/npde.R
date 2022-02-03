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
  .pt <- proc.time()
  .objName <- substitute(object)
  rxode2::.setWarnIdSort(FALSE)
  on.exit(rxode2::.setWarnIdSort(TRUE))
  if (missing(table)) table <- object$table
  if (any(names(object) == "NPDE")) {
    warning("already contains NPDE")
    return(object)
  }
  .malert("Add NPDE")
  if(missing(table)) {
    table <- object$table
  }
  table$npde <- TRUE
  .npde <- .calcNpde(object, dv=object$DV, table=table)
  .cls <- class(object)
  .new <- cbind(object, .npde[[2]])
  class(.new) <- .cls
  .env <- .new$env
  if (inherits(updateObject, "logical")) {
    if (updateObject) {
      .parent <- envir
      .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
        if (.cur == .objName && identical(.parent[[.cur]]$env, .env)) {
          return(.cur)
        }
        return(NULL)
      }))
      if (length(.bound) == 1) {
        if (exists(.bound, envir = .parent)) {
          assign(.bound, .new, envir = .parent)
        }
      }
    }
  }
  .env$time <- .data.frame(.env$time, npde = (proc.time() - .pt)["elapsed"], check.names = FALSE)
  .msuccess("done")
  .new
}
