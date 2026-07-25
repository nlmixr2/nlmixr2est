#' Try to fix a nlmixr2 fit
#'
#' Re-evaluates the model function against the current version of rxode2, for
#' fits created with an older nlmixr2/rxode2 version.
#'
#' @param fit nlmixr2 fit object from a different version of nlmixr2.
#'
#' @return A nlmixr2 fit that has been (possibly) adjusted to work
#'   with the current version of nlmixr2.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' \dontrun{
#'   # requires the qs package to read an older nlmixr2 v3 fit (qs is no
#'   # longer on CRAN); regenerates the rxode2 model so it works again
#'   # fit <- readRDS(system.file("testfit_nlmixr3.rds", package = "nlmixr2est"))
#'   # fit <- try(nlmixr2fix(fit))
#'   # if (!inherits(fit, "try-error")) rxSolve(fit)
#' }
nlmixr2fix <- function(fit) {
  message("# This function is meant to load nlmixr2 fits from other versions")
  message("# To reproduce the fit, you need to use the same version of nlmixr2")
  print(fit$env$sessioninfo)
  message("\n")
  message("# If all else fails you can try to install the version of nlmixr2 used to create the fit\n")
  .ui <- fit$env$ui$fun
  .ui <- suppressMessages(.ui())
  assign("ui", .ui, envir = fit$env)
  for (.v in ls(fit$env, all.names=TRUE)) {
    .raw <- get(.v, envir=fit$env)
    if (inherits(.raw, "raw")) {
      .type <- rxode2::rxGetSerialType_(.raw)
      .c <- if (.type == "qs2") {
        try(.legacyQs2(.raw, "qs_deserialize"), silent=TRUE)
      } else if (.type == "qdata") {
        try(.legacyQs2(.raw, "qd_deserialize"), silent=TRUE)
      } else {
        try(rxode2::rxDeserialize(.raw), silent=TRUE)
      }
      if (!inherits(.c, "try-error")) {
        assign(.v, .c, envir=fit$env)
      }
    }
  }
  fit
}

#' Deserialize a legacy qs2/qdata blob from an old fit
#'
#' Fits saved by nlmixr2est <= 7.0.1 may hold qs2- or qdata-serialized
#' objects; qs2 is no longer a dependency, so look it up dynamically.
#'
#' @param raw raw vector to deserialize
#' @param fun qs2 function name to use
#' @return the deserialized object
#' @noRd
.legacyQs2 <- function(raw, fun) {
  if (!requireNamespace("qs2", quietly=TRUE)) {
    stop("this object was saved with 'qs2'; install qs2 to read it",
         call.=FALSE)
  }
  getExportedValue("qs2", fun)(raw)
}
