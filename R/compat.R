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
  # repair raw slots first; `ui` itself may be serialized in an old fit
  for (.v in ls(fit$env, all.names=TRUE)) {
    .raw <- get(.v, envir=fit$env)
    if (inherits(.raw, "raw")) {
      .c <- try(.deserializeRaw(.raw), silent=TRUE)
      if (!inherits(.c, "try-error")) {
        assign(.v, .c, envir=fit$env)
      }
    }
  }
  .ui <- fit$env$ui$fun
  .ui <- suppressMessages(.ui())
  assign("ui", .ui, envir = fit$env)
  fit
}

#' Deserialize a raw fit component by its serialization tag
#'
#' Fits saved by nlmixr2est <= 7.0.1 may hold qs2-, qdata- or qs-serialized
#' objects; neither qs2 nor qs is a dependency any more, so look them up
#' dynamically when installed.
#'
#' @param raw raw vector to deserialize
#' @param type serialization tag (from `rxode2::rxGetSerialType_()`)
#' @return the deserialized object; errors when the format needs a package
#'   that is not installed
#' @noRd
.deserializeRaw <- function(raw, type=rxode2::rxGetSerialType_(raw)) {
  switch(type,
         qs2 = .legacyQsFn("qs2", "qs_deserialize")(raw),
         qdata = .legacyQsFn("qs2", "qd_deserialize")(raw),
         qs = .legacyQsFn("qs", "qdeserialize")(raw),
         xz = unserialize(memDecompress(raw, type="xz")),
         bzip2 = unserialize(memDecompress(raw, type="bzip2")),
         base = unserialize(raw),
         stop("unknown serialization type '", type, "'", call.=FALSE))
}

#' @param pkg legacy serialization package ("qs2" or "qs")
#' @param fun function name to look up in `pkg`
#' @return the function, when `pkg` is installed
#' @noRd
.legacyQsFn <- function(pkg, fun) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    stop("this object was saved with '", pkg, "'; install ", pkg,
         " to read it", call.=FALSE)
  }
  getExportedValue(pkg, fun)
}
