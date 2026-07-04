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
    if (inherits(.v, "raw")) {
      ## Try reading in with qs2 if it doesn't work try with qs
      .c <- try(qs2::qs_deserialize(get(.v, envir=fit$env)))
      if (inherits(.c, "try-error")) {
        rxode2::rxReq("qs")
        .c <- rxode2::rxOldQsDes(get(.v, envir=fit$env))
        assign(.v, .c, envir=fit$env)
      }
    }
  }
  fit
}
