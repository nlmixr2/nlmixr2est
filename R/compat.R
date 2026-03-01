#' Try to fix a nlmixr2 fit
#'
#' Currently this re-evaluates the function in the current version of rxode2.
#'
#' @param fit nlmixr2 fit object from a different version of nlmixr2.
#'
#' @return A nlmixr2 fit that has been (possibly) adjusted to work
#'   with the current version of nlmixr2.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#'
#'
#'   # This is a nlmixr2 v3 fit and requires the qs package to read in
#'   fit <- system.file("testfit_nlmixr3.rds", package = "nlmixr2est")
#'   fit <- readRDS(fit)
#'
#'   # While it prints well, it can't be used in all functions because
#'   # Language features (like +var()) are not supported in the v3 version
#'
#'   print(fit)
#'
#'   try(rxSolve(fit)) # should error, but with try it will just display the error
#'
#'   # This function attempts to fix it by regenerating the rxode2 model with the
#'   # new features
#'
#'   # This function also prints out the information on how this fit was created
#'
#'   fit <- try(nlmixr2fix(fit))
#'
#'   # Now solving and other functions work
#'   if (!inherits(fit, "try-error")) {
#'     rxSolve(fit)
#'   }
#'
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
