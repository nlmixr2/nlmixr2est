#' nlmixrCbind
#'
#' `cbind` for `nlmixr` objects that preserve the fit information
#'
#' @param fit nlmixr fit
#' @param extraData data to merge to nlmixr fit
#' @return fit expanded with extra values
#' @author Matthew L. Fidler
#' @export
nlmixrCbind <- function(fit, extra) {
  if (!inherits(fit, "nlmixr2FitCore")) stop("'fit' must be a nlmixr2 fit", call.=FALSE)
  .cls <- class(fit)
  .dat <- fit
  class(.dat) <- "data.frame"
  .dat <- cbind(.dat, extra)
  class(.dat) <- .cls
  .dat
}
