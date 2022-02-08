#' Update the nlmix2 object with new fit information
#'
#' @param fit nlmixr2 fit object to update in the environment
#' @param objName Name of the object
#' @param envir Environment to search
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @export
nlmixrUpdateObject <- function(fit, objName, envir) {
  .parent <- envir
  .env <- fit$env
  .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
    if (.cur == objName && identical(.parent[[.cur]]$env, .env)) {
      return(.cur)
    }
    return(NULL)
  }))
  if (length(.bound) == 1) {
    if (exists(.bound, envir = .parent)) {
      assign(.bound, fit, envir = .parent)
    }
  }
  invisible()
}
