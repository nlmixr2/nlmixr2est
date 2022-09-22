#' Update the nlmixr2 object with new fit information
#'
#' @param fit nlmixr2 fit object to update in the environment
#' @param objName Name of the object
#' @param envir Environment to search
#' @param origFitEnv Original fit$env to compare, otherwise simply use fit$env
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
nlmixrUpdateObject <- function(fit, objName, envir, origFitEnv=NULL) {
  .parent <- envir
  if (is.environment(origFitEnv)) {
    .env <- origFitEnv
  } else {
    .env <- fit$env
  }
  .ls <- ls(.parent, all.names = TRUE)
  .bound <- do.call("c", lapply(.ls, function(.cur) {
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
