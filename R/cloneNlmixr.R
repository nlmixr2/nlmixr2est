.cloneEnv <- function(env) {
  .cls <- attr(env, "class")
  .env <- new.env(parent = emptyenv())
  for (.x in ls(env, all.names=TRUE)) {
    if (is.environment(get(.x, env))){
      assign(.x, .cloneEnv(get(.x, env)), .env)
    } else {
      assign(.x, get(.x, env), .env)
    }
  }
  attr(.env, "class") <- .cls
  return(.env)
}

#' Clone nlmixr environment
#'
#' @param x nlmixr fit
#' @return cloned nlmixr environment
#' @author Matthew L. Fidler
#' @examples
#' \dontrun{
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
#' nlmixrClone(f)
#'
#' }
#' @export
nlmixrClone <- function(x) {
  assertNlmixrFit(x)
  if (inherits(x, "nlmixr2FitData")) {
    .cls <- class(x)
    .df <- x
    class(.df) <- "data.frame"
    .df <- cbind(.df)
    .env <- .cloneEnv(attr(.cls,".foceiEnv"))
    attr(.cls, ".foceiEnv") <- .env
    class(.df) <- .cls
    .df
  } else {
    .cloneEnv(x)
  }
}
