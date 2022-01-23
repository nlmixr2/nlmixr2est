##' Generic for nlmixr2 estimation methods
##'
##' @param env Environment for the nlmixr2 estimation routines.
##'
##' This needs to have:
##'
##' - rxode2 ui object in `$ui`
##'
##' - data to fit in the estimation routine in `$data`
##'
##' - control for the estimation routine's control options in `$ui`
##'
##' @return nlmixr2 fit object
##'
##' @author Matthew Fidler
##'
##' @details
##'
##' This is a S3 generic that allows others to use the nlmixr2
##'   environment to do their own estimation routines
##'
##' @export
nlmixr2Est <- function(env, ...) {
  if (!exists("ui", envir=env)) {
    stop("need 'ui' object", call.=FALSE)
  } else if (!inherits(env$ui, "rxUi")) {
    stop("'ui' is not an rxode2 object", call.=FALSE)
  }
  if (!exists("data", envir=env)) {
    stop("need 'data' object", call.=FALSE)
  } else if (!inherits(env$data, "data.frame")) {
    stop("'data' is not a data.frame", call.=FALSE)
  }
  env$data <- as.data.frame(env$data)
  if (!exists("control", envir=env)) {
    stop("need 'control' object", call.=FALSE)
  } else if (is.null(env$control)) {
  }
  if (!exists("table", envir=env)) {
    stop("need 'table' object", call.=FALSE)
  } else if (is.null(env$table)) {
  }
  UseMethod("nlmixr2Est")
}
#' Call nlmixr2Est wrapped to collect the warnings
#'
#'
#' @param env nlmixr2 estimate call
#' @param ... Other parameters
#' @return nlmixr2 object
#' @author Matthew L. Fidler
#' @noRd
nlmixr2Est0 <- function(env, ...) {
  .lst <- .collectWarnings(nlmixr2Est(env, ...), lst = TRUE)
  .ret <- .lst[[1]]
  assign("warnings", .lst[[2]], .ret$env)
  lapply(.lst[[2]], warning, call.=FALSE)
  .ret
}
