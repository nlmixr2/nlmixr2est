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
##' @export
nlmixr2Est <- function(env, ...) {
  if (!exists(env, "ui")) {
    stop("need 'ui' object", call.=FALSE)
  } else if (!inherits(env$ui, "rxUi")) {
    stop("'ui' is not an rxode2 object", call.=FALSE)
  }
  if (!exists(env, "data")) {
    stop("need 'data' object", call.=FALSE)
  } else if (!inherits(env$data, "data.frame")) {
    stop("'data' is not a data.frame", call.=FALSE)
  }
  env$data <- as.data.frame(env$data)
  if (!exists(env, "control")) {
    stop("need 'control' object", call.=FALSE)
  } else if (is.null(env$control)) {
  } else if (!inherits(env$control, "list")) {
    stop("'control' object should be a list-like object", call.=FALSE)
  }
  if (!exists(env, "table")) {
    stop("need 'table' object", call.=FALSE)
  } else if (is.null(env$table)) {
  } else if (!inherits(env$table, "list")) {
    stop("'table' object should be a list-like object", call.=FALSE)
  }
  UseMethod("nlmixr2Est")
}
