#' Generic for nlmixr2 estimation methods
#'
#' @param env Environment for the nlmixr2 estimation routines.
#'
#' This needs to have:
#'
#' - rxode2 ui object in `$ui`
#'
#' - data to fit in the estimation routine in `$data`
#'
#' - control for the estimation routine's control options in `$ui`
#'
#' @param ... Other arguments provided to `nlmixr2Est()` provided for
#'   flexibility but not currently used inside nlmixr
#'
#' @return nlmixr2 fit object
#'
#' @author Matthew Fidler
#'
#' @details
#'
#' This is a S3 generic that allows others to use the nlmixr2
#'   environment to do their own estimation routines
#'
#' @export
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
.tablePassthrough <- c("addDosing", "subsetNonmem", "cores", "keep", "drop")

#' Call nlmixr2Est wrapped to collect the warnings
#'
#'
#' @param env nlmixr2 estimate call
#' @param ... Other parameters
#' @return nlmixr2 object
#' @author Matthew L. Fidler
#' @noRd
nlmixr2Est0 <- function(env, ...) {
  rxode2::rxUnloadAll()
  if (!exists("missingTable", envir=env)) {
    env$missingTable <- FALSE
  }
  if (!exists("missingControl", envir=env)) {
    env$missingControl <- FALSE
  }
  if (!exists("missingEst", envir=env)) {
    env$missingEst <- FALSE
  }
  if (env$missingTable) {
    .meta <- env$ui$meta
    if (is.null(env$table)) {
      env$table <- tableControl()
    }
    .table <- env$table
    for (.elt in .tablePassthrough) {
      if (exists(.elt, envir=.meta)) {
        .table[[.elt]] <- .meta[[.elt]]
      }
    }
    env$table <- .table
  }
  .envReset <- new.env(parent=emptyenv())
  if (!getOption("nlmixr2.resetCache", TRUE)) {
    .envReset$ret <- .collectWarnings(nlmixr2Est(env, ...), lst = TRUE)
  } else{
    .envReset$reset <- TRUE
    .envReset$env <- new.env(parent=emptyenv())
    lapply(ls(envir = env, all.names = TRUE), function(item) {
      assign(item, get(item, envir = env), envir = .envReset$env)
    })
    .envReset$cacheReset <- FALSE
    .envReset$unload <- FALSE
    while (.envReset$reset) {
      .envReset$reset <- FALSE
      .envReset$ret <-try(.collectWarnings(nlmixr2Est(env, ...), lst = TRUE))
      if (inherits(.envReset$ret, "try-error")) {
        .msg <- attr(.envReset$ret, "condition")$message
        if (regexpr("not provided by package", .msg) != -1) {
          if (.envReset$cacheReset) {
            .malert("unsuccessful cache reset; try manual reset with 'rxode2::rxClean()'")
            stop(.msg, call.=FALSE)
          } else {
            # reset
            rm(list=ls(envir = env, all.names = TRUE), envir=env)
            lapply(ls(envir = .envReset$env, all.names = TRUE), function(item) {
              assign(item, get(item, envir = .envReset$env), envir = env)
            })
            gc()
            .minfo("try resetting cache")
            rxode2::rxClean()
            .envReset$cacheReset <- TRUE
            .envReset$reset <- TRUE
            .msuccess("done")
          }
        } else if (regexpr("maximal number of DLLs reached", .msg) != -1) {
          if (.envReset$unload) {
            .malert("Could not unload rxode2 models, try restarting R")
            stop(.msg, call.=FALSE)
          } else {
            # reset
            rm(list=ls(envir = env, all.names = TRUE), envir=env)
            lapply(ls(envir = .envReset$env, all.names = TRUE), function(item) {
              assign(item, get(item, envir = .envReset$env), envir = env)
            })
            gc()
            .minfo("try resetting cache and unloading all rxode2 models")
            try(rxode2::rxUnloadAll())
            rxode2::rxClean()
            .envReset$unload <- TRUE
            .envReset$reset <- TRUE
            .msuccess("done")
          }
        } else {
          stop(.msg, call.=FALSE)
        }
      }
    }
  }
  .lst <- .envReset$ret
  .ret <- .lst[[1]]
  assign("warnings", .lst[[2]], .ret$env)
  lapply(.lst[[2]], warning, call.=FALSE)
  .ret
}
