## Estimation interceptors: let a downstream package claim an estimation before
## the standard method runs, based on the model/environment.  Used by nlmixr2nn
## so that fitting a model that contains an nn() term with a standard method
## (est="focei", "saem", ...) transparently trains the embedded neural network.

.nlmixr2EstInterceptors <- new.env(parent = emptyenv())
## re-entrancy guard: while a claimed interceptor runs, interceptors are
## suppressed so its own inner nlmixr2() calls dispatch to the standard methods.
.nlmixr2EstInterceptState <- new.env(parent = emptyenv())
.nlmixr2EstInterceptState$depth <- 0L

#' Register / remove an estimation interceptor
#'
#' An interceptor is a `function(env)` consulted in [nlmixr2Est()] before the
#' standard estimation method is dispatched.  It receives the fully set-up
#' estimation environment (`env$ui`, `env$data`, `env$control`, the est name in
#' `class(env)[1]`, and any extra `nlmixr2()` arguments in `env$.nlmixr2Dots`).
#' If it returns a non-`NULL` value that value becomes the fit -- the interceptor
#' has *claimed* the estimation; if it returns `NULL` it declines and the next
#' interceptor (or the standard method) runs.  While a claimed interceptor runs,
#' interceptors are suppressed, so it may call `nlmixr2()` internally and get the
#' ordinary methods.
#'
#' @param name unique interceptor name.
#' @param fun `function(env)` returning a finalized fit to claim, or `NULL`.
#' @return invisibly the function (register) or `TRUE`/`FALSE` (remove).
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
registerEstInterceptor <- function(name, fun) {
  checkmate::assertCharacter(name, len = 1, min.chars = 1)
  checkmate::assertFunction(fun, args = "env")
  assign(name, fun, envir = .nlmixr2EstInterceptors)
  invisible(fun)
}

#' @rdname registerEstInterceptor
#' @export
removeEstInterceptor <- function(name) {
  if (exists(name, envir = .nlmixr2EstInterceptors, inherits = FALSE)) {
    rm(list = name, envir = .nlmixr2EstInterceptors)
    invisible(TRUE)
  } else {
    invisible(FALSE)
  }
}

## Run the registered interceptors on `env`; returns a claimed fit or NULL.
.nlmixr2RunEstInterceptors <- function(env) {
  if (.nlmixr2EstInterceptState$depth > 0L) return(NULL)         # nested: skip
  .fns <- ls(envir = .nlmixr2EstInterceptors, all.names = TRUE)
  if (length(.fns) == 0L) return(NULL)
  for (.n in .fns) {
    .fn <- get(.n, envir = .nlmixr2EstInterceptors, inherits = FALSE)
    .nlmixr2EstInterceptState$depth <- .nlmixr2EstInterceptState$depth + 1L
    .ret <- tryCatch(.fn(env),
                     finally = {
                       .nlmixr2EstInterceptState$depth <-
                         .nlmixr2EstInterceptState$depth - 1L
                     })
    if (!is.null(.ret)) return(.ret)                             # claimed
  }
  NULL
}
