## Cache of lowered simulation models, keyed on the fitted model
##
## `rxode2::rxSolve()` on a fit asks the fit for `$simulationModel` on every
## call, and lowering the fit to an rxode2 model is neither cheap nor free of
## side effects: for an ODE model it is ~50 ms and a couple of MB of process
## memory per call, and the memory is not given back by `gc()`,
## `rxode2::rxUnloadAll()` or `rxode2::rxClean()`.  Simulating from a fit in a
## loop therefore grew without bound (nlmixr2/rxode2#1289).
##
## The cache is keyed on the `ui` object held in the fit environment, and an
## entry holds that `ui` -- never the fit -- so a cached entry cannot keep a
## fit (and its data) alive.  The key is the entire model, structure *and*
## final estimates, so two fits that share a key lower to the same simulation
## model and can safely share the entry.  The lowered model is a pure function
## of the `ui`, and piping/updating a fit builds a new `ui`, which is what
## invalidates the entry.
.simModelCache <- new.env(parent = emptyenv())
.simModelCache$entries <- list()

## Number of lowered models kept.  Each entry is only as large as the model
## itself (a couple hundred KB for a typical PK model), and simulation loops
## work on one fit at a time, so a small cache is enough to remove the
## repeated lowering while bounding what the package holds onto.
.simModelCacheMax <- 5L

#' Key used to cache the lowered simulation model of a fit
#'
#' @param env object handed to `nmObjGet.simulationModel()`, that is the
#'   environment of a nlmixr2 fit.  Anything else (including a bare `rxUi`,
#'   which is an environment without a `ui` binding) is not cached.
#' @return the fit's `ui` object, or `NULL` when the value should not be cached
#' @noRd
.simModelCacheKey <- function(env) {
  if (!isTRUE(getOption("nlmixr2.simModelCache", TRUE))) {
    return(NULL)
  }
  if (!is.environment(env)) {
    return(NULL)
  }
  if (!exists("ui", envir = env, inherits = FALSE)) {
    return(NULL)
  }
  .ui <- get("ui", envir = env, inherits = FALSE)
  ## a fit holds a compressed ui, which has value semantics: replacing the
  ## model replaces the object, and that is what invalidates the entry.  An
  ## uncompressed ui is an environment, which could be edited in place behind
  ## the cache's back, so it is not a key that can be trusted.
  if (is.environment(.ui)) {
    return(NULL)
  }
  .ui
}

#' Get a cached lowered simulation model
#'
#' @param ui key from [.simModelCacheKey()]
#' @return the cached model, or `NULL` when this model has not been lowered yet
#' @noRd
.simModelCacheGet <- function(ui) {
  .entries <- .simModelCache$entries
  for (.i in seq_along(.entries)) {
    if (identical(.entries[[.i]]$ui, ui)) {
      if (.i != 1L) {
        # most recently used first, so the least recently used falls off the end
        .simModelCache$entries <- c(.entries[.i], .entries[-.i])
      }
      return(.entries[[.i]]$model)
    }
  }
  NULL
}

#' Cache a lowered simulation model
#'
#' @param ui key from [.simModelCacheKey()]
#' @param model lowered rxode2 simulation model
#' @return `model`, invisibly
#' @noRd
.simModelCacheSet <- function(ui, model) {
  .entries <- .simModelCache$entries
  .keep <- vapply(seq_along(.entries),
    function(i) !identical(.entries[[i]]$ui, ui),
    logical(1),
    USE.NAMES = FALSE
  )
  .entries <- c(list(list(ui = ui, model = model)), .entries[.keep])
  if (length(.entries) > .simModelCacheMax) {
    .entries <- .entries[seq_len(.simModelCacheMax)]
  }
  .simModelCache$entries <- .entries
  invisible(model)
}

#' Empty the lowered simulation model cache
#'
#' Only needed when measuring the cache itself; the cache is keyed on the
#' fitted model, so a stale entry cannot be returned for a changed model.
#'
#' @return nothing, called for the side effect
#' @noRd
.simModelCacheReset <- function() {
  .simModelCache$entries <- list()
  invisible()
}
