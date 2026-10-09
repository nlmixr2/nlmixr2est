# The fit's covariance store: the finite-difference ingredients each covariance step
# computed (R, S and the steps they were taken at), keyed by the hand-off (the estimates
# the step started from, env$covHandoff) and the settings they depend on.  A later
# setCov() request with the same key reads them instead of computing them again, and
# computes only what is missing (for example the S of "r,s" after "r").

#' Start a fit with an empty covariance store
#'
#' A fit's store holds what its own covariance steps computed, so every fit starts with
#' none (`covStore` is `NULL`), and with none of a previous fit's hand-off or stored
#' ingredients.  A `setCov()` refit hands in its fit's (`.covRefitInputs`, from
#' `.setCovRefit()`), which are unpacked here for the C++ covariance step.
#' @param env the estimation environment
#' @return invisibly `env`
#' @noRd
.covStoreFitStart <- function(env) {
  env$covStore <- NULL
  for (.n in c("covHandoff", ".fdFullStore", "covThetaStore")) {
    if (exists(.n, envir = env, inherits = FALSE)) rm(list = .n, envir = env)
  }
  .in <- get0(".covRefitInputs", envir = env, inherits = FALSE)
  if (is.list(.in)) {
    for (.n in c("covHandoff", ".fdFullStore", "covThetaStore")) {
      if (!is.null(.in[[.n]])) assign(.n, .in[[.n]], envir = env)
    }
    rm(list = ".covRefitInputs", envir = env)
  }
  invisible(env)
}

# Control elements the stored R, S and steps depend on, beyond the estimates: the step
# searches and differences (the covariance stage's, and the estimation's, which a separate
# full stage uses), the derivative method and the covariance step's tolerances.
.covStoreKeyFields <- c(
  "hessEps",
  "hessEpsLlik",
  "gillKcov",
  "gillKcovLlik",
  "gillStepCov",
  "gillStepCovLlik",
  "gillFtolCov",
  "gillFtolCovLlik",
  "covGillF",
  "rmatNorm",
  "rmatNormLlik",
  "smatNorm",
  "smatNormLlik",
  "covDerivMethod",
  "gillK",
  "gillStep",
  "gillFtol",
  "gillRtol",
  "covSolveTol",
  "covInnerTol"
)

#' Key of a covariance step in the fit's covariance store
#'
#' @param env fit environment (its `covHandoff`)
#' @param control the `foceiControl()` the covariance step runs with
#' @return list(handoff, settings), or `NULL` when the fit has no hand-off
#' @noRd
.covStoreKey <- function(env, control) {
  .h <- env$covHandoff
  if (!is.list(.h)) {
    return(NULL)
  }
  .settings <- lapply(stats::setNames(.covStoreKeyFields, .covStoreKeyFields), function(.f) control[[.f]])
  # the inner budget of the covariance legs: 0 holds the ETAs fixed (a conditional
  # covariance); a refit at a fit's estimates re-optimizes them (covMaxInnerIterations)
  .inner <- control$maxInnerIterations
  if (!checkmate::testNumber(.inner, lower = 1)) {
    .inner <- control$covMaxInnerIterations
  }
  .settings$covEtaLegs <- if (checkmate::testNumber(.inner, lower = 1)) as.integer(.inner) else 0L
  # a saved fit reloaded under other versions computes with other numerics
  .versions <- c(
    nlmixr2est = as.character(utils::packageVersion("nlmixr2est")),
    rxode2 = as.character(utils::packageVersion("rxode2"))
  )
  list(handoff = .h, settings = .settings, versions = .versions)
}

#' Whether a refit's settings are all ones the store's key covers
#'
#' A refit given anything else (`getVarCov(force = TRUE, ...)` can pass any control
#' element) neither reads nor writes the store.
#' @param args the refit's control arguments
#' @return single logical
#' @noRd
.covStoreRefitOk <- function(args) {
  all(names(args) %in% c("covMethod", "covFull", "covSmall", "covFallback", .covStoreKeyFields))
}

#' Index of the store entry for a key
#' @param store list of entries
#' @param key from `.covStoreKey()`
#' @return integer index, 0 when there is none
#' @noRd
.covStoreIndex <- function(store, key) {
  for (.i in seq_along(store)) {
    if (identical(store[[.i]]$key, key)) {
      return(.i)
    }
  }
  0L
}

#' Stored ingredients for a key
#' @param env fit environment
#' @param key from `.covStoreKey()`
#' @return the entry (list with `full` and `theta`), or `NULL`
#' @noRd
.covStoreGet <- function(env, key) {
  .store <- env$covStore
  if (is.null(key) || !is.list(.store)) {
    return(NULL)
  }
  .i <- .covStoreIndex(.store, key)
  if (.i == 0L) {
    return(NULL)
  }
  .store[[.i]]
}

#' Add what a covariance step computed to the fit's store
#'
#' The full stage's R, steps, point and S (`.fdFullR`, `.fdFullH`, `.fdFullX0`,
#' `.fdFullS`) and a separate theta-only stage's steps (`covSteps`), finite-difference R
#' (`R.0`) and S (`S0`, `Sper`, `SHasZero`).  An entry already under the key keeps what the new step did not compute, when
#' both were taken at the same steps.
#' @param env fit environment whose store is updated
#' @param key from `.covStoreKey()`
#' @param src environment the covariance step ran in (the fit's own, or a refit's)
#' @param fd whether the theta-only R was a finite-difference one
#' @return invisibly whether an entry was written
#' @noRd
.covStoreRecord <- function(env, key, src = env, fd = TRUE) {
  if (is.null(key)) {
    return(invisible(FALSE))
  }
  .store <- env$covStore
  if (!is.list(.store)) {
    .store <- list()
  }
  .i <- .covStoreIndex(.store, key)
  .e <- if (.i > 0L) .store[[.i]] else list(key = key, full = NULL, theta = NULL)
  .R <- get0(".fdFullR", envir = src, inherits = FALSE)
  .h <- get0(".fdFullH", envir = src, inherits = FALSE)
  .x0 <- get0(".fdFullX0", envir = src, inherits = FALSE)
  if (is.matrix(.R) && is.numeric(.h) && is.numeric(.x0)) {
    .S <- get0(".fdFullS", envir = src, inherits = FALSE)
    .old <- .e$full
    if (
      is.null(.S) &&
        !is.null(.old) &&
        identical(.old$R, .R) &&
        identical(.old$h, .h) &&
        identical(.old$x0, .x0)
    ) {
      .S <- .old$S
    }
    .e$full <- list(R = .R, h = .h, x0 = .x0, S = .S)
  }
  .steps <- get0("covSteps", envir = src, inherits = FALSE)
  if (is.list(.steps)) {
    .old <- .e$theta
    .same <- !is.null(.old) && identical(.old$steps, .steps)
    .R0 <- if (fd) get0("R.0", envir = src, inherits = FALSE)
    .S0 <- get0("S0", envir = src, inherits = FALSE)
    .Sper <- get0("Sper", envir = src, inherits = FALSE)
    .SHasZero <- get0("SHasZero", envir = src, inherits = FALSE)
    if (.same && is.null(.R0)) {
      .R0 <- .old$R0
    }
    if (.same && is.null(.S0)) {
      .S0 <- .old$S0
      .Sper <- .old$Sper
      .SHasZero <- .old$SHasZero
    }
    .e$theta <- list(steps = .steps, R0 = .R0, S0 = .S0, Sper = .Sper, SHasZero = .SHasZero)
  }
  if (is.null(.e$full) && is.null(.e$theta)) {
    return(invisible(FALSE))
  }
  if (.i > 0L) {
    .store[[.i]] <- .e
  } else {
    .store[[length(.store) + 1L]] <- .e
  }
  assign("covStore", .store, envir = env)
  invisible(TRUE)
}
