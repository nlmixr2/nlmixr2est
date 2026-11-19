# Precursors of the finite-difference covariance (covPrecursor): what a fit already
# holds about the curvature at the estimates.  Its diagonal seeds the full stage's
# step searches, and the verified shortcut (covShortcut) may skip measuring the
# off-diagonals when its correlations explain the curvature along fixed directions.

#' The precursor sources a covariance step can read
#'
#' `"fd"`: a full-stage finite-difference R this fit's covariance store holds for
#' other settings; `"analytic"`: the analytic covariance this fit holds.
#' @noRd
.covPrecursorSources <- c("fd", "analytic")

#' Check a `covPrecursor` value
#'
#' @param covPrecursor `NULL` (no precursors), or the sources to read, in order of
#'   preference
#' @return the value, `character(0)` for none
#' @noRd
.covPrecursorCheck <- function(covPrecursor) {
  if (is.null(covPrecursor) || (is.character(covPrecursor) && length(covPrecursor) == 0L)) {
    return(character(0))
  }
  .allowed <- paste0("\"", .covPrecursorSources, "\"", collapse = ", ")
  if (!is.character(covPrecursor) || anyNA(covPrecursor)) {
    stop("'covPrecursor' must be NULL or some of ", .allowed, call. = FALSE)
  }
  .bad <- setdiff(covPrecursor, .covPrecursorSources)
  if (length(.bad) > 0L) {
    stop(
      "'covPrecursor' has unknown source(s) ",
      paste0("\"", .bad, "\"", collapse = ", "),
      "; use some of ",
      .allowed,
      ", or NULL for none",
      call. = FALSE
    )
  }
  if (anyDuplicated(covPrecursor)) {
    stop("'covPrecursor' names a source more than once", call. = FALSE)
  }
  covPrecursor
}

#' The hint for a refit's full finite-difference stage
#'
#' The first source in `control$covPrecursor` that the fit holds, as the R matrix
#' (information) over the full stage's parameters.  The C++ stage checks that it
#' fits (names, finite) and records how it served.
#' @param env fit environment
#' @param control the refit's `foceiControl`
#' @param key the refit's covariance-store key, or `NULL`
#' @return list(R, source, shortcut), or `NULL` when no source is held
#' @noRd
.covPrecursorHint <- function(env, control, key) {
  .src <- control$covPrecursor
  if (length(.src) == 0L) {
    return(NULL)
  }
  .nm <- tryCatch(.foceiFdFullParams(env)$names, error = function(e) NULL)
  if (length(.nm) == 0L) {
    return(NULL)
  }
  for (.s in .src) {
    .R <- switch(
      .s,
      fd = .covPrecursorFd(env, key, .nm),
      analytic = .covPrecursorAnalytic(env, .nm)
    )
    if (is.matrix(.R)) {
      return(list(R = .R, source = .s, shortcut = isTRUE(control$covShortcut)))
    }
  }
  NULL
}

#' A full-stage R the fit's store holds for other settings at the same estimates
#'
#' The most recent entry with the key's hand-off and versions, under different
#' settings.
#' @param env fit environment
#' @param key the refit's covariance-store key
#' @param nm the full stage's parameter names
#' @return the R matrix, or `NULL`
#' @noRd
.covPrecursorFd <- function(env, key, nm) {
  .store <- env$covStore
  if (is.null(key) || !is.list(.store)) {
    return(NULL)
  }
  for (.e in rev(.store)) {
    .R <- .e$full$R
    if (
      identical(.e$key$handoff, key$handoff) &&
        identical(.e$key$versions, key$versions) &&
        !identical(.e$key$settings, key$settings) &&
        is.matrix(.R) &&
        identical(rownames(.R), nm)
    ) {
      return(.R)
    }
  }
  NULL
}

#' The information of the analytic covariance the fit holds
#'
#' The installed covariance or a `covList` entry labelled `"analytic"`, over the full
#' stage's parameters.
#' @param env fit environment
#' @param nm the full stage's parameter names
#' @return `solve(cov)` with the parameter names, or `NULL`
#' @noRd
.covPrecursorAnalytic <- function(env, nm) {
  .cands <- .covCacheGet(env)
  if (is.matrix(env$cov) && checkmate::testString(env$covMethod)) {
    .cands <- c(stats::setNames(list(env$cov), env$covMethod), .cands)
  }
  for (.n in names(.cands)) {
    .cov <- .cands[[.n]]
    if (.covBaseName(.n) != "analytic" || !is.matrix(.cov) || !identical(rownames(.cov), nm)) {
      next
    }
    .R <- tryCatch(solve(.cov), error = function(e) NULL)
    if (is.matrix(.R) && all(is.finite(.R))) {
      dimnames(.R) <- list(nm, nm)
      return(.R)
    }
  }
  NULL
}

#' Keep how a precursor served a covariance step
#'
#' `env$covPrecursorUsed[[label]]` is list(source, shortcut, checks) from the
#' refit's full stage (`.fdFullPrecursor`); a label computed without one has no entry.
#' @param env fit environment
#' @param label installed covariance label
#' @param src environment the covariance step ran in
#' @return invisibly `TRUE`
#' @noRd
.covPrecursorRecord <- function(env, label, src) {
  .rec <- get0(".fdFullPrecursor", envir = src, inherits = FALSE)
  .all <- env$covPrecursorUsed
  if (!is.list(.all)) {
    .all <- list()
  }
  .all[[label]] <- if (is.list(.rec)) .rec
  assign("covPrecursorUsed", .all, envir = env)
  invisible(TRUE)
}

#' One line describing how a precursor served a covariance
#' @param rec entry of `env$covPrecursorUsed`
#' @return character, or `NULL` when there is no record
#' @noRd
.covPrecursorLine <- function(rec) {
  if (!is.list(rec) || !checkmate::testString(rec$source)) {
    return(NULL)
  }
  .sc <- if (checkmate::testString(rec$shortcut) && rec$shortcut != "off") {
    paste0("; shortcut ", rec$shortcut)
  } else {
    ""
  }
  paste0("from the \"", rec$source, "\" precursor (seeded steps", .sc, ")")
}
