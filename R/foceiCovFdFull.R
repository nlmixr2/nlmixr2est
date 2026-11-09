# covFull=TRUE finite-difference full theta+sigma+Omega covariance: C++
# (foceiCalcRFdFull) stashes the Hessian inverse `.fdFullCov` (Rinv_full) and the OPG
# cross-product `.fdFullS` (Sfull); these helpers enumerate the parameters and install
# the cov per covMethod (the FD counterpart to .foceiInstallAnalyticCov).

#' Parameter enumeration for the C++ FD-full covariance (foceiCalcRFdFull).
#'
#' From the cov-step environment `e`, the free structural+residual theta positions
#' (0-based) and free `Omega` lower-triangle elements (1-based), plus matching
#' `om.<theta>` / `cov.<theta>.<theta>` names.  `same()` copies are not parameters:
#' `cpA`/`cpB`/`cpOf` list them with the element they move with.  `NULL` if there is
#' nothing to do.
#' @param e focei cov-step environment
#' @noRd
.foceiFdFullParams <- function(e) {
  ui <- get("ui", e)
  # bounded-parameter transforms put thetas on an internal scale the theta-sized
  # Jacobian hook cannot correct for a full cov; bow out and keep the native cov.
  if (!is.null(ui$boundedTransforms) && length(ui$boundedTransforms) > 0L) {
    return(NULL)
  }
  ini <- ui$iniDf
  isFix <- if (is.null(ini$fix)) rep(FALSE, nrow(ini)) else ini$fix
  isFix[is.na(isFix)] <- FALSE
  thRows <- which(!is.na(ini$ntheta) & !isFix) # structural + residual thetas
  if (length(thRows) == 0L) {
    return(NULL)
  }
  thPos <- as.integer(ini$ntheta[thRows] - 1L) # 0-based position in fullTheta
  thNames <- ini$name[thRows]
  Om <- get("omega", e)
  pairs <- .omegaDeclaredPairs(Om, ini) # Omega lower-triangle (a>=b)
  if (is.null(pairs)) {
    pairs <- matrix(integer(0), 0L, 2L)
  }
  map <- .foceiEtaThetaMap(ui)
  tie <- .foceiFdFullSameTie(
    pairs,
    .omegaFixed(ini, pairs),
    ui$omegaSameMap,
    map$etaNames,
    ini$name[!is.na(ini$ntheta)]
  )
  pairs <- pairs[tie$free, , drop = FALSE]
  if (nrow(pairs) == 0L) {
    return(list(thPos = thPos, omA = integer(0), omB = integer(0), names = thNames))
  }
  list(
    thPos = thPos,
    omA = as.integer(pairs[, 1]),
    omB = as.integer(pairs[, 2]),
    cpA = tie$cpA,
    cpB = tie$cpB,
    cpOf = tie$cpOf,
    names = c(thNames, .foceiOmegaCovNames(pairs, tie$onm))
  )
}

#' Tie `same()` Omega copies to the element they repeat
#'
#' A `same()` block (an IOV occasion's copy) shares its master's parameters, so it
#' gets no direction of its own: the FD stage moves each copy with its master.
#' @param pairs declared Omega element pairs (`.omegaDeclaredPairs()`)
#' @param fixed which `pairs` are fixed
#' @param same `omegaSameMap` (each eta's master eta, 0 for none), or `NULL`
#' @param onm eta names
#' @param thetas the model's theta names, to name an IOV group
#' @return list: `free` (rows of `pairs` that are parameters), `cpA`/`cpB` (the
#'   copies of a free element) with `cpOf` (0-based index of their parameter
#'   among the free pairs), and `onm` with each master named by its IOV parameter
#' @noRd
.foceiFdFullSameTie <- function(pairs, fixed, same, onm, thetas = character(0)) {
  .none <- list(free = !fixed, cpA = integer(0), cpB = integer(0), cpOf = integer(0), onm = onm)
  same <- as.integer(same)
  same[is.na(same)] <- 0L
  if (length(same) != length(onm) || !any(same > 0L) || any(same < 0L | same > length(same))) {
    return(.none)
  }
  # follow a chained map to the block every copy repeats
  for (.i in seq_along(same)) {
    .c <- same > 0L & same[pmax(same, 1L)] > 0L
    if (!any(.c)) {
      break
    }
    same[.c] <- same[same[.c]]
  }
  if (any(same > 0L & same[pmax(same, 1L)] > 0L)) {
    return(.none) # a cycle
  }
  .ma <- same[pairs[, 1]]
  .mb <- same[pairs[, 2]]
  .isCp <- .ma > 0L & .mb > 0L
  if (any(.isCp & (pairs[, 1] == pairs[, 2]) != (.ma == .mb))) {
    return(.none) # two copies of one eta inside a block
  }
  .free <- !.isCp & !fixed
  .key <- function(a, b) paste(pmax(a, b), pmin(a, b))
  .of <- match(.key(.ma, .mb), .key(pairs[.free, 1], pairs[.free, 2]))
  # a copy moves with a free master; a fixed master's copy stays put
  .mv <- .isCp & !is.na(.of)
  for (.m in unique(same[same > 0L])) {
    onm[.m] <- .foceiSameGroupName(onm[.m], onm[same == .m], thetas)
  }
  list(
    free = .free,
    cpA = as.integer(pairs[.mv, 1]),
    cpB = as.integer(pairs[.mv, 2]),
    cpOf = as.integer(.of[.mv] - 1L),
    onm = onm
  )
}

#' The IOV parameter a `same()` group repeats
#'
#' Occasion copies are named `rx.<iov>.<level>`, and the expanded model keeps
#' `<iov>` as its magnitude theta; the group takes the longest such theta every
#' copy is named after.  Any other `same()` group keeps its master's name.
#' @param master master eta name
#' @param copies its copies' eta names
#' @param thetas the model's theta names
#' @return name
#' @noRd
.foceiSameGroupName <- function(master, copies, thetas) {
  .nm <- c(master, copies)
  .hit <- vapply(thetas, function(v) all(startsWith(.nm, paste0("rx.", v, "."))), logical(1))
  if (!any(.hit)) {
    return(master)
  }
  .v <- thetas[.hit]
  .v[which.max(nchar(.v))]
}

#' The three FD-full covariances, each checked by `.covGuard()`
#'
#' The sandwich is usable only when `Rinv` and `solve(S)` both are: built from
#' an indefinite `Rinv` it is still positive semi-definite, so it can pass the
#' check on its own.
#' @param Rinv full Hessian inverse (`.fdFullCov`)
#' @param S full score cross-product (`.fdFullS`), or `NULL`
#' @return named list ("r", "s", "r,s") of `.covGuard()` results, the matrices
#'   dimnamed like `Rinv`
#' @noRd
.foceiFdFullShapes <- function(Rinv, S) {
  .dn <- dimnames(Rinv)
  .shapes <- list(
    r = Rinv,
    s = if (is.matrix(S)) tryCatch(solve(S), error = function(e) NULL),
    "r,s" = if (is.matrix(S)) Rinv %*% S %*% Rinv
  )
  for (.n in names(.shapes)) {
    if (is.matrix(.shapes[[.n]])) {
      dimnames(.shapes[[.n]]) <- .dn
    }
    .shapes[[.n]] <- .covGuard(.shapes[[.n]])
  }
  if (!.shapes$r$ok || !.shapes$s$ok) {
    .shapes[["r,s"]] <- list(ok = FALSE, reason = paste0("needs a positive-definite ", if (.shapes$r$ok) "S" else "R"))
  }
  .shapes
}

#' Install the C++ FD-full covariance as `fit$cov` (and `fit$covR/covS/covRS`) when
#' `covFull = TRUE`, routing on the requested covMethod: "r,s" -> the sandwich
#' `Rinv %*% S %*% Rinv`, "s" -> `solve(S)`, "r" -> `Rinv`; an "r,s" with an
#' indefinite `Rinv` installs `solve(S)` with a warning.  The native cov is
#' kept -- with a warning when the requested shape is not usable, silently when
#' covMethod is not an FD method or the pieces were not computed.  Every usable
#' shape, native or full, is cached for `setCov()`; an unusable one is neither
#' stored nor cached.  FD counterpart to [.foceiInstallAnalyticCov].
#' @param .ret focei fit environment
#' @return invisibly TRUE when the full covariance was installed
#' @noRd
.foceiInstallFdFullCov <- function(.ret) {
  if (!exists(".fdFullCov", envir = .ret, inherits = FALSE)) {
    return(invisible(FALSE))
  }
  # The fit env's covMethod records what the NATIVE theta-only step produced, not what was
  # asked for -- C++ downgrades it to "s" when that step's "r" fails.  The full R/S pieces
  # used below are computed independently of that step, so route on the REQUESTED control;
  # otherwise a requested "r,s" silently installs solve(S) with a usable .Rinv in hand.
  .env <- if (.covIsName(.ret$covMethod)) .ret$covMethod else ""
  .cty <- tryCatch(rxode2::rxGetControl(.ret$ui, "covType", "fd"), error = function(e) "fd")
  .req <- tryCatch(rxode2::rxGetControl(.ret$ui, "covMethod", NA_integer_), error = function(e) NA_integer_)
  .req <- if (identical(.cty, "analytic")) "" else .covMethodFromSlot(.req)
  .type <- .covFdType(if (nzchar(.req)) .req else .env)
  if (!nzchar(.type)) {
    return(invisible(FALSE))
  } # analytic / failed / "" / boundary -> keep native
  .S <- get0(".fdFullS", envir = .ret, inherits = FALSE)
  .full <- .foceiFdFullShapes(get(".fdFullCov", envir = .ret), .S)
  if (.type != "r" && is.null(.S) && .full$r$ok) {
    return(invisible(FALSE))
  } # no S computed -> keep native
  # an indefinite R is not a minimum; fall back to S like the native step does (#1152)
  if (.type == "r,s" && !.full$r$ok && .full$s$ok) {
    .covWarn(sprintf("\"%s\" needs a positive-definite R; using \"%s\"", .covFullName("r,s"), .covFullName("s")))
    .type <- "s"
  }
  # the native theta-only pieces, cached so setCov() can swap to that shape without
  # recomputing anything (they are already in hand)
  .nat <- stats::setNames(mget(c("covR", "covS", "covRS"), envir = .ret, ifnotfound = list(NULL)), c("r", "s", "r,s"))
  .installed <- .full[[.type]]$ok
  if (.installed) {
    # covMethod="s"/"r" write only e["cov"] -- the chosen covariance is not always
    # mirrored into covR/covS/covRS -- so cache the native about to be replaced too
    .envType <- .covFdType(.env)
    if (nzchar(.envType) && is.null(.nat[[.envType]])) {
      .nat[[.envType]] <- .ret$cov
    }
    # Keep the reported covMethod consistent with what was installed: routing on the
    # requested control can install a sandwich where the env still says "s".  Only
    # rewrite the TYPE when it differs, so the env's "r+"/"|r|" decorations survive when
    # they agree; either way the name carries the " (full)" scope suffix.
    .covInstall(
      .ret,
      .full[[.type]]$cov,
      .covFullName(if (identical(.type, .envType)) .env else .type),
      stash = FALSE,
      refresh = "none"
    )
    for (.n in names(.full)) {
      .slot <- c(r = "covR", s = "covS", "r,s" = "covRS")[[.n]]
      if (.full[[.n]]$ok) {
        assign(.slot, .full[[.n]]$cov, envir = .ret)
      } else if (.n != "s" && !.full$r$ok && exists(.slot, envir = .ret, inherits = FALSE)) {
        # the native theta-only piece is cached below; do not leave it beside a full cov
        rm(list = .slot, envir = .ret)
      }
    }
  } else {
    .covRejectWarn(.ret, .covFullName(.type), .full[[.type]]$reason)
  }
  for (.n in names(.nat)) {
    if (.covGuard(.nat[[.n]])$ok) .covCacheAdd(.ret, .n, .nat[[.n]])
  }
  for (.n in names(.full)) {
    if (.full[[.n]]$ok) .covCacheAdd(.ret, .covFullName(.n), .full[[.n]]$cov)
  }
  .covCacheDrop(.ret, .ret$covMethod)
  .covCacheDrop(.ret, .covFullName(.type))
  # Report the swap: the SEs the C++ step derived from the native theta-only
  # covariance describe a matrix that is no longer $cov, so the caller must
  # refresh the parameter table.
  invisible(.installed)
}
