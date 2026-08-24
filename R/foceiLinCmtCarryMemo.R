# Session memo + persistent cache-directory sidecar for the linCmt() carry
# candidate pairs.  Layers (fastest first): session env -> the sidecar file
# in rxode2's cache directory (the same directory and digest key
# rxUiGet.foceiModel already uses for the built model shapes, so a warm
# cache reloads instead of recomputing) -> full symbolic detection.  Only
# the data-INDEPENDENT candidate result is stored; every data-dependent
# conclusion is re-derived per fit (foceiLinCmtCarryData.R).

#' Session memo for the data-independent candidate pairs, keyed by the
#' focei model digest (which already covers the model text, iniDf,
#' linCmtSensCarry, the rxode2 carry capability and covsInterpolation --
#' everything the symbolic detection depends on).  Only the `data = NULL`
#' candidate result is memoized; every data-dependent conclusion (varying,
#' ss/evid, jump-regimen checks) is re-derived per fit.
#' @noRd
.carryPairsMemo <- new.env(parent = emptyenv())

#' Drop every memoized session entry, preserving the counters; with
#' `files = TRUE` also remove those entries' sidecar files.  `files =
#' FALSE` simulates a fresh session against a warm cache directory.
#' (tests; a session never needs this itself)
#' @noRd
.foceiLinCmtCarryMemoClear <- function(files = TRUE) {
  .keys <- setdiff(ls(envir = .carryPairsMemo, all.names = FALSE),
                   c("hits", "misses", "fileHits"))
  rm(list = .keys, envir = .carryPairsMemo)
  if (files) {
    .dir <- tryCatch(rxode2::rxTempDir(), error = function(e) NULL)
    if (!is.null(.dir)) {
      unlink(file.path(.dir, paste0("focei-carry-", .keys, ".rds")))
    }
  }
  invisible(NULL)
}

#' Memo hit/miss counters (mechanism observable for tests); fileHits counts
#' results served from the persistent sidecar in rxode2's cache directory
#' @noRd
.foceiLinCmtCarryMemoStats <- function(reset = FALSE) {
  .st <- c(hits = get0("hits", envir = .carryPairsMemo, ifnotfound = 0L),
           misses = get0("misses", envir = .carryPairsMemo, ifnotfound = 0L),
           fileHits = get0("fileHits", envir = .carryPairsMemo, ifnotfound = 0L))
  if (reset) {
    assign("hits", 0L, envir = .carryPairsMemo)
    assign("misses", 0L, envir = .carryPairsMemo)
    assign("fileHits", 0L, envir = .carryPairsMemo)
  }
  .st
}

#' Sidecar file for the candidate pairs in rxode2's cache directory,
#' following the focei model cache convention (focei-<digest>.rds there).
#' With no persistent cache (rxCreateCache() never run) rxTempDir() is a
#' session temp dir, so the sidecar quietly degrades to session scope.
#' @noRd
.foceiLinCmtCarryCacheFile <- function(key) {
  tryCatch(file.path(rxode2::rxTempDir(), paste0("focei-carry-", key, ".rds")),
           error = function(e) NULL)
}

#' Read the sidecar; NULL unless it exists and was written by THIS
#' nlmixr2est build (nlmixr2.md5 -- a package upgrade invalidates; an
#' rxode2 upgrade wipes the whole cache directory itself)
#' @noRd
.foceiLinCmtCarryCacheRead <- function(key) {
  .f <- .foceiLinCmtCarryCacheFile(key)
  if (is.null(.f) || !file.exists(.f)) return(NULL)
  .x <- tryCatch(readRDS(.f), error = function(e) NULL)
  if (!is.list(.x) || !identical(.x$md5, nlmixr2.md5)) return(NULL)
  .x$pairs
}

#' Best-effort sidecar write (a read-only cache directory is not an error)
#' @noRd
.foceiLinCmtCarryCacheWrite <- function(key, pairs) {
  .f <- .foceiLinCmtCarryCacheFile(key)
  if (is.null(.f)) return(invisible(NULL))
  tryCatch(saveRDS(list(md5 = nlmixr2.md5, pairs = pairs), .f),
           error = function(e) NULL)
  invisible(NULL)
}

#' Memo/sidecar lookup: session env first, then the cache-directory
#' sidecar (hydrating the session env on a file hit); NULL on a miss
#' @noRd
.foceiLinCmtCarryMemoGet <- function(key) {
  if (is.null(key)) return(NULL)
  .bump <- function(ctr) {
    assign(ctr, .foceiLinCmtCarryMemoStats()[[ctr]] + 1L, envir = .carryPairsMemo)
  }
  if (exists(key, envir = .carryPairsMemo, inherits = FALSE)) {
    .bump("hits")
    return(get(key, envir = .carryPairsMemo, inherits = FALSE))
  }
  .cached <- .foceiLinCmtCarryCacheRead(key)
  if (is.null(.cached)) return(NULL)
  .bump("fileHits")
  assign(key, .cached, envir = .carryPairsMemo)
  .cached
}

#' Store a computed result in the session env (bounded) and the sidecar
#' @noRd
.foceiLinCmtCarryMemoPut <- function(key, pairs) {
  if (is.null(key)) return(invisible(NULL))
  assign("misses", .foceiLinCmtCarryMemoStats()[["misses"]] + 1L,
         envir = .carryPairsMemo)
  # bound the memo (a session rarely fits more than a handful of models)
  .keys <- setdiff(ls(envir = .carryPairsMemo, all.names = FALSE),
                   c("hits", "misses", "fileHits"))
  if (length(.keys) >= 64L) {
    rm(list = .keys, envir = .carryPairsMemo)
  }
  assign(key, pairs, envir = .carryPairsMemo)
  .foceiLinCmtCarryCacheWrite(key, pairs)
  invisible(NULL)
}
