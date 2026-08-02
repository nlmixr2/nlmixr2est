#' Compile a generated model under a stable, role-tagged artifact name.
#'
#' rxode2 names a model's generated `.c`/`.so` from `.rxPre(model, modName)`, which for
#' an ANONYMOUS model is `rx_<parsed_md5>_<arch>_` -- the parsed model text alone.  The
#' emitted C also depends on inputs the parsed text cannot see (most importantly the
#' event-sensitivity code, generated afterwards and injected), so two builds of one
#' model text that differ in those inputs land on a single `.so`.  The later build wins
#' for both, and because entry points resolve BY NAME (`R_GetCCallable`) a model object
#' bound to the earlier one silently starts executing the replacement.
#'
#' Measured: an augmented FOCEi sensitivity model built at 193856 bytes was replaced by
#' a 167792-byte build under the same name, after which its `calc_lhs` computed 4 of the
#' 29 columns it declares.  Every analytic outer gradient then came back non-finite and
#' the fit silently finished on finite differences.  The same failure mode is what the
#' "not provided by package" cache-reset recovery in `nlmixr2Est.R` papers over.
#' See nlmixr2/rxode2#1171.
#'
#' Tagging every generated model with its ROLE keeps nlmixr2est's own models from
#' colliding with each other or with an anonymous build of the same text.  The parsed
#' md5 stays in the name, so two genuinely different models still never share an
#' artifact -- a bare role name would be worse than the default, not better.
#'
#' @param model model text (or anything `rxode2::rxModelVars()` accepts)
#' @param role short role tag identifying which generated model this is, e.g.
#'   `"rxInner"`, `"rxOuter"`, `"rxPred"`.  Distinct roles must not share a name.
#' @param ... passed to `rxode2::rxode2()` (e.g. `eventSens`)
#' @return the compiled rxode2 model
#' @noRd
.nlmixr2estRxode2 <- function(model, role, ...) {
  .nm <- .nlmixr2estModName(model, role, ...)
  .wd <- .nlmixr2estModDir()
  ## The fallback must still build in OUR directory.  Dropping back to rxode2's own
  ## naming AND its own directory would put the model right back where two builds of
  ## one text overwrite each other, which is the failure this exists to prevent.
  if (is.null(.nm)) return(rxode2::rxode2(model, wd = .wd, ...))
  rxode2::rxode2(model, modName = .nm, wd = .wd, ...)
}

#' Stable artifact name for a generated model: role + everything that changes its code.
#'
#' `eventSens` is folded in as well as the role, since it selects which
#' event-sensitivity code rxode2 emits and would otherwise let a "jump" and an "fd"
#' build of one model text share an artifact.
#' @inheritParams .nlmixr2estRxode2
#' @return modName string, or NULL to keep rxode2's default naming
#' @noRd
.nlmixr2estModName <- function(model, role, ...) {
  if (is.null(role) || !nzchar(role)) return(NULL)
  .md5 <- tryCatch(rxode2::rxModelVars(model)$md5[["parsed_md5"]], error = function(e) NULL)
  if (is.null(.md5) || !nzchar(.md5)) return(NULL)
  .dots <- list(...)
  ## eventSens can arrive as an un-evaluated match.arg default, i.e. c("jump", "fd"):
  ## take the first element, as match.arg would, so the name stays length 1.  A
  ## zero-length value must fall back to "" rather than produce character(0), which
  ## would sail past an is.null() check and reach rxode2 as an empty modName.
  .es <- .dots$eventSens
  .es <- if (is.null(.es) || length(.es) == 0L) "" else gsub("\\W", "", as.character(.es)[1L])
  if (is.na(.es)) .es <- ""
  ## "_" separates the parts: without it role "rxA" + es "bc" and role "rxAb" + es "c"
  ## would produce one name.  The md5 is last and fixed-width, so distinct model texts
  ## can never collide.
  .nm <- paste0(role, "_", .es, "_", .md5)
  if (length(.nm) != 1L || is.na(.nm) || !nzchar(.nm)) return(NULL)
  .nm
}

#' Where generated models are built.
#'
#' Three constraints, and only one directory satisfies all of them.
#'
#' NOT `getwd()`: a NAMED model builds there by default, which would scatter
#' `<modName>.d` directories through the user's working directory.
#'
#' NOT `rxTempDir()`: that directory is rxode2's, and something in a session
#' invalidates artifacts inside it -- measured as a generated model's `.so` being
#' replaced mid-run by a build emitting different event-sensitivity code (193856 bytes
#' -> 167792, after which its `calc_lhs` computed 4 of 29 declared lhs).  Role-tagged
#' names do NOT protect against this, because the replacement arrives under the same
#' name; it was still reproducible with every generated model uniquely named.
#'
#' R's session temp directory is neither: nothing else clears it, no consent is needed
#' to use it, and it goes away when R exits.  The cost is that artifacts are not shared
#' with an opted-in `rxCreateCache()`, so each session rebuilds them once -- worth it
#' while nlmixr2/rxode2#1171 is open, since the alternative is silently wrong
#' gradients.  Revisit when that is fixed.
#' @return build directory (created on demand)
#' @noRd
.nlmixr2estModDir <- function() {
  .d <- file.path(tempdir(), "nlmixr2estSens")
  if (!dir.exists(.d)) {
    dir.create(.d, recursive = TRUE, showWarnings = FALSE)
    ## dir.create() is quiet about failure here (another process may have won the
    ## race, which is fine), so confirm rather than hand back a path that does not
    ## exist -- rxode2 would then fail deep in the compile with a confusing error.
    if (!dir.exists(.d)) return(tempdir())
  }
  .d
}
