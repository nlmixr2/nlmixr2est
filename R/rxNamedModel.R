#' Compile a generated model.
#'
#' A thin wrapper over `rxode2::rxode2()`, kept because every generated model in the
#' package goes through one place.
#'
#' It used to do much more.  rxode2 named a model's `.c`/`.so` from the parsed model
#' text alone, while the emitted C also depends on the event-sensitivity code generated
#' afterwards -- so two builds of one text whose event-sensitivity code differed landed
#' on a single `.so`, the later build won for both, and because entry points resolve BY
#' NAME (`R_GetCCallable`) a model object bound to the earlier one silently began
#' executing the replacement.  Measured: an augmented FOCEi sensitivity model built at
#' 193856 bytes was replaced by a 167792-byte build under the same name, after which its
#' `calc_lhs` computed 4 of the 29 columns it declares; every analytic outer gradient
#' then came back non-finite and the fit silently finished on finite differences.
#'
#' This package worked around that by naming every artifact `role_eventSens_md5` and
#' building into `tempdir()/nlmixr2estSens` rather than the rxode2 cache.  rxode2 now
#' folds the event-sensitivity code into its own cache key (nlmixr2/rxode2#1171), so the
#' variants separate upstream and the models belong back in the rxode2 cache directory
#' under rxode2's naming.  Requires rxode2 >= 5.1.6.
#'
#' @param model model text (or anything `rxode2::rxModelVars()` accepts)
#' @param role short role tag, e.g. `"rxInner"`.  No longer shapes the compiled
#'   artifact -- kept so the call sites still read as self-documenting.
#' @param ... passed to `rxode2::rxode2()` (e.g. `eventSens`)
#' @return the compiled rxode2 model
#' @noRd
.nlmixr2estRxode2 <- function(model, role, ...) {
  ## COMPATIBILITY LAYER -- remove with the follow-up issue, once rxode2 with the
  ## eventSensCode cache key (nlmixr2/rxode2#1171) is the minimum version.
  ##
  ## rxode2 < 5.1.6 keyed a model's generated .c/.so on the PARSED MODEL TEXT alone,
  ## while the emitted C also depends on the event-sensitivity code injected afterwards.
  ## Two builds of one text whose event-sensitivity code differed therefore landed on a
  ## single .so, the later build won for both, and because entry points resolve BY NAME
  ## (R_GetCCallable) a model object bound to the earlier one silently began executing
  ## the replacement -- measured as an augmented model declaring 29 lhs whose calc_lhs
  ## computed 4, which drops the analytic gradient for a whole fit.
  ##
  ## Where rxode2 fixes that itself, let it: models go in the rxode2 cache directory
  ## under rxode2's own naming.  Where it does not, fall back to what this package used
  ## to do -- tag the artifact role_eventSens_md5 and build in our own directory.
  if (.nlmixr2estRxHasSensKey()) {
    return(.nlmixr2estStripModelSrcCompat(rxode2::rxode2(model, ...)))
  }
  .nm <- .nlmixr2estModName(model, role, ...)
  .wd <- .nlmixr2estModDir()
  ## The fallback must still build in OUR directory: rxode2's own naming AND its own
  ## directory is exactly where two builds of one text overwrite each other.
  if (is.null(.nm)) {
    return(.nlmixr2estStripModelSrcCompat(rxode2::rxode2(model, wd = .wd, ...)))
  }
  .nlmixr2estStripModelSrcCompat(rxode2::rxode2(model, modName = .nm, wd = .wd, ...))
}

#' Does the installed rxode2 fold eventSensCode into its compiled-model cache key?
#'
#' nlmixr2/rxode2#1171.  Probed once per session.  COMPATIBILITY -- remove with the
#' follow-up issue.
#' @noRd
.nlmixr2estRxHasSensKey <- local({
  .v <- NULL
  function() {
    if (!is.null(.v)) return(.v)
    .f <- tryCatch(get(".rxEventSensKey", envir = asNamespace("rxode2")),
                   error = function(e) NULL)
    .v <<- is.function(.f)
    .v
  }
})

#' Strip the srcrefs that pin a session onto a compiled model.
#'
#' rxode2 >= 5.1.6 does this itself while building the model (rxStripModelSrc), so this
#' only has to act on older rxode2.  COMPATIBILITY -- remove with the follow-up issue.
#' @noRd
.nlmixr2estStripModelSrcCompat <- function(mod) {
  if (!is.environment(mod)) return(mod)
  if (is.function(tryCatch(get("rxStripModelSrc", envir = asNamespace("rxode2")),
                           error = function(e) NULL))) {
    return(mod)                       # rxode2 already stripped it during the build
  }
  .strip <- function(env, nm) {
    if (!is.environment(env)) return(invisible(NULL))
    .f <- tryCatch(get(nm, envir = env, inherits = FALSE), error = function(e) NULL)
    if (is.function(.f)) {
      tryCatch(assign(nm, removeSource(.f), envir = env), error = function(e) NULL)
    }
    invisible(NULL)
  }
  .ap <- tryCatch(get("assignPtr", envir = mod, inherits = FALSE), error = function(e) NULL)
  if (is.function(.ap)) .strip(environment(.ap), ".f")
  .dll <- tryCatch(get(".rxDll", envir = mod, inherits = FALSE), error = function(e) NULL)
  if (is.list(.dll) && is.function(.dll$.call)) .strip(environment(.dll$.call), ".badBuild")
  mod
}

#' Role-tagged artifact name (compatibility fallback; see .nlmixr2estRxode2)
#' @noRd
.nlmixr2estModName <- function(model, role, ...) {
  if (is.null(role) || !nzchar(role)) return(NULL)
  .md5 <- tryCatch(rxode2::rxModelVars(model)$md5[["parsed_md5"]], error = function(e) NULL)
  if (is.null(.md5) || !nzchar(.md5)) return(NULL)
  .dots <- list(...)
  ## eventSens can arrive as an un-evaluated match.arg default, i.e. c("jump", "fd"):
  ## take the first element, as match.arg would, so the name stays length 1.
  .es <- .dots$eventSens
  .es <- if (is.null(.es) || length(.es) == 0L) "" else gsub("\\W", "", as.character(.es)[1L])
  if (is.na(.es)) .es <- ""
  ## "_" separates the parts, and the fixed-width md5 is last, so distinct model texts
  ## can never collide.
  .nm <- paste0(role, "_", .es, "_", .md5)
  if (length(.nm) != 1L || is.na(.nm) || !nzchar(.nm)) return(NULL)
  .nm
}

#' Private build directory (compatibility fallback; see .nlmixr2estRxode2)
#' @noRd
.nlmixr2estModDir <- function() {
  .d <- file.path(tempdir(), "nlmixr2estSens")
  if (!dir.exists(.d)) {
    dir.create(.d, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(.d)) return(tempdir())
  }
  .d
}
