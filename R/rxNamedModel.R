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
#' @param ui rxode2 ui whose data the model solves; when given, the model's
#'   string-literal numbering is pinned to the data's (#1299)
#' @return the compiled rxode2 model
#' @noRd
.nlmixr2estRxode2 <- function(model, role, ..., ui = NULL) {
  ## `role` no longer shapes the artifact, and the srcref strip that used to live here
  ## moved to rxode2 (rxStripModelSrc), applied while it builds the model.  Both need
  ## rxode2 >= 5.1.6 (DESCRIPTION).
  .strCmp <- if (is.null(ui)) list() else .nlmixr2estStrCmp(ui)
  model <- .nlmixr2estStrCmpSeed(model, .strCmp)
  .ret <- rxode2::rxode2(model, ...)
  .bad <- .nlmixr2estStrCmpBad(rxode2::rxModelVars(.ret), .strCmp)
  if (length(.strCmp) > 0L && length(.bad) > 0L) {
    stop("string literal order differs from the data for: ", paste(.bad, collapse = ", "), call. = FALSE)
  }
  ## Record which event-sensitivity mode built this model.  The focei model bundle is
  ## cached as rxNorm() TEXT and rehydrated with rxode2() (.foceiModelCacheInflate),
  ## and eventSens is NOT recoverable from a compiled model -- rxModelVars() does not
  ## carry it -- so without this tag a rehydrated sensitivity model is rebuilt in "fd"
  ## mode and every dosing-parameter (f/alag/rate/dur) sensitivity silently becomes
  ## zero.  Only set when the caller passed one: rxode2() resolves an absent eventSens
  ## through getOption("rxode2.eventSens"), and `missing(eventSens)` changes what it
  ## builds, so "not passed" has to replay as "not passed".
  .es <- list(...)$eventSens
  if (!is.null(.es)) {
    attr(.ret, "nlmixr2estEventSens") <- .es
  }
  .ret
}

#' String-literal table the estimation data is coded against
#'
#' `.foceiPreProcessData()` codes string-compared covariates with the normModel's
#' `strCmpParams` levels.  rxode2 numbers literals in first-seen order per compiled
#' model, so a generated model whose terms were reordered decodes them differently
#' (#1299).
#' @param ui rxode2 ui
#' @return named list of factors (rxode2 `strCmpParams`), possibly empty
#' @noRd
.nlmixr2estStrCmp <- function(ui) {
  .ret <- tryCatch(ui$mv0$strCmpParams, error = function(e) NULL)
  if (is.null(.ret)) list() else .ret
}

#' Covariates whose literal numbering in `mv` differs from `strCmp`
#' @param mv model variables of the generated model
#' @param strCmp the data's literal table (`.nlmixr2estStrCmp()`)
#' @return character vector of covariate names (model spelling) that differ
#' @noRd
.nlmixr2estStrCmpBad <- function(mv, strCmp) {
  .cur <- mv$strCmpParams
  if (length(.cur) == 0L) {
    return(character(0))
  }
  .canon <- tolower(names(strCmp))
  .bad <- vapply(
    names(.cur),
    function(.n) {
      .i <- match(tolower(.n), .canon)
      if (is.na(.i)) {
        return(TRUE)
      }
      .want <- levels(strCmp[[.i]])
      !identical(levels(.cur[[.n]])[seq_along(.want)], .want)
    },
    logical(1),
    USE.NAMES = FALSE
  )
  names(.cur)[.bad]
}

#' Pin a generated model's string-literal numbering to the data's
#'
#' Adds one suppressed `rx_strcmp_<i>_ ~ ...` statement per covariate, right after
#' the `params()` statement, that compares the covariate to every literal in the
#' data's order, so rxode2 numbers them that way.  Text whose numbering already
#' matches is returned unchanged.
#' @param model generated model text
#' @param strCmp the data's literal table (`.nlmixr2estStrCmp()`)
#' @return model text
#' @noRd
.nlmixr2estStrCmpSeed <- function(model, strCmp) {
  if (length(strCmp) == 0L || !is.character(model)) {
    return(model)
  }
  .txt <- paste(model, collapse = "\n")
  .bad <- .nlmixr2estStrCmpBad(rxode2::rxModelVars(.txt), strCmp)
  if (length(.bad) == 0L) {
    return(model)
  }
  .canon <- tolower(names(strCmp))
  .seed <- vapply(
    seq_along(.bad),
    function(.i) {
      .lv <- levels(strCmp[[match(tolower(.bad[.i]), .canon)]])
      paste0(
        "rx_strcmp_",
        .i,
        "_ ~ ",
        paste0("(", .bad[.i], " == ", encodeString(.lv, quote = "\""), ")", collapse = " + ")
      )
    },
    character(1)
  )
  .lines <- strsplit(.txt, "\n", fixed = TRUE)[[1]]
  .w <- which(grepl("^[ \t]*params?[ \t]*\\(", .lines))
  .at <- 0L
  if (length(.w) > 0L) {
    # after the line that closes params(), which may span lines
    .n <- function(ch) lengths(regmatches(.lines, gregexpr(ch, .lines, fixed = TRUE)))
    .depth <- cumsum(.n("(") - .n(")"))
    .before <- c(0L, .depth)[.w[1]]
    .close <- which(seq_along(.lines) >= .w[1] & .depth <= .before)
    .at <- if (length(.close) > 0L) .close[1] else .w[1]
  }
  paste(append(.lines, .seed, after = .at), collapse = "\n")
}
