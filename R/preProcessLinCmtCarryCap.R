#' Carry columns a `linCmt()` model's eta-gradient carry would get wrong
#'
#' Zero when the model cannot be past `cap` (no `linCmt()`, no covariate, or
#' too few etas to fill the columns), so only a model that may be past the
#' cap pays for the symbolic eligibility check.  Past the cap, the build
#' either stops or drops pairs (`.rxFoceiLinCmtCarryCapKeep()`); a dropped
#' pair whose covariates are constant within every subject keeps an exact
#' standard gradient, so only a dropped pair the data shows (or may show)
#' varying makes the model need the ODE translation.
#'
#' @param ui rxode2 ui
#' @param cap the carry columns the loaded rxode2 has
#' @param data the fit's data, or NULL
#' @return integer, the columns the carry needs when the build would stop or
#'   approximate a varying pair, else 0
#' @noRd
#' @author Matthew L. Fidler
.rxFoceiLinCmtCarryColumnsNeeded <- function(ui, cap, data = NULL) {
  if (rxode2::.rxLinNcmt(ui)["numLin"] <= 0L || length(ui$allCovs) == 0L) {
    return(0L)
  }
  .iniDf <- ui$iniDf
  .nEta <- sum(!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2)
  # one column per pair, one tracker per slot (7 slots), one amounts and one
  # lag column at most
  if (.nEta + min(.nEta, 7L) + 2L <= cap) {
    return(0L)
  }
  .pairs <- tryCatch(.foceiLinCmtCarryPairs(ui), error = function(e) NULL)
  if (isFALSE(attr(.pairs, "modelOk"))) {
    return(0L)
  }
  .pairs <- .rxFoceiLinCmtCarryJumpFilter(.pairs)
  if (is.null(.pairs)) {
    return(0L)
  }
  .need <- as.integer(.rxFoceiLinCmtCarryColumns(.pairs))
  .keep <- .rxFoceiLinCmtCarryCapKeep(.pairs, cap)
  if (is.null(.keep)) {
    return(.need)
  }
  .vary <- vapply(
    strsplit(.pairs$covs[!.keep], ",", fixed = TRUE),
    .linCmtCarryMayVary,
    logical(1),
    data = data
  )
  if (any(.vary)) .need else 0L
}

#' Could any of these covariates vary within a subject in the data?
#'
#' A covariate the data cannot settle (no data, or no matching column) may.
#' @param covs covariate names
#' @param data the fit's data, or NULL
#' @return logical
#' @noRd
#' @author Matthew L. Fidler
.linCmtCarryMayVary <- function(covs, data) {
  if (!is.data.frame(data)) {
    return(TRUE)
  }
  .n <- names(data)
  any(vapply(
    covs,
    function(.cv) {
      .col <- .n[.n == .cv]
      if (length(.col) != 1L) {
        .col <- .n[tolower(.n) == tolower(.cv)]
      }
      if (length(.col) != 1L) {
        return(TRUE)
      }
      !isFALSE(.rxFoceiCarryCovVaries(data, .col))
    },
    logical(1)
  ))
}

#' The compartment each number names, in a ui's numbering
#'
#' States by name, then the endpoints that are not states by their position
#' in `ui$predDf` (`linToOde()` keeps their order but renames a bare
#' `linCmt()` endpoint).
#' @param ui rxode2 ui
#' @return character, indexed by compartment number (NA where unnamed)
#' @noRd
#' @author Matthew L. Fidler
.linCmtCarryCmtNames <- function(ui) {
  .state <- rxode2::rxModelVars(ui)$state
  .pd <- ui$predDf
  .ret <- rep(NA_character_, max(c(length(.state), .pd$cmt)))
  .ret[.pd$cmt] <- paste0("(endpoint ", seq_len(nrow(.pd)), ")")
  .ret[seq_along(.state)] <- .state
  .ret
}

#' Renumber a numeric `cmt` data column from one ui's numbering to another's
#'
#' `linToOde()` puts the peripheral compartments right after depot/central,
#' ahead of the endpoints, so a numeric `cmt` naming an endpoint would point
#' at a peripheral compartment.
#' @param data the fit's data
#' @param from the ui the data was written against
#' @param to the translated ui
#' @return `data` with its numeric `cmt` column renumbered
#' @noRd
#' @author Matthew L. Fidler
.linCmtCarryRemapCmt <- function(data, from, to) {
  .w <- which(tolower(names(data)) == "cmt")
  if (length(.w) != 1L || !is.numeric(data[[.w]])) {
    return(data)
  }
  .from <- .linCmtCarryCmtNames(from)
  .to <- .linCmtCarryCmtNames(to)
  .c <- data[[.w]]
  .k <- abs(.c)
  .new <- match(.from[.k], .to)
  .ok <- !is.na(.c) & .k >= 1L & .k <= length(.from) & !is.na(.new)
  .c[.ok] <- sign(.c[.ok]) * .new[.ok]
  data[[.w]] <- .c
  data
}

#' Solve a `linCmt()` model past the carry cap as ODEs
#'
#' The eta-gradient carry has `.rxFoceiLinCmtCarryMaxPairs()` columns.  The
#' etas past them would keep the standard gradient, which is wrong when a
#' covariate changes a `linCmt()` parameter within a subject and makes the
#' objective depend on the order the etas are declared in (#1318).  The ODE
#' translation has exact sensitivities for every eta.  Named to run after
#' `.uiApplyIov()` and `.uiApplyMu2()`, so it sees the model the build does.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the translated ui and renumbered data, or `NULL` when
#'   nothing to do
#' @export
#' @author Matthew L. Fidler
.uiApplyOdeLinCmtCarryCap <- function(ui, est, data, control) {
  if (!(est %in% .linCmtFoceiEstFamily)) {
    return(NULL)
  }
  # .preProcessLinCmtOde() translates these
  if (.uiIsMixedLinCmtOde(ui)) {
    return(NULL)
  }
  .ctl <- if (is.list(control)) control else list()
  .carry <- .ctl$linCmtSensCarry
  if (is.null(.carry)) {
    .carry <- "auto"
  }
  if (!.rxFoceiLinCmtCarryOn(.carry, .ctl$rxControl$covsInterpolation)) {
    return(NULL)
  }
  .cap <- .rxFoceiLinCmtCarryMaxPairs()
  .need <- .rxFoceiLinCmtCarryColumnsNeeded(ui, .cap, data)
  if (.need <= .cap) {
    return(NULL)
  }
  .ui <- try(rxode2::linToOde(ui), silent = TRUE)
  if (inherits(.ui, "try-error")) {
    return(NULL)
  }
  warning(
    sprintf("linCmt() solved as ODEs: etas need %d of %d carry columns", .need, .cap),
    call. = FALSE
  )
  if (is.data.frame(data)) {
    list(ui = .ui, data = .linCmtCarryRemapCmt(data, ui, .ui))
  } else {
    list(ui = .ui)
  }
}

preProcessHooksAdd(".uiApplyOdeLinCmtCarryCap", .uiApplyOdeLinCmtCarryCap)
