#' Column of a dataset by case-insensitive name
#'
#' @param data data frame
#' @param name lower-case column name
#' @return the column, or `NULL` when absent
#' @noRd
.nlmixrDataCol <- function(data, name) {
  .w <- which(tolower(names(data)) == name)
  if (length(.w) == 0L) {
    NULL
  } else {
    data[[.w[1]]]
  }
}

#' Observation records, as etTrans() reads them: EVID (less MDV == 1), else
#' MDV, else AMT
#'
#' @param data data frame
#' @return logical vector, one per row
#' @noRd
.nlmixrDataObsRows <- function(data) {
  .evid <- .nlmixrDataCol(data, "evid")
  .mdv <- .nlmixrDataCol(data, "mdv")
  .amt <- .nlmixrDataCol(data, "amt")
  if (!is.null(.evid)) {
    # a missing EVID with no dose is an observation
    .noDose <- if (is.null(.amt)) TRUE else is.na(.amt) | .amt %in% 0
    .ret <- .evid %in% 0 | (is.na(.evid) & .noDose)
    if (!is.null(.mdv)) {
      .ret <- .ret & !(.mdv %in% 1)
    }
    return(.ret)
  }
  if (!is.null(.mdv)) {
    return(.mdv %in% 0)
  }
  if (!is.null(.amt)) {
    return(is.na(.amt) | .amt == 0)
  }
  rep(TRUE, nrow(data))
}

#' Observation records of one endpoint
#'
#' A single endpoint scores every observation, whatever its cmt.  Otherwise a
#' record belongs to the endpoint named by its `DVID` (rank-converted, as
#' etTrans() does, when it exceeds the endpoint count), or by `CMT` (number or
#' name) when the record has no usable `DVID`.
#'
#' @param data data frame
#' @param predDf the ui's `predDf`
#' @param i row of `predDf`
#' @param obs observation rows (`.nlmixrDataObsRows()`)
#' @return logical vector, one per row
#' @noRd
.nlmixrDataEndpointRows <- function(data, predDf, i, obs = .nlmixrDataObsRows(data)) {
  if (nrow(predDf) <= 1L) {
    return(obs)
  }
  .pred1 <- predDf[i, ]
  .n <- nrow(data)
  .m <- rep(FALSE, .n)
  .useDvid <- rep(FALSE, .n)
  .dvid <- .nlmixrDataCol(data, "dvid")
  if (!is.null(.dvid)) {
    if (is.character(.dvid) || is.factor(.dvid)) {
      .dvid <- as.character(.dvid)
      .useDvid <- !is.na(.dvid)
      .m <- .useDvid & .dvid == .pred1$cond
    } else {
      # etTrans() skips the conversion when any DVID is missing
      if (length(.dvid) > 0L && !anyNA(.dvid) && max(.dvid) > nrow(predDf)) {
        .dvid <- match(.dvid, sort(unique(.dvid)))
      }
      # an out-of-range DVID falls back to CMT
      .useDvid <- !is.na(.dvid) & .dvid > 0 & .dvid <= nrow(predDf)
      .m <- .useDvid & .dvid == .pred1$dvid
    }
  }
  .cmt <- .nlmixrDataCol(data, "cmt")
  if (!is.null(.cmt)) {
    .mc <- if (is.character(.cmt) || is.factor(.cmt)) {
      !is.na(.cmt) & as.character(.cmt) == .pred1$cond
    } else {
      !is.na(.cmt) & .cmt == .pred1$cmt
    }
    .m <- .m | (!.useDvid & .mc)
  }
  obs & .m
}

#' Does this endpoint's transform need a positive DV?
#'
#' `lnorm()` always does; `boxCox()` unless its lambda is fixed at 1.
#'
#' @param ui rxode2 ui
#' @param pred1 one row of the ui's `predDf`
#' @return logical
#' @noRd
.nlmixrEndpointNeedsPositiveDv <- function(ui, pred1) {
  .tr <- as.character(pred1$transform)
  if (identical(.tr, "lnorm")) {
    return(TRUE)
  }
  if (!identical(.tr, "boxCox")) {
    return(FALSE)
  }
  .lam <- suppressWarnings(as.numeric(as.character(pred1$lambda)))
  if (isTRUE(.lam == 1)) {
    return(FALSE)
  }
  .iniDf <- ui$iniDf
  .w <- which(.iniDf$err %in% "boxCox" & .iniDf$condition %in% pred1$cond)
  !(length(.w) == 1L && isTRUE(.iniDf$fix[.w]) && isTRUE(.iniDf$est[.w] == 1))
}

#' Warn when a log/Box-Cox endpoint has non-positive observed DV
#'
#' `lnorm()` and `boxCox()` clamp a non-positive DV to a tiny positive value,
#' silently distorting the objective.  Censored rows (`CENS != 0`) are
#' skipped since their DV holds the limit.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the ui (unmodified; this hook only warns)
#' @noRd
.preProcessDvDomainWarn <- function(ui, est, data, control) {
  if (
    !is.data.frame(data) ||
      isTRUE(est %in% c("rxSolve", "simulate", "simulation", "predict", "output"))
  ) {
    return(list(ui = ui))
  }
  .predDf <- tryCatch(ui$predDf, error = function(e) NULL)
  if (is.null(.predDf) || nrow(.predDf) == 0L) {
    return(list(ui = ui))
  }
  .dv <- .nlmixrDataCol(data, "dv")
  if (is.null(.dv) || !is.numeric(.dv)) {
    return(list(ui = ui))
  }
  .multi <- nrow(.predDf) > 1L
  if (.multi && is.null(.nlmixrDataCol(data, "cmt")) && is.null(.nlmixrDataCol(data, "dvid"))) {
    return(list(ui = ui))
  }
  .obs <- .nlmixrDataObsRows(data) & !is.na(.dv)
  .cens <- .nlmixrDataCol(data, "cens")
  if (!is.null(.cens)) {
    .obs <- .obs & (is.na(.cens) | .cens == 0)
  }
  for (.i in seq_len(nrow(.predDf))) {
    .pred1 <- .predDf[.i, ]
    if (!.nlmixrEndpointNeedsPositiveDv(ui, .pred1)) {
      next
    }
    .n <- sum(.nlmixrDataEndpointRows(data, .predDf, .i, .obs) & .dv <= 0)
    if (.n > 0L) {
      warning(
        .n,
        " DV <= 0 on ",
        as.character(.pred1$transform),
        " endpoint '",
        .pred1$cond,
        "'",
        call. = FALSE
      )
    }
  }
  list(ui = ui)
}
preProcessHooksAdd(".preProcessDvDomainWarn", .preProcessDvDomainWarn)
