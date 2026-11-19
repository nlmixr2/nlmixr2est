#' Warn when a combined error model's additive SD starts negligibly small
#'
#' A tiny additive start relative to the data scale can trap a combined
#' (add + prop/pow) model in a basin where the additive term collapses (#1383).
#' The initial estimate is never changed.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the ui (unmodified; this hook only warns)
#' @noRd
NULL

#' Observation rows' DV for each endpoint of `predDf`
#'
#' @param predDf the ui's `predDf`
#' @inheritParams nlmixr2
#' @return a list of numeric DV vectors, one per `predDf` row (`NULL` when the
#'   rows of that endpoint cannot be identified)
#' @noRd
.addErrScaleEndpointDv <- function(predDf, data) {
  .nms <- tolower(names(data))
  .col <- function(name) {
    .w <- which(.nms == name)
    if (length(.w) == 1L) data[[.w]] else NULL
  }
  .dv <- suppressWarnings(as.numeric(.col("dv")))
  if (length(.dv) == 0L) {
    return(vector("list", nrow(predDf)))
  }
  .obs <- !is.na(.dv)
  .evid <- .col("evid")
  .mdv <- .col("mdv")
  .amt <- .col("amt")
  if (!is.null(.evid)) {
    .obs <- .obs & !is.na(.evid) & .evid == 0
  } else if (!is.null(.mdv)) {
    .obs <- .obs & !is.na(.mdv) & .mdv == 0
  } else if (!is.null(.amt)) {
    .obs <- .obs & (is.na(.amt) | .amt == 0)
  }
  if (nrow(predDf) == 1L) {
    return(list(.dv[.obs]))
  }
  .dvid <- .col("dvid")
  .cmt <- .col("cmt")
  lapply(seq_len(nrow(predDf)), function(i) {
    if (!is.null(.dvid)) {
      .w <- as.character(.dvid) %in% c(as.character(predDf$dvid[i]), predDf$cond[i])
    } else if (!is.null(.cmt)) {
      .w <- as.character(.cmt) %in% c(as.character(predDf$cmt[i]), predDf$cond[i])
    } else {
      return(NULL)
    }
    .dv[.obs & .w]
  })
}

.preProcessAddErrScaleWarn <- function(ui, est, data, control) {
  if (isTRUE(est %in% c(.preProcessHooksSimulationEst, "posthoc", "output")) || !is.data.frame(data)) {
    return(list(ui = ui))
  }
  .predDf <- tryCatch(ui$predDf, error = function(e) NULL)
  if (is.null(.predDf) || nrow(.predDf) == 0L) {
    return(list(ui = ui))
  }
  .comb <- which(
    .predDf$distribution == "norm" &
      .predDf$transform == "untransformed" &
      .predDf$errType %in% c("add + prop", "add + pow")
  )
  if (length(.comb) == 0L) {
    return(list(ui = ui))
  }
  .iniDf <- ui$iniDf
  .dvs <- .addErrScaleEndpointDv(.predDf, data)
  for (i in .comb) {
    .w <- which(.iniDf$condition == .predDf$cond[i] & .iniDf$err == "add" & !is.na(.iniDf$ntheta))
    if (length(.w) != 1L || isTRUE(.iniDf$fix[.w])) {
      next
    }
    .dv <- .dvs[[i]]
    if (length(.dv) == 0L || !any(.dv > 0)) {
      next
    }
    .med <- stats::median(abs(.dv))
    .sd <- .iniDf$est[.w]
    if (isTRUE(.predDf$variance[i])) {
      .sd <- sqrt(abs(.sd))
    }
    if (!is.finite(.med) || .med <= 0 || !is.finite(.sd) || .sd >= 0.01 * .med) {
      next
    }
    warning(
      paste0(
        "add sd '",
        .iniDf$name[.w],
        "'=",
        signif(.sd, 2),
        " tiny vs median |DV| ",
        signif(.med, 3),
        "; try ~",
        signif(0.1 * .med, 2)
      ),
      call. = FALSE
    )
  }
  list(ui = ui)
}
preProcessHooksAdd(".preProcessAddErrScaleWarn", .preProcessAddErrScaleWarn)
