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
#' Rows are assigned to endpoints by rxode2's own `etTrans()` (DVID/CMT ->
#' the endpoint's `predDf$cmt`), so each endpoint is judged on its own DV scale.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return a list of numeric DV vectors, one per `predDf` row (`NULL` when the
#'   rows of that endpoint cannot be identified)
#' @noRd
.addErrScaleEndpointDv <- function(ui, data) {
  .predDf <- ui$predDf
  .none <- vector("list", nrow(.predDf))
  .et <- tryCatch(
    suppressWarnings(suppressMessages(as.data.frame(rxode2::etTrans(
      data,
      rxode2::rxModelVars(ui$simulationModel),
      addCmt = TRUE,
      dropUnits = TRUE,
      allTimeVar = TRUE
    )))),
    error = function(e) NULL
  )
  if (is.null(.et) || !all(c("EVID", "DV", "CMT") %in% names(.et))) {
    return(.none)
  }
  .obs <- .et$EVID == 0 & !is.na(.et$DV)
  if (nrow(.predDf) == 1L) {
    return(list(.et$DV[.obs]))
  }
  lapply(seq_len(nrow(.predDf)), function(i) {
    .et$DV[.obs & .et$CMT == .predDf$cmt[i]]
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
  .dvs <- .addErrScaleEndpointDv(ui, data)
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
        "'",
        .predDf$cond[i],
        "' add sd '",
        .iniDf$name[.w],
        "'=",
        signif(.sd, 2),
        " < 1% median |DV| ",
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
