#' Preprocess Covariates needed (or other data items)
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the ui (possibly modified)
#' @export
#' @author Matthew L. Fidler
.nlmixr0preProcessCovariatesPresent <- function(ui, est, data, control) {
  # Could later support stacking data or custom DV/IDV columns.

  # optional=TRUE: an unregistered method name (e.g. foceiLikLoad(est=) naming
  # a consumer whose package is not loaded) has no attributes to consult, which
  # means "no covPresent requirement", not an error (#939) -- the same
  # treatment .nlmixr2PriorSupport() gives an unknown method
  if (!missing(data) &&
        length(data) > 0L &&
        isTRUE(attr(utils::getS3method("nlmixr2Est", est, optional = TRUE),
                    "covPresent"))) {
    .covNames <- ui$covariates
    .newNames <- .nmUpcaseNonCov(names(data), .covNames)
    colnames(data) <- .newNames
    requiredCols <- c("TIME", .covNames)
    checkmate::assert_names(.newNames, must.include = requiredCols)
  }
  NULL
}

preProcessHooksAdd(".nlmixr0preProcessCovariatesPresent", .nlmixr0preProcessCovariatesPresent)
