#' Preprocess Covariates needed (or other data items)
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the ui (possibly modified)
#' @export
#' @author Matthew L. Fidler
.nlmixr0preProcessCovariatesPresent <- function(ui, est, data, control) {
  # Could possibly use to stack data or use an DV or IDV different
  # than what is present in the data

  if (!missing(data) &&
        length(data) > 0L &&
         isTRUE(attr(utils::getS3method("nlmixr2Est", est), "covPresent"))) {
    .covNames <- ui$covariates
    .newNames <- vapply(names(data), function(x) {
      if (any(x == .covNames)) {
        x
      } else {
        toupper(x)
      }
    }, character(1))
    colnames(data) <- .newNames
    requiredCols <- c("TIME", .covNames)
    checkmate::assert_names(.newNames, must.include = requiredCols)
  }
  NULL
}

preProcessHooksAdd(".nlmixr0preProcessCovariatesPresent", .nlmixr0preProcessCovariatesPresent)
