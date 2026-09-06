#' Theta-sensitivity model for SAEM's non-mu (phi0) theta refinement
#'
#' `refinePhi0Lik()` moves the non-mu thetas by a derivative-free search that
#' spends a full population solve on every objective evaluation and has no
#' derivative information at all.  nlmixr2 already emits exact symbolic
#' sensitivities, so one solve of this model yields `d(f)/d(theta)` for every
#' such theta at once -- see `src/nonMuThetaGrad.h` for why that is both cheaper
#' than the search and better directed, and for the schedule the two share.
#'
#' The model itself is FOCEi codegen (`THETA[]`/`ETA[]` parameters), exactly as
#' `saemPhi1Inner`'s `innerHess2`/`predNoLhs` peers are, so it is driven from
#' SAEM's own phi matrix through the same kind of column map.
#'
#' @param x rxode2 ui, in a list
#' @return list with the compiled `thetaSens` model and its phi-column map, or
#'   `NULL` when the model shape is out of scope (the refinement then falls back
#'   to the search alone, which is what it did before this existed)
#' @noRd
#' @author Matthew L. Fidler
rxUiGet.saemThetaSens <- function(x, ...) {
  .ui <- x[[1]]
  if (!isTRUE(tryCatch(as.logical(rxode2::rxGetControl(.ui, "nonMuThetaGrad", TRUE)),
                       error = function(e) FALSE))) {
    return(NULL)
  }
  ## the same shape restrictions saemPhi1Inner's map imposes: a covariate
  ## mu-group or a non-mu ETA means the phi columns are not a plain
  ## one-per-parameter map and the column translation below would be wrong
  if (length(.ui$nonMuEtas) > 0) return(NULL)
  .cov <- tryCatch(rxUiGet.saemMuRefCovariateDataFrame(list(.ui)),
                   error = function(e) NULL)
  if (is.null(.cov) || length(.cov$covariateParameter) > 0) return(NULL)
  .mod <- tryCatch(.impmapThetaSensModel(.ui), error = function(e) NULL)
  if (is.null(.mod)) return(NULL)
  .map <- .saemThetaSensMap(.ui)
  if (is.null(.map)) return(NULL)
  c(list(thetaSens = .mod), .map)
}
attr(rxUiGet.saemThetaSens, "rstudio") <- emptyenv()

#' Map each sensitivity output to the SAEM phi0 column it differentiates
#'
#' `.impmapThetaSensModel()` emits one `rx__sens_rx_pred__BY_THETA_j___` per
#' estimated non-mu theta, in `.impmapEstTheta(ui)$all` order.  SAEM optimizes
#' phi0 COLUMN values, so the refinement needs each output's phi0 column --
#' `-1` for an output that is not a phi0 column at all (it is then skipped
#' rather than misapplied).
#'
#' @param ui rxode2 ui
#' @return list with `sensPhi0Col` (0-based, `-1` = not phi0) and `sensTheta`
#'   (1-based `ntheta` of each output), or `NULL`
#' @noRd
#' @author Matthew L. Fidler
.saemThetaSensMap <- function(ui) {
  .est <- tryCatch(.impmapEstTheta(ui)$all, error = function(e) NULL)
  if (is.null(.est) || length(.est) == 0L) return(NULL)
  .iniDf <- ui$iniDf
  .parsAll <- tryCatch(rxUiGet.saemParamsToEstimateCov(list(ui)),
                       error = function(e) NULL)
  if (is.null(.parsAll)) return(NULL)
  .muRef <- ui$muRefDataFrame
  .phi0Names <- .parsAll[!(.parsAll %in% .muRef$theta)]
  .col <- integer(length(.est))
  for (.i in seq_along(.est)) {
    .nm <- .iniDf$name[!is.na(.iniDf$ntheta) & .iniDf$ntheta == .est[.i]]
    if (length(.nm) != 1L) return(NULL)
    .m <- match(.nm, .phi0Names)
    .col[.i] <- if (is.na(.m)) -1L else (as.integer(.m) - 1L)
  }
  ## nothing to refine by gradient if no output lands on a phi0 column
  if (all(.col < 0L)) return(NULL)
  list(sensPhi0Col = as.integer(.col), sensTheta = as.integer(.est), ok = TRUE)
}
