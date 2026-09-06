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
#' @export
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
  .par <- .saemThetaSensParMap(.ui, .mod)
  if (is.null(.par)) return(NULL)
  c(list(thetaSens = .mod), .map, .par)
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

#' THETA[k]/ETA[k] -> SAEM phi-column map for the sensitivity peer
#'
#' The sensitivity model is FOCEi codegen, so its parameters are `THETA[]` and
#' `ETA[]` while SAEM carries a phi matrix.  Driving it from SAEM's own phi
#' therefore needs the same translation `.saemPhi1TargetMap()` builds for the
#' phi1 peers -- which only exists for a general-likelihood fit, and a plain
#' normal model wants the gradient refinement just as much.
#'
#' `kind` is 1 for a mu-referenced (phi1) parameter, 0 for a non-mu (phi0) one
#' and -1 for a theta that is fixed and so never appears in phi at all (its
#' `ini()` value is carried in `fixedVal`).
#'
#' @param ui rxode2 ui
#' @param mod the compiled sensitivity model
#' @return list of `thetaKind`, `thetaCol`, `thetaFixedVal`, `etaCol`, `dvCol`,
#'   or `NULL` when the model's parameter shape is not one this can translate
#' @noRd
#' @author Matthew L. Fidler
.saemThetaSensParMap <- function(ui, mod) {
  .pars <- rxode2::rxParam(mod)
  ## anything beyond THETA[]/ETA[]/DV means the phi translation below would be
  ## guessing at what to put in the extra slots
  .other <- .pars[!(grepl("^(THETA|ETA)\\[", .pars) | .pars == "DV")]
  if (length(.other) > 0) return(NULL)
  .dvCol <- match("DV", .pars) - 1L
  if (is.na(.dvCol)) return(NULL)
  .iniDf <- ui$iniDf
  .parsAll <- tryCatch(rxUiGet.saemParamsToEstimateCov(list(ui)),
                       error = function(e) NULL)
  if (is.null(.parsAll)) return(NULL)
  .muRef <- ui$muRefDataFrame
  .isPhi1 <- .parsAll %in% .muRef$theta
  .phi1Names <- .parsAll[.isPhi1]
  .phi0Names <- .parsAll[!.isPhi1]
  .nTheta <- length(grep("^THETA\\[", .pars))
  .thetaKind <- integer(.nTheta)
  .thetaCol <- integer(.nTheta)
  .thetaFixedVal <- numeric(.nTheta)
  for (.k in seq_len(.nTheta)) {
    .nm <- .iniDf$name[!is.na(.iniDf$ntheta) & .iniDf$ntheta == .k]
    if (length(.nm) != 1L) return(NULL)
    if (.nm %in% .phi1Names) {
      .thetaKind[.k] <- 1L
      .thetaCol[.k] <- match(.nm, .phi1Names) - 1L
    } else if (.nm %in% .phi0Names) {
      .thetaKind[.k] <- 0L
      .thetaCol[.k] <- match(.nm, .phi0Names) - 1L
    } else {
      .thetaKind[.k] <- -1L
      .est <- .iniDf$est[.iniDf$name == .nm]
      if (length(.est) != 1L || is.na(.est)) return(NULL)
      .thetaFixedVal[.k] <- .est
    }
  }
  .nEta <- length(grep("^ETA\\[", .pars))
  .etaDiag <- !is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2
  .etaNames <- .iniDf$name[.etaDiag][order(.iniDf$neta1[.etaDiag])]
  .etaCol <- integer(.nEta)
  for (.k in seq_len(.nEta)) {
    if (.k > length(.etaNames)) return(NULL)
    ## SAEM's phi1 columns are the mu-referenced parameters, in that order; an
    ## eta maps to its own parameter's column
    .th <- .muRef$theta[match(.etaNames[.k], .muRef$eta)]
    if (is.na(.th)) return(NULL)
    .m <- match(.th, .phi1Names)
    if (is.na(.m)) return(NULL)
    .etaCol[.k] <- as.integer(.m) - 1L
  }
  list(thetaKind = as.integer(.thetaKind), thetaCol = as.integer(.thetaCol),
       thetaFixedVal = as.numeric(.thetaFixedVal),
       etaCol = as.integer(.etaCol), dvCol = as.integer(.dvCol))
}
