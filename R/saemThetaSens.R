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
  if (!isTRUE(tryCatch(as.logical(rxode2::rxGetControl(.ui, "nonMuThetaGrad", FALSE)),
                       error = function(e) FALSE))) {
    return(NULL)
  }
  ## A covariate mu-group splits one parameter across several phi columns, so
  ## the one-per-parameter translation below would be wrong.  A non-mu ETA is
  ## NOT excluded here: `.saemPhi1Split()` resolves those through the same map
  ## SAEM itself uses, and they are exactly the case this refinement exists for
  ## -- a `dist()`-declared eta has no additive `theta + eta` form and so is
  ## always classified `nonMuEta`.
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
  .split <- .saemPhi1Split(ui)
  if (is.null(.split)) return(NULL)
  .phi0Names <- .split$parsAll[!.split$isPhi1]
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
  ## DV is a model INPUT only for a general-likelihood model, where the
  ## likelihood expression names it.  A normal-error model (`prop()`, `add()`,
  ## ...) predicts rx_pred_/rx_r_ and the observation comes from SAEM's own
  ## data, so no DV parameter exists -- -1 tells the accumulate loop to take y
  ## from there rather than writing a parameter slot that is not present.
  .dvCol <- match("DV", .pars) - 1L
  if (is.na(.dvCol)) .dvCol <- -1L
  .iniDf <- ui$iniDf
  .split <- .saemPhi1Split(ui)
  if (is.null(.split)) return(NULL)
  .parsAll <- .split$parsAll
  .phi1Names <- .parsAll[.split$isPhi1]
  .phi0Names <- .parsAll[!.split$isPhi1]
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
  if (.nEta > length(.split$etaPhi1Col)) return(NULL)
  ## ETA[k] is the k'th eta in neta order, which is the order .saemPhi1Split()
  ## returns its phi1 columns in
  .etaCol <- .split$etaPhi1Col[seq_len(.nEta)]
  ## Which ETA[k] carries its parameter's whole value rather than a deviation
  ## from a THETA.
  ##
  ## The pooled setter reproduces a mu-referenced parameter by putting the
  ## combined phi value in THETA[k] and 0 in ETA[k].  That is only right when
  ## the parameter HAS a theta.  A dist()-declared eta is a nonMuEta: the
  ## parameter IS the eta, no THETA[] slot refers to it, and zeroing ETA[k]
  ## silently evaluates the model at a latent eta of 0 for every subject.
  ## Flag those so the setter puts the phi value in ETA[k] instead.
  .etaNonMu <- vapply(seq_along(.etaCol), function(.k) {
    !any(.thetaKind == 1L & .thetaCol == .etaCol[.k])
  }, logical(1))
  list(thetaKind = as.integer(.thetaKind), thetaCol = as.integer(.thetaCol),
       thetaFixedVal = as.numeric(.thetaFixedVal),
       etaCol = as.integer(.etaCol), dvCol = as.integer(.dvCol),
       etaNonMu = as.integer(.etaNonMu))
}
#' Split SAEM's parameter vector into its phi1 and phi0 groups
#'
#' SAEM itself makes this split from the omega diagonal (`covstruct`):
#' `i1 <- grep(1, diag(covstruct))` in `.configsaem()`.  A parameter is phi1
#' exactly when it carries a random effect.
#'
#' The obvious proxy -- `parsAll %in% muRefDataFrame$theta` -- gets that right
#' only for mu-referenced etas.  A `dist()`-declared eta has no additive
#' `theta + eta` form, so rxode2's mu-ref scanner classifies it as a
#' `nonMuEta`, it has no `muRefDataFrame` row, and the proxy files it under
#' phi0 even though SAEM samples it in phi1.  `saemEtaTrans` is the map SAEM
#' actually uses (`.saemEtaTrans()`), and it resolves both kinds, so derive
#' membership from it instead.
#'
#' @param ui rxode2 ui
#' @return list with `parsAll`, `isPhi1` (logical over `parsAll`), `etaPhi1Col`
#'   (0-based phi1 column per eta, in neta order) and `etaNames`, or `NULL`
#'   when the map does not resolve
#' @noRd
#' @author Matthew L. Fidler
.saemPhi1Split <- function(ui) {
  .parsAll <- tryCatch(rxUiGet.saemParamsToEstimateCov(list(ui)),
                       error = function(e) NULL)
  if (is.null(.parsAll)) return(NULL)
  .trans <- tryCatch(rxUiGet.saemEtaTrans(list(ui)), error = function(e) NULL)
  if (is.null(.trans) || anyNA(.trans)) return(NULL)
  ## a negative index is .saemEtaTrans()'s nonMu=TRUE (pred-model) encoding;
  ## rxUiGet.saemEtaTrans is the nonMu=FALSE form and must not produce one
  if (any(.trans < 1L) || any(.trans > length(.parsAll))) return(NULL)
  .isPhi1 <- seq_along(.parsAll) %in% .trans
  .iniDf <- ui$iniDf
  .etaDiag <- !is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2
  .etaNames <- .iniDf$name[.etaDiag][order(.iniDf$neta1[.etaDiag])]
  if (length(.etaNames) != length(.trans)) return(NULL)
  .etaPhi1Col <- match(.trans, which(.isPhi1)) - 1L
  if (anyNA(.etaPhi1Col)) return(NULL)
  list(parsAll = .parsAll, isPhi1 = .isPhi1,
       etaPhi1Col = as.integer(.etaPhi1Col), etaNames = .etaNames)
}

#' Prediction-only peer for a NORMAL model, so its solves share the peer layout
#'
#' The exact-gradient refinement needs the theta-sensitivity model to be
#' solvable alongside SAEM's own.  It is not, on the ordinary path: `rxSolve_`
#' lays out ONE parameter vector -- the solved model's -- and the two models do
#' not share one.  The sensitivity peer declares `THETA[k]`/`ETA[k]`; SAEM's own
#' model declares native names (`lclm`, `rxz.eta.cl`, ...).  Sizing the pool for
#' the peer and switching to SAEM's model makes it read the peer's slots as its
#' own parameters (measured: an immediate segfault).
#'
#' A general-likelihood fit does not have this problem, because it already
#' routes SAEM's own likelihood read through `predNoLhs` -- a FOCEi-codegen peer
#' that shares the `THETA[]`/`ETA[]` declaration.  `predNoLhs` is ordinary FOCEi
#' codegen and exists for a normal model too; only `.saemPhi1TargetMap()`'s
#' hard requirement of a `DV` parameter kept it out of reach, and a normal model
#' has none (it predicts, and the observation comes from the data).
#'
#' So build the same peer here.  With it, every model in the pool speaks
#' `THETA[]`/`ETA[]`, the pool can be sized by whichever measures widest, and
#' the sensitivity peer becomes solvable -- which is what lets the refinement
#' engage at all.
#'
#' `dvCol = -1` is deliberate and load-bearing: it is what keeps
#' `_saemPhi1PoolReady` FALSE in the C++, so this turns on the pooled SOLVE
#' routing without also turning on the phi1 theta refinement, which is a
#' general-likelihood step and is not wanted here.
#'
#' @param x rxode2 ui, in a list
#' @return list with `predNoLhs` and its phi-column map, or `NULL` when out of
#'   scope (SAEM then takes its original single-model path unchanged)
#' @noRd
#' @author Matthew L. Fidler
#' @export
rxUiGet.saemOwnPred <- function(x, ...) {
  .ui <- x[[1]]
  ## only for a NORMAL model -- a general-likelihood fit already has this peer
  ## through saemPhi1Inner, and building a second one would fight it
  if (.saemGeneralLik(.ui)) return(NULL)
  ## no point paying for the peer unless the sensitivity model it exists to
  ## make solvable actually resolved
  if (!isTRUE(tryCatch(as.logical(rxode2::rxGetControl(.ui, "nonMuThetaGrad", FALSE)),
                       error = function(e) FALSE))) {
    return(NULL)
  }
  .fm <- tryCatch(.ui$focei, error = function(e) NULL)
  if (is.null(.fm)) return(NULL)
  .pred <- .fm$predNoLhs
  if (is.null(.pred)) return(NULL)
  .par <- .saemThetaSensParMap(.ui, .pred)
  if (is.null(.par)) return(NULL)
  ## the pooled read drives every phi1 column through ETA[]; a mismatch here
  ## would silently mis-map columns rather than fail
  .nphi1 <- sum(.saemPhi1Split(.ui)$isPhi1)
  if (length(.par$etaCol) != .nphi1) return(NULL)
  c(list(predNoLhs = .pred, ok = TRUE), .par)
}
attr(rxUiGet.saemOwnPred, "rstudio") <- emptyenv()
