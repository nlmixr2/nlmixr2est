#' FOCEi inner model for SAEM's general-likelihood phi1 theta step
#'
#' Reuses FOCEi's own inner (eta-sensitivity-only) model, built by
#' `ui$focei` -- SAEM and FOCEi already read the same rxode2 ui object,
#' just via different `rxUiGet.*` accessors, so no separate model-generation
#' path is needed. The inner model's `rx_pred_` IS the log-density (mu and
#' eta are separate live model quantities -- eta is pushed in as a solved
#' parameter, not collapsed into phi the way `saemModel` flattens it), with
#' `rx__sens_rx_pred__BY_ETA_k___` (1st order, always available) and
#' `rx__d2pred_i_j__` (2nd order, when `$innerHess2` builds) giving the
#' eta-gradient and eta-Hessian a per-theta EBE-mode search and its Laplace
#' correction need.
#'
#' Deliberately reuses the INNER model, not the augmented OUTER sensitivity
#' model (`rxUiGet.foceiOuter`): the outer model's theta-direction
#' sensitivity chains exist solely to give FOCEi's own outer optimizer an
#' analytic gradient with respect to theta. SAEM's phi1 theta step is
#' optimized with bobyqa (derivative-free), so no theta-gradient is ever
#' needed -- only, for each theta bobyqa proposes, the inner (eta-only)
#' value/gradient/Hessian used to find the conditional mode and score it.
#'
#' `ui$focei$innerHess2` (the exact 2nd-order eta-Hessian model) is only built
#' when the ui's `saemControl(phi1Hessian=TRUE)` -- the default is `FALSE`
#' (an ablation check found the Laplace `log|H|` correction it feeds is not
#' load-bearing for convergence and can dominate/diverge for a heavy-tailed
#' `t()`/`cauchy()` endpoint, nlmixr2/nlmixr2est#999), so building it is
#' skipped by default to avoid the extra `foceiControl(fast=TRUE)` compile.
#' When requested, it requires `foceiControl(fast=TRUE)` and
#' `.foceiLLGradInScope(ui)`; this accessor temporarily installs
#' `foceiControl(fast=TRUE)` on the ui if it is not already FOCEi-shaped
#' (SAEM's own control has no `fast=` field), builds `ui$focei` under that
#' control, then restores the ui's original control -- `ui$focei` is a
#' cached accessor, so the built model list persists.
#'
#' `innerHess2` legitimately comes back `NULL` for some model shapes (a
#' `linCmt()` prediction is the known case -- FOCEi's own analytic-Hessian
#' build declines there too, per `.foceiMaybeAddHdEta2`'s own
#' `tryCatch(..., error = function(e) .s)` fallback). A caller must treat a
#' `NULL` `innerHess2` as "fall back to a finite-difference eta-Hessian for
#' this fit," not an error -- exactly FOCEi's own behavior; that fallback
#' uses `predNoLhs` (below), a bare no-sensitivity prediction-only model, the
#' same one `calcEtaHessian`'s own Shi(2021) finite-difference branch uses
#' (`odeSlotPred`).
#'
#' The correct log-density model is `$inner` for a literal `ll()` endpoint
#' (its `rx_pred_` is unconditionally the user's own log-density expression),
#' but `$innerLlik` for a `t()`/`cauchy()`/`dnorm()` *sugar* endpoint --
#' `$inner` there is FOCEi's closed-form Gaussian mean/variance
#' approximation (what FOCEi's own Newton search actually fits these
#' endpoints with; the true log-density form exists only for CWRES). This
#' accessor selects the right one automatically.
#'
#' @param x rxUiGet-style single-element list holding the ui
#' @param ... unused
#' @return `NULL` if not a general-likelihood endpoint. Otherwise a list
#'   with `inner` (the eta-sensitivity log-density model, always non-NULL),
#'   `innerHess2` (the exact eta-Hessian model, or `NULL` when out of scope
#'   for this model shape), `predNoLhs` (the bare, no-sensitivity
#'   prediction-only FD-fallback model, always non-NULL), and `ok` plus the
#'   `THETA[k]`/`ETA[k]` -> SAEM phi-column maps `thetaKind`/`thetaCol`/
#'   `thetaFixedVal`/`etaCol` (see `.saemPhi1TargetMap`) that Phase 4's C++
#'   theta step needs to drive `innerHess2`/`predNoLhs` from SAEM's own
#'   phi/theta state without any R call in its hot loop. `ok=FALSE` means the
#'   maps could not be built for this fit (a covariate on a phi0/phi1 theta,
#'   a non-mu-referenced eta, or a raw model covariate in the Hess2/pred
#'   model are all out of v1 scope) -- Phase 4 falls back to the historic
#'   SA-recursion for such fits rather than guessing.
#' @author Matthew L. Fidler
#' @noRd
#' @export
rxUiGet.saemPhi1Inner <- function(x, ...) {
  .ui <- x[[1]]
  if (!.saemGeneralLik(.ui)) return(NULL)
  # saemControl(phi1Hessian=FALSE) is the default (see its own docs -- an
  # ablation check found the Laplace log|H| correction was not what fixed a
  # diverging Gaussian twin, and it can instead dominate/diverge for a
  # heavy-tailed t()/cauchy() endpoint). innerHess2 needs foceiControl(fast
  # =TRUE) to build at all, which is otherwise wasted compile time -- only
  # force it, and only ask .saemPhi1TargetMap to resolve a Hessian-capable
  # map, when the Hessian correction is actually wanted.
  .wantHessian <- isTRUE(tryCatch(
    as.logical(rxode2::rxGetControl(.ui, "phi1Hessian", FALSE)),
    error = function(e) FALSE))
  .origControl <- .ui$control
  if (.wantHessian &&
        !isTRUE(tryCatch(as.logical(rxode2::rxGetControl(.ui, "fast", FALSE)), error = function(e) FALSE))) {
    assign("control", foceiControl(fast = TRUE), envir = .ui)
  }
  .fm <- tryCatch(.ui$focei, error = function(e) NULL)
  if (!is.null(.origControl)) {
    assign("control", .origControl, envir = .ui)
  } else if (exists("control", envir = .ui, inherits = FALSE)) {
    rm(list = "control", envir = .ui)
  }
  if (is.null(.fm)) return(NULL)
  .inner <- if (!is.null(.fm$innerLlik)) .fm$innerLlik else .fm$inner
  if (is.null(.inner)) return(NULL)
  .predNoLhs <- if (!is.null(.fm$predNoLhsLlik)) .fm$predNoLhsLlik else .fm$predNoLhs
  .innerHess2 <- if (.wantHessian) .fm$innerHess2 else NULL
  .map <- .saemPhi1TargetMap(.ui, .innerHess2, .predNoLhs)
  c(list(inner = .inner, innerHess2 = .innerHess2, predNoLhs = .predNoLhs),
    if (is.null(.map)) list(ok = FALSE) else .map)
}
attr(rxUiGet.saemPhi1Inner, "rstudio") <- emptyenv()

#' Build the THETA\[k\]/ETA\[k\] -> SAEM phi-column maps for the phi1 Laplace step
#'
#' `innerHess2`/`predNoLhs` are compiled by FOCEi's own codegen, so their
#' parameter vector is `THETA[k]`/`ETA[k]` (plus `DV`, a general-likelihood
#' model's own event-table covariate -- already supplied the same way SAEM's
#' own model gets it, so it needs no per-row handling here). SAEM's C++ theta
#' step needs, resolved ONCE and consumed as plain integer arrays (no R call
#' in the hot per-bobyqa-evaluation loop): for each `THETA[k]`, whether it is
#' a phi1 (mu-referenced, actively optimized) column, a phi0 (fixed-effect)
#' column read from SAEM's current state, or a genuinely FIXED value; and for
#' each `ETA[k]`, which phi1 column it pairs with.
#'
#' Deliberately conservative for v1: returns `NULL` (the caller falls back to
#' the historic SA-recursion) whenever a shape this map cannot represent
#' shows up -- a covariate on a phi0/phi1 theta, a non-mu-referenced
#' (`ui$nonMuEtas`) estimated eta, an eta not paired to a phi1 theta by
#' mu-referencing, a raw model covariate (anything in the compiled model's
#' parameter vector besides `THETA[.]`/`ETA[.]`/`DV`), or `innerHess2` and
#' `predNoLhs` disagreeing on THETA/ETA parameter order (they are compiled
#' from the same theta/eta declarations, so a mismatch means something else
#' is wrong and this must not guess).
#'
#' @param ui rxode2 ui object
#' @param hess2Mod the (possibly `NULL`) `innerHess2` model
#' @param predMod the `predNoLhs` model (`NULL` only if `hess2Mod` is also `NULL`)
#' @return `NULL`, or a list with `ok=TRUE`, `thetaKind` (integer vector,
#'   length = number of `THETA[k]` slots; -1 fixed, 0 phi0, 1 phi1),
#'   `thetaCol` (0-based column within that kind's own subset),
#'   `thetaFixedVal` (the fixed value when `thetaKind==-1`, ignored
#'   otherwise), and `etaCol` (integer vector, length = number of `ETA[k]`
#'   slots; the 0-based phi1 column each eta pairs with).
#' @author Matthew L. Fidler
#' @noRd
.saemPhi1TargetMap <- function(ui, hess2Mod, predMod) {
  .refMod <- if (!is.null(hess2Mod)) hess2Mod else predMod
  if (is.null(.refMod)) return(NULL)
  .parsH2 <- rxode2::rxParam(.refMod)
  if (!is.null(hess2Mod) && !is.null(predMod)) {
    .thEtH2 <- .parsH2[grepl("^(THETA|ETA)\\[", .parsH2)]
    .thEtPred <- rxode2::rxParam(predMod)
    .thEtPred <- .thEtPred[grepl("^(THETA|ETA)\\[", .thEtPred)]
    if (!identical(.thEtH2, .thEtPred)) return(NULL)
  }
  .other <- .parsH2[!(grepl("^(THETA|ETA)\\[", .parsH2) | .parsH2 == "DV")]
  if (length(.other) > 0) return(NULL)
  if (length(ui$nonMuEtas) > 0) return(NULL)
  .cov <- rxUiGet.saemMuRefCovariateDataFrame(list(ui))
  if (length(.cov$covariateParameter) > 0) return(NULL)

  .iniDf <- ui$iniDf
  .parsAll <- rxUiGet.saemParamsToEstimateCov(list(ui))
  .muRef <- ui$muRefDataFrame
  .isPhi1 <- .parsAll %in% .muRef$theta
  .phi1Names <- .parsAll[.isPhi1]
  .phi0Names <- .parsAll[!.isPhi1]

  .nTheta <- length(grep("^THETA\\[", .parsH2))
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

  .nEta <- length(grep("^ETA\\[", .parsH2))
  .etaCol <- integer(.nEta)
  .etaDiag <- !is.na(.iniDf$neta1) &
    (is.na(.iniDf$neta2) | .iniDf$neta1 == .iniDf$neta2)
  for (.k in seq_len(.nEta)) {
    .nm <- .iniDf$name[.etaDiag & .iniDf$neta1 == .k]
    if (length(.nm) != 1L) return(NULL)
    .thNm <- .muRef$theta[.muRef$eta == .nm]
    if (length(.thNm) != 1L || !(.thNm %in% .phi1Names)) return(NULL)
    .etaCol[.k] <- match(.thNm, .phi1Names) - 1L
  }

  list(ok = TRUE, thetaKind = .thetaKind, thetaCol = .thetaCol,
       thetaFixedVal = .thetaFixedVal, etaCol = .etaCol)
}
