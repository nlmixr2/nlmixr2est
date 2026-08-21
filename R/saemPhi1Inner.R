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
#' `ui$focei$innerHess2` (the exact 2nd-order eta-Hessian model) requires
#' `foceiControl(fast=TRUE)` and `.foceiLLGradInScope(ui)`; this accessor
#' temporarily installs `foceiControl(fast=TRUE)` on the ui if it is not
#' already FOCEi-shaped (SAEM's own control has no `fast=` field), builds
#' `ui$focei` under that control, then restores the ui's original control --
#' `ui$focei` is a cached accessor, so the built model list persists.
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
#'   for this model shape), and `predNoLhs` (the bare, no-sensitivity
#'   prediction-only FD-fallback model, always non-NULL).
#' @author Matthew L. Fidler
#' @noRd
#' @export
rxUiGet.saemPhi1Inner <- function(x, ...) {
  .ui <- x[[1]]
  if (!.saemGeneralLik(.ui)) return(NULL)
  .origControl <- .ui$control
  if (!isTRUE(tryCatch(as.logical(rxode2::rxGetControl(.ui, "fast", FALSE)), error = function(e) FALSE))) {
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
  list(inner = .inner, innerHess2 = .fm$innerHess2, predNoLhs = .predNoLhs)
}
attr(rxUiGet.saemPhi1Inner, "rstudio") <- emptyenv()
