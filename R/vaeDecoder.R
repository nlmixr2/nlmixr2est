# vaeDecoder.R -- VAE decoder: structural prediction f, residual variance R, and
# their eta-sensitivities df/deta, dR/deta, reusing the focei analytic
# sensitivity machinery (R/foceiCovAnalytic.R). The augmented rxode2 model emits
# rx_predf_ (=f), rx_f1_<eta> (=df/deta), rx_rvarf_ (=R), rx_rvar1_<eta>
# (=dR/deta) as solve columns, so no ODE-sensitivity code is re-implemented here.
#
# The decoder's contribution to the ELBO is the data term p(x|z) and its
# gradient d(p_x_z)/d(eta) -- which is exactly the encoder's dLoss/dz upstream
# (z = z_pop + eta, so d/dz = d/deta). It reuses focei's rho(f,R) partials.

#' Build the VAE decoder augmented model (directions = etas only)
#'
#' @param ui rxode2 ui object
#' @return list from .foceiAnalyticAugModelDirs (augMod + metadata), or NULL
#' @noRd
.vaeDecoderModel <- function(ui) {
  .neta <- length(ui$iniDf$name[!is.na(ui$iniDf$neta1) & ui$iniDf$neta1 == ui$iniDf$neta2])
  .dirs <- paste0("ETA_", seq_len(.neta), "_")
  .foceiAnalyticAugModelDirs(ui, .dirs)
}

#' Solve one subject's decoder model, returning f, R, df/deta, dR/deta per obs
#'
#' @param am augmented model from .vaeDecoderModel
#' @param th named theta vector (THETA_1_ = value, ...)
#' @param eta numeric eta vector for the subject
#' @param ev per-subject event table (dosing + covariates)
#' @param times observation times
#' @param tol solve tolerance
#' @param maxRecalc max tolerance-relaxation retries on a failed/non-finite solve
#' @param recalcFactor multiplicative tolerance loosening per retry
#' @param fdFallback when TRUE, if the analytic eta-sensitivities come back
#'   non-finite, replace them with central finite differences of the (robust)
#'   prediction/variance columns
#' @return list(f, R, a = df/deta [nObs x neta], aR = dR/deta [nObs x neta]) or NULL
#' @noRd
.vaeDecoderSolveSubject <- function(am, th, eta, ev, times, tol = 1e-10,
                                    maxRecalc = 5L, recalcFactor = 10^(0.5),
                                    fdFallback = TRUE) {
  .etav <- am$dirs
  .solve <- function(e, t) .foceiAnalyticSolveFA(am, c(th, setNames(e, .etav)), ev, times, tol = t)
  .allFin <- function(E) !is.null(E) && all(is.finite(E$f)) && all(is.finite(E$a)) &&
    all(is.finite(E$R)) && all(is.finite(E$aR))
  ## tolerance relaxation: retry with progressively looser tol (inner.cpp-style)
  .t <- tol; E <- NULL
  for (.attempt in 0:maxRecalc) {
    E <- .solve(eta, .t)
    if (.allFin(E)) return(list(f = E$f, R = E$R, a = E$a, aR = E$aR))
    .t <- .t * recalcFactor
  }
  ## can't get a usable prediction even relaxed -> give up
  if (is.null(E) || !all(is.finite(E$f)) || !all(is.finite(E$R))) return(NULL)
  if (!fdFallback) return(NULL)
  ## finite-difference fallback for the eta-sensitivities (the analytic sens
  ## states blew up but f/R from the same solve are finite)
  .neta <- length(eta); .no <- length(E$f); .h <- 1e-4
  a <- matrix(0, .no, .neta); aR <- matrix(0, .no, .neta)
  for (k in seq_len(.neta)) {
    ep <- eta; ep[k] <- ep[k] + .h; em <- eta; em[k] <- em[k] - .h
    Ep <- .solve(ep, .t); Em <- .solve(em, .t)
    if (is.null(Ep) || is.null(Em) || !all(is.finite(Ep$f)) || !all(is.finite(Em$f))) return(NULL)
    a[, k] <- (Ep$f - Em$f) / (2 * .h)
    aR[, k] <- (Ep$R - Em$R) / (2 * .h)
  }
  list(f = E$f, R = E$R, a = a, aR = aR)
}

#' ELBO data term p(x|z) and its eta-gradient for one subject
#'
#' p_x_z = sum_obs[ 0.5*res^2/R + 0.5*log(R) + 0.5*log(2 pi) ], res = y - f.
#' d(p_x_z)/d(eta_k) = sum_obs[ rf*df/deta_k + rR*dR/deta_k ],
#' rf = -res/R, rR = 0.5*(1/R - res^2/R^2)  (the focei rho(f,R) partials).
#'
#' @param E list from .vaeDecoderSolveSubject
#' @param y observed DV vector (on the rx_pred_ scale for TBS models)
#' @return list(pxz scalar, gEta = d(p_x_z)/d(eta) length neta)
#' @noRd
.vaeDecoderPxz <- function(E, y) {
  .res <- y - E$f
  .R <- E$R
  .rf <- -.res / .R
  .rR <- 0.5 * (1 / .R - .res^2 / .R^2)
  .pxz <- sum(0.5 * .res^2 / .R + 0.5 * log(.R) + 0.5 * log(2 * pi))
  .neta <- ncol(E$a)
  .gEta <- vapply(seq_len(.neta), function(k) sum(.rf * E$a[, k] + .rR * E$aR[, k]),
                  numeric(1))
  list(pxz = .pxz, gEta = .gEta)
}
