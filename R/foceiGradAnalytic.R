# Analytic FOCE/FOCEI outer (population) gradient from Almquist (2015) sensitivity
# equations -- the first-derivative precursor of the analytic observed-information
# R-matrix in foceiCovAnalytic.R.  Enabled by foceiControl(fast=TRUE); called from
# the C++ outer-gradient hook (numericGrad -> analyticOuterGrad).  Returns a
# length-npars gradient of the objective (-2*logLik) in the model-theta scale, or
# NULL to fall back to the finite-difference gradient.
#
# The outer gradient needs at most SECOND-order state sensitivities (Almquist Eqs
# 38-40) -- one order less than the covariance R-matrix -- so it reuses the
# direction set, error machinery, and per-subject H/Ht/N/dHtD/etaP assembly from
# foceiCovAnalytic.R without the 3rd-order (Ath) tier.

#' Analytic outer gradient of the FOCE/FOCEI objective (Almquist 2015 Eq 23).
#'
#' Called from C++ (`analyticOuterGrad`) when `foceiControl(fast=TRUE)`.  The fit
#' environment `e` carries the current theta/eta state; the C++ caller has already
#' evaluated the objective at the current theta so the inner solutions (eta*) are
#' current.
#'
#' @param e focei fit environment
#' @return numeric gradient vector (length npars, model-theta scale), or `NULL`
#'   to fall back to the finite-difference gradient
#' @noRd
.foceiCalcGradAnalytic <- function(e) {
  # TODO(phase-1): assemble the analytic gradient.  Stub returns NULL so the
  # C++ hook falls back to the finite-difference gradient while the plumbing is
  # validated.
  NULL
}
