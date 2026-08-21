#ifndef __TRUSTHESSIANUPDATE_H__
#define __TRUSTHESSIANUPDATE_H__
#if defined(__cplusplus)

// Shared hessianMethod= codes for every RcppTrust-backed Hessian in this
// package (est="trust"'s outer theta problem, src/nlm.cpp; FOCEi's inner eta
// problem for non-normal-endpoint models, src/inner.cpp's calcEtaHessian()).
#define trustHessFd     1
#define trustHessBfgs   2
#define trustHessSr1    3
#define trustHessBofill 4

// Considered and rejected for this Hessian (kept here so neither gets
// re-proposed without new infrastructure this codebase does not have):
//  - Curvature Propagation (Martens, Sutskever & Swersky, "Estimating the
//    Hessian by Back-propagating Curvature", arXiv:1206.6464, ICML 2012):
//    an unbiased stochastic Hessian estimator built by backward-propagating
//    injected noise through a general reverse-mode autodiff computational
//    graph. This codebase's gradients (nlmSolveGrad, lpInner) come from
//    hand-written ODE sensitivity equations compiled by rxode2, not a
//    traceable AD graph over arbitrary ops -- CP needs a general AD engine
//    this package does not have.
//  - O1NumHess (Zhao et al., "O1NumHess: A Fast and Accurate Seminumerical
//    Hessian Algorithm Using Only O(1) Gradients", arXiv:2508.07544, JCTC
//    2025): reduces a molecular Hessian from O(N_atom) gradients to O(1) by
//    exploiting the off-diagonal low-rank (ODLR) property that Hessian
//    blocks coupling two SPATIALLY DISTANT groups of atoms are low rank.
//    That property has no analog in a small (order 2-20) theta/eta vector
//    with no spatial locality, and the plain O(n) FD-of-gradient cost it is
//    built to avoid is not expensive at this size.
//
// Quasi-Newton Hessian update from two consecutive (point, gradient) pairs --
// s = x_k - x_{k-1}, y = g_k - g_{k-1} -- avoiding a fresh per-parameter
// finite-difference step-size derivation every iteration. H is updated in
// place; left unchanged when the secant pair is unreliable (near-zero
// denominator -- the same failure mode as an unstable finite-difference
// step, generalized).
static inline void trustHessianUpdate(int method, arma::mat &H,
                                       const arma::vec &s, const arma::vec &y) {
  arma::vec Hs = H * s;
  if (method == trustHessBfgs) {
    // Damped BFGS: Nocedal & Wright, Numerical Optimization, 2nd ed.
    // (2006), Procedure 18.2 / Eq. 18.16. trust does no line search, so the
    // curvature condition s'y>0 is not guaranteed -- Powell's damping keeps
    // the update positive-definite regardless.
    double sBs = arma::dot(s, Hs);
    if (sBs <= 0) return; // H itself not PD (should not happen after seeding); skip
    double sy = arma::dot(s, y);
    double thetaD = (sy >= 0.2 * sBs) ? 1.0 : (0.8 * sBs / (sBs - sy));
    arma::vec yBar = thetaD * y + (1.0 - thetaD) * Hs;
    double syBar = arma::dot(s, yBar);
    if (syBar <= 1e-10 * arma::norm(s) * arma::norm(yBar)) return; // skip
    H += (yBar * yBar.t()) / syBar - (Hs * Hs.t()) / sBs;
  } else if (method == trustHessSr1) {
    // Symmetric Rank-1: Nocedal & Wright Eq. 6.24, skip safeguard Eq. 6.26
    // (originally Murtagh & Sargent, Comput. J. 13, 185-194, 1970). N&W
    // recommend SR1 specifically for trust-region methods: unlike BFGS it
    // is not forced positive-definite, so it can represent the indefinite
    // curvature a trust-region step can legitimately encounter.
    arma::vec r = y - Hs;
    double denom = arma::dot(s, r);
    double rNorm = arma::norm(r), sNorm = arma::norm(s);
    if (std::fabs(denom) >= 1e-8 * sNorm * rNorm && rNorm > 0)
      H += (r * r.t()) / denom;
    // else: unreliable secant pair (e.g. a tiny reject-then-shrink step) --
    // skip rather than corrupt the running Hessian.
  } else if (method == trustHessBofill) {
    // Bofill, J. Comput. Chem. 15, 1-11 (1994): a phi-weighted blend of SR1
    // (Murtagh-Sargent) and Powell-Symmetric-Broyden (PSB), the standard
    // Berny/transition-state-search Hessian update -- phi (in [0,1]) is the
    // squared cosine between s and the SR1 correction r, so the update
    // smoothly favors whichever of SR1/PSB the data actually supports.
    arma::vec r = y - Hs;
    double ss = arma::dot(s, s), rr = arma::dot(r, r), sr = arma::dot(s, r);
    if (ss <= 0 || rr <= 0) return; // skip
    double phi = (sr * sr) / (ss * rr);
    arma::mat Hms = H + (r * r.t()) / sr;
    arma::mat Hpsb = H + (r * s.t() + s * r.t()) / ss - (sr * (s * s.t())) / (ss * ss);
    if (std::fabs(sr) < 1e-8 * std::sqrt(ss) * std::sqrt(rr)) {
      // SR1 term undefined (sr ~ 0) -- fall back to pure PSB.
      H = Hpsb;
    } else {
      H = phi * Hms + (1.0 - phi) * Hpsb;
    }
  }
}

#endif
#endif
