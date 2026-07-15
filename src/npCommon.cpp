// Shared numerical routines for the nonparametric engines.  See npCommon.h.
//
// npBurke is a direct port of Burke's interior-point method as implemented in
// the Pmetrics/pmcore Rust reference (routines/estimation/ipm.rs, crate pmcore
// 0.25.2), which itself follows Yamada 2021 Appendix A.  Kept close to the
// reference so the golden-fixture tests can compare element-wise.
#include <RcppArmadillo.h>
#include "npCommon.h"

using namespace arma;

arma::vec npBurke(const arma::mat& psiIn, double* obj) {
  mat psi = psiIn;
  // Coerce negatives to non-negative; non-finite entries are an error.
  for (uword j = 0; j < psi.n_elem; ++j) {
    double v = psi[j];
    if (!std::isfinite(v)) {
      Rcpp::stop("npBurke: the likelihood (psi) matrix must have finite entries");
    }
    psi[j] = std::fabs(v);
  }

  const uword nSub = psi.n_rows;
  const uword nPoint = psi.n_cols;
  const double eps = 1e-8;

  vec ecol(nPoint, fill::ones);   // sums over support points
  vec erow(nSub, fill::ones);     // sums over subjects

  vec lam(nPoint, fill::ones);
  vec plam = psi * ecol;          // per-subject mixture density (row sums)
  double sig = 0.0;

  vec w = 1.0 / plam;             // n_sub
  vec ptw = psi.t() * w;          // n_point

  double shrink = 2.0 * ptw.max();
  lam *= shrink;
  plam *= shrink;
  w /= shrink;
  ptw /= shrink;

  vec y = ecol - ptw;                 // n_point
  vec r = erow - (w % plam);          // n_sub
  double normR = norm(r, "inf");

  double sumLogPlam = accu(log(plam));
  double sumLogW = accu(log(w));
  double gap = std::fabs(sumLogW + sumLogPlam) / (1.0 + std::fabs(sumLogPlam));
  double mu = dot(lam, y) / (double)nPoint;

  while (mu > eps || normR > eps || gap > eps) {
    double smu = sig * mu;
    vec inner = lam / y;              // n_point
    vec wPlam = plam / w;            // n_sub

    // H = psi * diag(inner) * psi^T + diag(wPlam)   (n_sub x n_sub, SPD)
    mat psiInner = psi;
    psiInner.each_row() %= inner.t();
    mat H = psiInner * psi.t();
    H.diag() += wPlam;

    mat L;
    if (!chol(L, H, "lower")) {
      Rcpp::stop("npBurke: Cholesky decomposition failed (Newton matrix not "
                 "positive definite); usually model misspecification or numerical "
                 "issues");
    }

    vec smuyinv = smu * (ecol / y);              // n_point
    vec rhsdw = (erow / w) - psi * smuyinv;      // n_sub
    // Solve H dw = rhsdw via the Cholesky factor (H = L L^T).
    vec z = solve(trimatl(L), rhsdw);
    vec dw = solve(trimatu(L.t()), z);           // n_sub

    vec dy = -(psi.t() * dw);                    // n_point
    vec dlam = smuyinv - lam - (inner % dy);     // n_point

    double minRatioDlam = (dlam / lam).min();
    double alfpri = -1.0 / std::min(minRatioDlam, -0.5);
    alfpri = std::min(0.99995 * alfpri, 1.0);

    double minRatioDy = (dy / y).min();
    double minRatioDw = (dw / w).min();
    double alfdual = -1.0 / std::min(minRatioDy, -0.5);
    alfdual = std::min(alfdual, -1.0 / std::min(minRatioDw, -0.5));
    alfdual = std::min(0.99995 * alfdual, 1.0);

    lam += alfpri * dlam;
    w += alfdual * dw;
    y += alfdual * dy;

    mu = dot(lam, y) / (double)nPoint;
    plam = psi * lam;
    r = erow - (w % plam);
    ptw -= alfdual * dy;

    normR = norm(r, "inf");
    sumLogPlam = accu(log(plam));
    sumLogW = accu(log(w));
    gap = std::fabs(sumLogW + sumLogPlam) / (1.0 + std::fabs(sumLogPlam));

    if (mu < eps && normR > eps) {
      sig = 1.0;
    } else {
      double c1 = (1.0 - alfpri) * (1.0 - alfpri);
      double c2 = (1.0 - alfdual) * (1.0 - alfdual);
      double c3 = (normR - mu) / (normR + 100.0 * mu);
      sig = std::min(std::max(std::max(c1, c2), c3), 0.3);
    }
  }

  lam /= (double)nSub;
  if (obj != nullptr) {
    *obj = accu(log(psi * lam));
  }
  lam /= accu(lam);   // normalize to a probability vector
  return lam;
}
