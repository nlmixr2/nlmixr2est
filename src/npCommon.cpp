// Shared numerical routines for the nonparametric engines.  See npCommon.h.
//
// npBurke is a direct port of Burke's interior-point method as implemented in
// the Pmetrics/pmcore Rust reference (routines/estimation/ipm.rs, crate pmcore
// 0.25.2), which itself follows Yamada 2021 Appendix A.  Kept close to the
// reference so the golden-fixture tests can compare element-wise.
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <boost/random/sobol.hpp>
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

arma::mat npSobolGrid(int n, const arma::vec& lower, const arma::vec& upper) {
  int d = (int)lower.n_elem;
  arma::mat grid(n, d);
  if (n <= 0 || d <= 0) { grid.set_size(std::max(n, 0), std::max(d, 0)); return grid; }
  boost::random::sobol eng((std::size_t)d);
  for (int k = 0; k < n; ++k) {
    for (int j = 0; j < d; ++j) {
      uint64_t v = (uint64_t)eng();
      double u = std::ldexp((double)(v >> 11) + 0.5, -53);  // in (0,1)
      grid(k, j) = lower[j] + u * (upper[j] - lower[j]);
    }
  }
  return grid;
}

arma::uvec npCondenseWeights(const arma::vec& lambda, double ratio) {
  double thr = lambda.max() * ratio;
  std::vector<arma::uword> keep;
  keep.reserve(lambda.n_elem);
  for (arma::uword k = 0; k < lambda.n_elem; ++k) {
    if (lambda[k] > thr) keep.push_back(k);
  }
  return arma::uvec(keep);
}

arma::uvec npCondenseQR(const arma::mat& psi, double tol) {
  const int nsub = (int)psi.n_rows;
  const int npoint = (int)psi.n_cols;
  if (nsub == 0 || npoint == 0) return arma::uvec();
  // Row-normalize psi (each subject sums to 1); a zero-sum row is an error.
  Eigen::MatrixXd m(nsub, npoint);
  for (int i = 0; i < nsub; ++i) {
    double s = 0.0;
    for (int j = 0; j < npoint; ++j) s += psi(i, j);
    if (s == 0.0) Rcpp::stop("npCondenseQR: psi row %d sums to zero", i + 1);
    for (int j = 0; j < npoint; ++j) m(i, j) = psi(i, j) / s;
  }
  // Column-pivoted (rank-revealing) QR: columns ordered by decreasing pivot.
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(m);
  Eigen::MatrixXd R = qr.matrixR().triangularView<Eigen::Upper>();
  Eigen::VectorXi perm = qr.colsPermutation().indices();
  int keepN = std::min(nsub, npoint);
  std::vector<arma::uword> keep;
  keep.reserve(keepN);
  for (int i = 0; i < keepN; ++i) {
    double colNorm = R.col(i).norm();
    double diag = R(i, i);
    if (colNorm > 0.0 && std::fabs(diag / colNorm) >= tol) {
      keep.push_back((arma::uword)perm[i]);
    }
  }
  std::sort(keep.begin(), keep.end());
  return arma::uvec(keep);
}
