// SAEM burn-in convergence test (#1401): batch-means slope t-test per column.
#ifndef NLMIXR2EST_SAEM_BURN_TREND_H
#define NLMIXR2EST_SAEM_BURN_TREND_H

#include "armahead.h"
#include <cmath>
#include <limits>

// Flags each column of `x` (iterations x parameters) that still trends: the
// trailing rows are cut into `nBatch` equal batches and the slope of the batch
// means is t-tested (nBatch - 2 df) against `tcrit`.  Batch means, not raw
// iterates, keep the false-trend rate near nominal under autocorrelation.
// A constant column never trends; one with other non-finite values always does.
// `tstat` (optional) receives each column's |t| (Inf for a flagged degenerate one).
static inline arma::uvec saemBurnTrend(const arma::mat &x, int nBatch, double tcrit,
                                       arma::vec *tstat = nullptr) {
  const arma::uword nc = x.n_cols;
  arma::uvec ret(nc, arma::fill::zeros);
  if (tstat) tstat->zeros(nc);
  if (nBatch < 3) return ret;
  const arma::uword m = x.n_rows / (arma::uword)nBatch;
  if (m == 0) return ret;
  const arma::uword r0 = x.n_rows - m * (arma::uword)nBatch;
  const double eps = std::numeric_limits<double>::epsilon();
  const double xbar = (nBatch + 1) / 2.0;
  double sxx = 0.0;
  for (int b = 1; b <= nBatch; ++b) sxx += (b - xbar) * (b - xbar);
  arma::vec y(nBatch);
  for (arma::uword c = 0; c < nc; ++c) {
    const double *col = x.colptr(c) + r0;
    const arma::uword n = m * (arma::uword)nBatch;
    bool constant = true, finite = true;
    for (arma::uword i = 0; i < n; ++i) {
      if (!std::isfinite(col[i])) finite = false;
      if (!(col[i] == col[0])) constant = false;
    }
    if (constant) continue;
    if (!finite) {
      ret(c) = 1;
      if (tstat) (*tstat)(c) = R_PosInf;
      continue;
    }
    double ybar = 0.0, scale = 0.0;
    for (int b = 0; b < nBatch; ++b) {
      double s = 0.0;
      for (arma::uword i = 0; i < m; ++i) s += col[b * m + i];
      y(b) = s / (double)m;
      ybar += y(b);
      scale = std::max(scale, std::fabs(y(b)));
    }
    ybar /= nBatch;
    double sxy = 0.0;
    for (int b = 0; b < nBatch; ++b) sxy += (b + 1 - xbar) * (y(b) - ybar);
    const double slope = sxy / sxx;
    // a drift below rounding of the values themselves is no trend
    if (std::fabs(slope) * (nBatch - 1) <= 64.0 * eps * scale) continue;
    double rss = 0.0;
    for (int b = 0; b < nBatch; ++b) {
      const double e = y(b) - ybar - slope * (b + 1 - xbar);
      rss += e * e;
    }
    const double se = std::sqrt(rss / (nBatch - 2) / sxx);
    const double t = se > 0.0 ? std::fabs(slope) / se : R_PosInf;
    if (tstat) (*tstat)(c) = t;
    if (t > tcrit) ret(c) = 1;
  }
  return ret;
}

#endif
