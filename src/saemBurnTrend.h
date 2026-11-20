// SAEM burn-in convergence test (#1401): batch-means slope t-test per column.
#ifndef NLMIXR2EST_SAEM_BURN_TREND_H
#define NLMIXR2EST_SAEM_BURN_TREND_H

#include "armahead.h"
#include <algorithm>
#include <cmath>
#include <limits>

// 0 = constant, 1 = varying and finite, 2 = non-constant with a non-finite value.
static inline int saemBurnTrendColKind(const double *col, arma::uword n) {
  bool constant = true, finite = true;
  for (arma::uword i = 0; i < n; ++i) {
    if (!std::isfinite(col[i])) finite = false;
    if (!(col[i] == col[0])) constant = false;
  }
  if (constant) return 0;
  return finite ? 1 : 2;
}

// Fills `y` with the `nBatch` means of `m` consecutive rows of `col`; returns
// the largest |mean|.
static inline double saemBurnTrendBatchMeans(const double *col, arma::uword m,
                                             int nBatch, arma::vec &y) {
  double scale = 0.0;
  for (int b = 0; b < nBatch; ++b) {
    double s = 0.0;
    for (arma::uword i = 0; i < m; ++i) s += col[b * m + i];
    y(b) = s / static_cast<double>(m);
    scale = std::max(scale, std::fabs(y(b)));
  }
  return scale;
}

// |t| of the slope of the batch means `y`; 0 when there is no drift to test.
static inline double saemBurnTrendSlopeT(const arma::vec &y, int nBatch,
                                         double sxx, double scale) {
  const double eps = std::numeric_limits<double>::epsilon();
  const double xbar = (nBatch + 1) / 2.0;
  double ybar = 0.0;
  for (int b = 0; b < nBatch; ++b) ybar += y(b);
  ybar /= nBatch;
  double sxy = 0.0;
  for (int b = 0; b < nBatch; ++b) sxy += (b + 1 - xbar) * (y(b) - ybar);
  const double slope = sxy / sxx;
  // a drift below rounding of the values themselves is no trend
  if (std::fabs(slope) * (nBatch - 1) <= 64.0 * eps * scale) return 0.0;
  double rss = 0.0;
  for (int b = 0; b < nBatch; ++b) {
    const double e = y(b) - ybar - slope * (b + 1 - xbar);
    rss += e * e;
  }
  const double se = std::sqrt(rss / (nBatch - 2) / sxx);
  return se > 0.0 ? std::fabs(slope) / se : R_PosInf;
}

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
  arma::vec t(nc, arma::fill::zeros);
  const arma::uword m = nBatch < 3 ? 0 : x.n_rows / static_cast<arma::uword>(nBatch);
  if (m > 0) {
    const arma::uword n = m * static_cast<arma::uword>(nBatch);
    const double xbar = (nBatch + 1) / 2.0;
    double sxx = 0.0;
    for (int b = 1; b <= nBatch; ++b) sxx += (b - xbar) * (b - xbar);
    arma::vec y(nBatch);
    for (arma::uword c = 0; c < nc; ++c) {
      const double *col = x.colptr(c) + (x.n_rows - n);
      const int kind = saemBurnTrendColKind(col, n);
      if (kind == 2) {
        t(c) = R_PosInf;
        ret(c) = 1;
      } else if (kind == 1) {
        const double scale = saemBurnTrendBatchMeans(col, m, nBatch, y);
        t(c) = saemBurnTrendSlopeT(y, nBatch, sxx, scale);
        ret(c) = t(c) > tcrit ? 1 : 0;
      }
    }
  }
  if (tstat) *tstat = t;
  return ret;
}

#endif
