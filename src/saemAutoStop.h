// Noise-aware windowed early stopping for SAEM's EM phase (#1411).
#ifndef NLMIXR2EST_SAEM_AUTO_STOP_H
#define NLMIXR2EST_SAEM_AUTO_STOP_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

struct saemAutoStop {
  int window = 50;
  int consec = 4;
  double rtol = 1e-3;
  double atol = 1e-4;
  int count = 0;      // consecutive stationary iterations
  int n = 0;          // iterates pushed
  arma::mat ring;     // window x p ring buffer of the iterates
  arma::vec drift, mcse, tol;

  void init(int w, int c, double r, double a, int p) {
    window = std::max(w, 4);
    consec = std::max(c, 1);
    rtol = r;
    atol = a;
    count = 0;
    n = 0;
    ring.set_size(window, p);
    drift.set_size(p); drift.fill(NA_REAL);
    mcse.set_size(p);  mcse.fill(NA_REAL);
    tol.set_size(p);   tol.fill(NA_REAL);
  }

  // Add one iterate; true once every coordinate has been stationary for
  // `consec` consecutive full windows.
  bool push(const arma::vec &x) {
    if (x.n_elem != ring.n_cols) return false;
    ring.row(n % window) = x.t();
    n++;
    if (n < window) return false;
    count = stationary() ? count + 1 : 0;
    return count >= consec;
  }

  // Compare the older and newer half means of the window.  The Monte-Carlo
  // error comes from the window detrended by a least-squares line, so a trend
  // shows up in the drift and never inflates its own tolerance; lag-1
  // autocorrelation of the residuals discounts the effective sample size.
  bool stationary() {
    const int h = window / 2, m = 2 * h, p = (int)ring.n_cols;
    arma::mat w(m, p);
    for (int i = 0; i < m; ++i) w.row(i) = ring.row((n - m + i) % window);
    if (!w.is_finite()) {
      drift.fill(NA_REAL); mcse.fill(NA_REAL); tol.fill(NA_REAL);
      return false;
    }
    arma::vec t = arma::linspace<arma::vec>(-(m - 1) / 2.0, (m - 1) / 2.0, m);
    const double stt = arma::dot(t, t);
    bool ok = true;
    for (int j = 0; j < p; ++j) {
      arma::vec y = w.col(j);
      const double mA = arma::mean(y.head(h)), mB = arma::mean(y.tail(h));
      const double yb = arma::mean(y);
      arma::vec r = y - yb - (arma::dot(t, y - yb) / stt) * t;
      const double ss = arma::dot(r, r);
      double rho = 0.0;
      if (ss > 0.0) rho = arma::dot(r.head(m - 1), r.tail(m - 1)) / ss;
      rho = std::min(std::max(rho, 0.0), 0.9);
      const double neff = std::max(1.0, h * (1.0 - rho) / (1.0 + rho));
      const double s2 = ss / (m - 2);
      drift(j) = std::fabs(mB - mA);
      mcse(j) = std::sqrt(2.0 * s2 / neff);
      tol(j) = std::max({atol, rtol * std::max(1.0, std::fabs(mA)), 2.0 * mcse(j)});
      if (!(drift(j) <= tol(j))) ok = false;
    }
    return ok;
  }
};

// a plain R vector rather than arma's one-column matrix
inline Rcpp::NumericVector saemAutoStopVec(const arma::vec &v) {
  return Rcpp::NumericVector(v.begin(), v.end());
}

#endif // NLMIXR2EST_SAEM_AUTO_STOP_H
