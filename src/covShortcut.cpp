#define STRICT_R_HEADER
#include "covShortcut.h"
#include <cmath>

// Directions checked: rows 0..3 of the Sylvester-Hadamard matrix of order
// 2^ceil(log2 n), restricted to the n coordinates.  Row r's sign at column j is
// (-1)^popcount(r & j); row 0 is all +1, and together the rows weigh every
// coordinate pair with both signs.
static const int covShortcutNdir = 4;

// A direction is explained when measured and predicted agree to 1%.  One wrong
// correlation rho_ij moves the second difference by about 2*rho_ij/n of it, so this
// catches an error of about 0.035 in a pair of 7 parameters.  Finite-difference noise
// larger than this makes the check fail, which only costs measuring the off-diagonals.
static const double covShortcutRelTol = 0.01;

static inline double covShortcutSign(int r, int j) {
  return (__builtin_popcount((unsigned)(r & j)) & 1) ? -1.0 : 1.0;
}

// 5-point second difference of the objective along v about x (f0 = f(x)), in the
// objective's own units: f''(x)[v, v]; NA_REAL when a probe is not finite.
static double covShortcutD2(FdHessObj &obj, double *x, int n, double f0,
                            const std::vector<double> &v) {
  static const double off[4] = {2, 1, -1, -2};
  std::vector<double> x0(x, x + n);
  double fv[4];
  for (int k = 0; k < 4; ++k) {
    for (int j = 0; j < n; ++j) x[j] = x0[j] + off[k]*v[j];
    fv[k] = obj.f(x);
    if (!R_FINITE(fv[k])) {
      std::copy(x0.begin(), x0.end(), x);
      return NA_REAL;
    }
  }
  std::copy(x0.begin(), x0.end(), x);
  return (-fv[0] + 16*fv[1] - 30*f0 + 16*fv[2] - fv[3])/12.0;
}

// Direction v_j = s_j*h[j] (the coordinate steps, signed by a Hadamard row); H (fac
// times the objective's Hessian) predicts v' H v / fac for its second difference.  A
// direction is explained when they agree to covShortcutRelTol.
CovShortcutResult covShortcutVerify(FdHessObj &obj, double *x, int n, double f0,
                                    const double *h, const arma::mat &H, double fac) {
  CovShortcutResult res;
  int nDir = n < 2 ? 0 : std::min(covShortcutNdir, 1 << (int)std::ceil(std::log2((double)n)));
  res.checks.set_size(nDir, 3);
  std::vector<double> v(n);
  res.status = 1;
  for (int r = 0; r < nDir; ++r) {
    for (int j = 0; j < n; ++j) v[j] = covShortcutSign(r, j)*h[j];
    arma::vec va(v.data(), n, false, true);
    double pred = arma::as_scalar(va.t() * H * va)/fac;
    double meas = covShortcutD2(obj, x, n, f0, v);
    double allow = covShortcutRelTol*std::fabs(pred);
    res.checks(r, 0) = meas;
    res.checks(r, 1) = pred;
    res.checks(r, 2) = allow;
    if (!R_FINITE(meas)) {
      res.status = -1;
      res.checks.resize(r + 1, 3);
      break;
    }
    if (std::fabs(meas - pred) > allow) {
      res.status = 0;
      res.checks.resize(r + 1, 3);
      break;
    }
  }
  obj.restore(x);
  return res;
}
