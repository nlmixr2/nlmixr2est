// Maximum-likelihood M-step for a declared random-effect distribution.
//
// In the (y, eta) augmentation the complete-data likelihood factors as
//   log p(y | eta) + log p(eta | theta_dist)
// and theta_dist appears ONLY in the second term, so its M-step is a pure
// distribution fit to the etas -- no data term, no ODE solve -- exactly as the
// residual step fits accumulated residuals rather than re-solving.
// rxEtaDistExpand() breaks that by rewriting eta = Q(phi(z)) with z ~ N(0,1),
// which moves theta_dist into the DATA likelihood and lands it in a
// derivative-free search over ODE solves.  EM lets the augmentation be chosen
// freely: sample in z-space, take this M-step in eta-space.
//
// Not circular: prior draws of z would make Q(phi(z)) exactly family(theta_old)
// and return theta_old, but these are POSTERIOR draws -- which is why their
// measured spread is below 1 rather than equal to it -- so they carry data
// information.
//
// Estimated in the family's NATIVE parameters: the objective and the simplex
// are both C++, with no R on the hot loop.  Positive parameters are optimized
// on the log scale because nelder_fn is unbounded.
#include "etaDistFam.h"
#include <algorithm>

typedef void (*fn_ptr) (double *, double *);
extern "C" void nelder_fn(fn_ptr func, int n, double *start, double *step,
                          int itmax, double ftol_rel, double rcoef, double ecoef,
                          double ccoef, int *iconv, int *it, int *nfcall,
                          double *ynewlo, double *xmin, int *iprint);

// nelder_fn takes a plain function pointer, so the objective's data has to be
// reachable at file scope.  Both estimators call the M-step from serial code
// (saem between MCMC sweeps, imp between E-steps), so this is not shared across
// threads.
static std::vector<double> gEtaDistVals;   // etas being fit
static std::vector<double> gEtaDistW;      // matching weights (empty => unit)
static int gEtaDistFam = 0;
static int gEtaDistNa = 0;
static int gEtaDistPos = 0;

static inline void gEtaDistUnpack(const double *p, double *a) {
  for (int i = 0; i < gEtaDistNa; ++i) {
    a[i] = (gEtaDistPos & (1 << i)) ? std::exp(p[i]) : p[i];
  }
}

static double gEtaDistObj(const double *p) {
  double a[4];
  gEtaDistUnpack(p, a);
  for (int i = 0; i < gEtaDistNa; ++i) if (!std::isfinite(a[i])) return 1e300;
  double nll = 0.0;
  const size_t n = gEtaDistVals.size();
  const bool wtd = !gEtaDistW.empty();
  for (size_t i = 0; i < n; ++i) {
    double wi = wtd ? gEtaDistW[i] : 1.0;
    if (wi == 0.0) continue;
    double l = rxEtaDistLogD(gEtaDistFam, gEtaDistVals[i], a);
    if (!std::isfinite(l)) return 1e300;
    nll -= wi*l;
  }
  return std::isfinite(nll) ? nll : 1e300;
}
static void gEtaDistNmFn(double *p, double *fx) { *fx = gEtaDistObj(p); }

bool rxEtaDistMleW(int fam, const std::vector<double> &vals,
                   const std::vector<double> *w, double *a0) {
  int na = rxEtaDistNarg(fam);
  if (na <= 0 || vals.size() < 2) return false;
  gEtaDistFam = fam; gEtaDistNa = na; gEtaDistPos = rxEtaDistPosMask(fam);
  gEtaDistVals = vals;
  gEtaDistW.clear();
  if (w != nullptr) {
    if (w->size() != vals.size()) return false;
    // A weight vector that sums to nothing carries no information; refusing is
    // the honest answer, and leaves the caller's parameters where they were.
    double sw = 0.0;
    for (size_t i = 0; i < w->size(); ++i) {
      double wi = (*w)[i];
      if (!std::isfinite(wi) || wi < 0.0) return false;
      sw += wi;
    }
    if (!(sw > 0.0)) return false;
    gEtaDistW = *w;
  }
  std::vector<double> st(na), stp(na), xm(na);
  for (int i = 0; i < na; ++i) {
    double v = (gEtaDistPos & (1 << i)) ? std::log(a0[i]) : a0[i];
    if (!std::isfinite(v)) return false;
    st[i] = v; xm[i] = v;
    // nelder_fn derives nothing from the start, so give every coordinate a
    // usable step even when it starts at zero
    stp[i] = (std::fabs(v) > 1e-8) ? 0.1*std::fabs(v) : 0.1;
  }
  int iconv, it, nfcall, iprint = 0;
  double ynewlo;
  nelder_fn(gEtaDistNmFn, na, st.data(), stp.data(), 200*na, 1e-8,
            1.0, 2.0, 0.5, &iconv, &it, &nfcall, &ynewlo, xm.data(), &iprint);
  if (!std::isfinite(ynewlo) || ynewlo >= 1e300) return false;
  double a[4];
  gEtaDistUnpack(xm.data(), a);
  for (int i = 0; i < na; ++i) {
    if (!std::isfinite(a[i])) return false;
    a0[i] = a[i];
  }
  return true;
}

bool rxEtaDistMle(int fam, const std::vector<double> &vals, double *a0) {
  return rxEtaDistMleW(fam, vals, nullptr, a0);
}

bool rxEtaDistSpreadOk(const std::vector<double> &w, double lo, double hi,
                       double *sdOut, const std::vector<double> *wt) {
  if (sdOut != nullptr) *sdOut = NA_REAL;
  if (wt != nullptr && wt->size() < w.size()) return false;
  size_t n = 0;
  double sw = 0.0, m = 0.0;
  for (size_t i = 0; i < w.size(); ++i) {
    if (!std::isfinite(w[i])) continue;
    double wi = (wt == nullptr) ? 1.0 : (*wt)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    m += wi*w[i]; sw += wi; n++;
  }
  if (n < 2 || !(sw > 0.0)) return false;
  m /= sw;
  double v = 0.0;
  for (size_t i = 0; i < w.size(); ++i) {
    if (!std::isfinite(w[i])) continue;
    double wi = (wt == nullptr) ? 1.0 : (*wt)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    double d = w[i] - m; v += wi*d*d;
  }
  // reliability weights: the unbiased normalizer is sum(w) - sum(w^2)/sum(w),
  // which reduces to (n - 1) when every weight is 1
  double sw2 = 0.0;
  for (size_t i = 0; i < w.size(); ++i) {
    if (!std::isfinite(w[i])) continue;
    double wi = (wt == nullptr) ? 1.0 : (*wt)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    sw2 += wi*wi;
  }
  double den = sw - sw2/sw;
  if (!(den > 0.0)) den = sw;
  v /= den;
  double sd = std::sqrt(v > 0.0 ? v : 0.0);
  if (sdOut != nullptr) *sdOut = sd;
  return std::isfinite(sd) && sd >= lo && sd <= hi;
}

// Closed-form M-step for a Gaussian copula's correlation.
//
// The SAMPLE CORRELATION, not the raw product-moment mean(z1*z2).  The latent
// pair is bivariate normal with unit variances, so the constrained MLE is the
// product-moment -- but ONLY when the draws actually have unit variance, and
// they do not: an unmixed chain gives a pooled spread of 2.15, whence
// mean(z^2) ~ 4.6 and a true correlation of 0.3 comes out as 0.3*4.6 = 1.38,
// which clamps to 0.999.  It then STAYS there, because 0.999*mean(z^2) is
// itself ~0.999 -- and a pinned correlation makes the copula partner's latent
// numerically equal to its partner's (w_k = rho*z_j + sqrt(1-rho^2)*z_k with
// rho ~ 1), so BOTH declared families end up fitted to the same draws.
//
// Measured on Bauer's gamma model: the correlation was pinned at 0.999 by the
// first M-step of every diverging run and never moved again.
//
// Normalizing by the observed spreads costs nothing, is bounded in [-1, 1] by
// construction, and agrees with the product-moment exactly when the draws do
// have unit variance -- so it is strictly the safer estimator here.
double rxEtaDistCorMleW(const std::vector<double> &z1,
                        const std::vector<double> &z2,
                        const std::vector<double> *w) {
  size_t n = std::min(z1.size(), z2.size());
  if (n < 2) return NA_REAL;
  if (w != nullptr && w->size() < n) return NA_REAL;
  double sw = 0.0, m1 = 0.0, m2 = 0.0; size_t m = 0;
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(z1[i]) || !std::isfinite(z2[i])) continue;
    double wi = (w == nullptr) ? 1.0 : (*w)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    m1 += wi*z1[i]; m2 += wi*z2[i]; sw += wi; m++;
  }
  if (m < 2 || !(sw > 0.0)) return NA_REAL;
  m1 /= sw; m2 /= sw;
  double s11 = 0.0, s22 = 0.0, s12 = 0.0;
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(z1[i]) || !std::isfinite(z2[i])) continue;
    double wi = (w == nullptr) ? 1.0 : (*w)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    double d1 = z1[i] - m1, d2 = z2[i] - m2;
    s11 += wi*d1*d1; s22 += wi*d2*d2; s12 += wi*d1*d2;
  }
  if (!(s11 > 0.0) || !(s22 > 0.0)) return NA_REAL;
  double r = s12 / std::sqrt(s11*s22);
  if (!std::isfinite(r)) return NA_REAL;
  return std::max(std::min(r, 0.999), -0.999);
}

double rxEtaDistCorMle(const std::vector<double> &z1,
                       const std::vector<double> &z2) {
  return rxEtaDistCorMleW(z1, z2, nullptr);
}
