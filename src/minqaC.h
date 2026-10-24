#ifndef NLMIXR2EST_MINQAC_H
#define NLMIXR2EST_MINQAC_H

// Derivative-free optimizers callable from C++ without a round-trip through R:
// RcppTrust's thread-safe NEWUOA/BOBYQA C API plus a 1-D Brent minimizer.
// The RcppTrust pointers are defined once by the iniRcppTrust expansion in
// inner.cpp and filled by .iniRcppTrustPtr() at load.

#include <RcppTrust.h>
#include <cfloat>
#include <cmath>
#include <functional>
#include <vector>

typedef std::function<double(const double *)> minqaCFn;

struct minqaCRes {
  std::vector<double> x;
  double value = NA_REAL;
  int convergence = -42;
  int feval = 0;
};

static inline int minqaCTramp(int n, const double *x, double *f, void *ud) {
  (void)n;
  *f = (*static_cast<const minqaCFn *>(ud))(x);
  return 0;
}

// minqa's defaults for unset rhobeg/rhoend: rhobeg = min(0.95, 0.2*max|par|).
static inline double minqaCRhobeg(int n, const double *par) {
  double mx = 0.0;
  for (int i = 0; i < n; ++i) mx = std::max(mx, std::fabs(par[i]));
  double rb = std::min(0.95, 0.2 * mx);
  return rb > 0.0 ? rb : 0.1;
}

static inline void minqaCFill(minqaCRes &out, int n, minqa_result_t &res, int ierr) {
  out.convergence = ierr;
  out.feval = res.feval;
  if (ierr >= 0 && res.par != NULL) {
    out.x.assign(res.par, res.par + n);
    out.value = res.fval;
  } else {
    out.value = NA_REAL;
  }
}

// rhobeg/rhoend <= 0 select minqa's defaults; npt <= 0 selects 2n+1.
static inline minqaCRes newuoaC(int n, const double *par, const minqaCFn &fn,
                                double rhobeg, double rhoend, int maxfun,
                                int npt = 0) {
  minqaCRes out;
  out.x.assign(par, par + n);
  if (newuoa_solve_c_ptr == NULL || minqa_result_free_ptr == NULL) return out;
  minqa_options_t opts = minqa_options_default(n, par);
  opts.rhobeg = rhobeg > 0.0 ? rhobeg : minqaCRhobeg(n, par);
  opts.rhoend = rhoend > 0.0 ? rhoend : 1.0e-6 * opts.rhobeg;
  if (opts.rhoend > opts.rhobeg) opts.rhoend = opts.rhobeg;
  opts.npt = npt > 0 ? npt : 2 * n + 1;
  opts.maxfun = maxfun;
  minqa_result_t res;
  res.par = NULL; res.feval = 0; res.fval = NA_REAL;
  int ierr = newuoa_solve_c_ptr(n, par, minqaCTramp, (void *)&fn, &opts, &res);
  minqaCFill(out, n, res, ierr);
  minqa_result_free_ptr(&res);
  return out;
}

// Mirrors minqa::bobyqa()'s R-side setup (npt = 2n+1 as .bobyqa() sets it,
// rhobeg shrunk to 0.2*min(upper-lower) for a tight box).
static inline minqaCRes bobyqaC(int n, const double *par, const double *lower,
                                const double *upper, const minqaCFn &fn,
                                double rhobeg, double rhoend, int maxfun) {
  minqaCRes out;
  out.x.assign(par, par + n);
  if (bobyqa_solve_c_ptr == NULL || minqa_result_free_ptr == NULL) return out;
  minqa_options_t opts = minqa_options_default(n, par);
  opts.rhobeg = rhobeg > 0.0 ? rhobeg : minqaCRhobeg(n, par);
  opts.rhoend = rhoend > 0.0 ? rhoend : 1.0e-6 * opts.rhobeg;
  if (opts.rhoend > opts.rhobeg) opts.rhoend = opts.rhobeg;
  // like minqa, the shrink comes after the rhoend check, so rhoend > rhobeg can result
  double minRng = R_PosInf;
  for (int i = 0; i < n; ++i) minRng = std::min(minRng, upper[i] - lower[i]);
  if (minRng < 2.0 * opts.rhobeg) opts.rhobeg = 0.2 * minRng;
  opts.npt = 2 * n + 1;
  opts.maxfun = maxfun;
  minqa_result_t res;
  res.par = NULL; res.feval = 0; res.fval = NA_REAL;
  int ierr = bobyqa_solve_c_ptr(n, par, lower, upper, minqaCTramp, (void *)&fn,
                                &opts, &res);
  minqaCFill(out, n, res, ierr);
  minqa_result_free_ptr(&res);
  return out;
}

// Port of R's Brent_fmin() (stats::optimize); non-finite values become DBL_MAX
// as optimize() does.
static inline double brentFmin(double ax, double bx,
                               const std::function<double(double)> &f0,
                               double tol, int *nEval = NULL) {
  auto f = [&](double x) {
    if (nEval) (*nEval)++;
    double v = f0(x);
    return std::isfinite(v) ? v : DBL_MAX;
  };
  const double c = (3. - std::sqrt(5.)) * .5;
  double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
  eps = std::sqrt(DBL_EPSILON);
  a = ax; b = bx;
  v = a + c * (b - a);
  w = v; x = v;
  d = 0.; e = 0.;
  fx = f(x);
  fv = fx; fw = fx;
  tol3 = tol / 3.;
  for (;;) {
    xm = (a + b) * .5;
    tol1 = eps * std::fabs(x) + tol3;
    t2 = tol1 * 2.;
    if (std::fabs(x - xm) <= t2 - (b - a) * .5) break;
    p = 0.; q = 0.; r = 0.;
    if (std::fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if (q > 0.) p = -p; else q = -q;
      r = e;
      e = d;
    }
    if (std::fabs(p) >= std::fabs(q * .5 * r) ||
        p <= q * (a - x) || p >= q * (b - x)) {
      if (x < xm) e = b - x; else e = a - x;
      d = c * e;
    } else {
      d = p / q;
      u = x + d;
      if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm) d = -d;
      }
    }
    if (std::fabs(d) >= tol1) u = x + d;
    else if (d > 0.) u = x + tol1;
    else u = x - tol1;
    fu = f(u);
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; w = x; x = u;
      fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  return x;
}

// Golden-ratio downhill expansion from x0 until the minimum is bracketed;
// returns false (leaving the widest bracket tried) if it never turns up.
static inline bool bracketMin1D(const std::function<double(double)> &f,
                                double x0, double step, double &lo, double &hi,
                                int *nEval = NULL, int maxIt = 60) {
  const double gold = 1.618034;
  double a = x0, b = x0 + step;
  double fa = f(a), fb = f(b);
  if (nEval) *nEval += 2;
  if (!std::isfinite(fa)) fa = DBL_MAX;
  if (!std::isfinite(fb)) fb = DBL_MAX;
  if (fb > fa) { std::swap(a, b); std::swap(fa, fb); }
  double c = b + gold * (b - a), fc = f(c);
  if (nEval) (*nEval)++;
  if (!std::isfinite(fc)) fc = DBL_MAX;
  int it = 0;
  while (fb > fc && it++ < maxIt) {
    a = b; fa = fb;
    b = c; fb = fc;
    c = b + gold * (b - a);
    fc = f(c);
    if (nEval) (*nEval)++;
    if (!std::isfinite(fc)) fc = DBL_MAX;
  }
  lo = std::min(a, c);
  hi = std::max(a, c);
  return fb <= fc;
}

// Unbounded 1-D minimization: bracket from x0, then Brent on the bracket.
static inline double min1D(const std::function<double(double)> &f, double x0,
                           double tol, int *nEval = NULL) {
  double step = std::max(0.1 * std::fabs(x0), 0.1);
  double lo, hi;
  bracketMin1D(f, x0, step, lo, hi, nEval);
  return brentFmin(lo, hi, f, tol, nEval);
}

// Port of .boundedResidOpt(): infinite bounds widened to a finite box around
// the start; bobyqa for n >= 2, bounded Brent for n == 1.
static inline minqaCRes boundedResidOptC(int n, const double *par,
                                         const double *lower, const double *upper,
                                         const minqaCFn &fn, double rhoend = 0.0,
                                         int maxfun = 10000) {
  std::vector<double> lo(n), hi(n), p(n);
  for (int i = 0; i < n; ++i) {
    double span = std::max(std::fabs(par[i]), 1.0) * 1e3;
    lo[i] = std::isfinite(lower[i]) ? lower[i] : par[i] - span;
    hi[i] = std::isfinite(upper[i]) ? upper[i] : par[i] + span;
    double w = hi[i] - lo[i];
    p[i] = std::min(std::max(par[i], lo[i] + 1e-8 * w), hi[i] - 1e-8 * w);
  }
  if (n == 1) {
    minqaCRes out;
    double x = brentFmin(lo[0], hi[0],
                         [&](double v) { return fn(&v); },
                         std::pow(DBL_EPSILON, 0.25));
    out.x.assign(1, x);
    out.value = fn(&x);
    out.convergence = 0;
    return out;
  }
  return bobyqaC(n, p.data(), lo.data(), hi.data(), fn, 0.0, rhoend, maxfun);
}

#endif // NLMIXR2EST_MINQAC_H
