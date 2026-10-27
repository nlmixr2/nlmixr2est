#ifndef NLMIXR2EST_MINQAC_H
#define NLMIXR2EST_MINQAC_H

// Derivative-free optimizers callable from C++ without a round-trip through R:
// RcppTrust's thread-safe NEWUOA/BOBYQA C API plus a 1-D Brent minimizer.
// The RcppTrust pointers are defined once by the iniRcppTrust expansion in
// inner.cpp and filled by .iniRcppTrustPtr() at load.

#include <RcppTrust.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <vector>

// Objective of the minqa-style C++ optimizers: value at an n-vector
// @noRd
typedef std::function<double(const double *)> minqaCFn;

// Result of newuoaC(), bobyqaC() and boundedResidOptC()
//
// `x` is the final point (the start when the solver is unavailable),
// `value` its objective (NA on failure), `convergence` RcppTrust's ierr
// (negative on failure) and `feval` the number of evaluations.
// @noRd
struct minqaCRes {
  std::vector<double> x;
  double value = NA_REAL;
  int convergence = -42;
  int feval = 0;
};

// C callback handed to RcppTrust; forwards to the wrapped minqaCFn
// @param n dimension (unused; the minqaCFn knows it)
// @param x point
// @param f objective value (output)
// @param ud the minqaCFn
// @return 0
// @noRd
static inline int minqaCTramp(int n, const double *x, double *f, void *ud) {
  (void)n;
  *f = (*static_cast<const minqaCFn *>(ud))(x);
  return 0;
}

// minqa's default rhobeg
// @param n dimension
// @param par start
// @return min(0.95, 0.2*max|par|), or 0.1 when par is all zero
// @noRd
static inline double minqaCRhobeg(int n, const double *par) {
  double mx = 0.0;
  for (int i = 0; i < n; ++i) mx = std::max(mx, std::fabs(par[i]));
  double rb = std::min(0.95, 0.2 * mx);
  return rb > 0.0 ? rb : 0.1;
}

// Copy an RcppTrust minqa result into a minqaCRes
// @param out result to fill
// @param n dimension
// @param res RcppTrust result (still owned by the caller)
// @param ierr solver return code
// @return nothing
// @noRd
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

// NEWUOA through RcppTrust's C API
// @param n dimension (>= 2)
// @param par start
// @param fn objective
// @param rhobeg initial trust radius; <= 0 selects minqa's default
// @param rhoend final trust radius; <= 0 selects 1e-6*rhobeg
// @param maxfun evaluation limit
// @param npt interpolation points; <= 0 selects 2n+1
// @return a minqaCRes
// @noRd
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

// BOBYQA through RcppTrust's C API
//
// Mirrors minqa::bobyqa()'s R-side setup: npt = 2n+1 as .bobyqa() sets it,
// and rhobeg shrunk to 0.2*min(upper-lower) for a tight box.
// @param n dimension (>= 2)
// @param par start, inside the box
// @param lower lower bounds
// @param upper upper bounds
// @param fn objective
// @param rhobeg initial trust radius; <= 0 selects minqa's default
// @param rhoend final trust radius; <= 0 selects 1e-6*rhobeg
// @param maxfun evaluation limit
// @return a minqaCRes
// @noRd
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

// State of the brentFmin() search
//
// Bracket [a, b], best point x, runner-up w, previous w in v (with their
// objectives), and the last two steps d and e.
// @noRd
struct brentState {
  double a, b, d, e, v, w, x, fv, fw, fx;
};

// Next brentFmin() trial point
//
// A parabolic step when it is acceptable, otherwise golden section.
// @param st search state (d and e are updated)
// @param xm bracket midpoint
// @param tol1 current tolerance
// @param t2 twice tol1
// @return the trial point u
// @noRd
static inline double brentStep(brentState &st, double xm, double tol1, double t2) {
  const double c = (3. - std::sqrt(5.)) * .5;
  double p = 0., q = 0., r = 0.;
  if (std::fabs(st.e) > tol1) {
    r = (st.x - st.w) * (st.fx - st.fv);
    q = (st.x - st.v) * (st.fx - st.fw);
    p = (st.x - st.v) * q - (st.x - st.w) * r;
    q = (q - r) * 2.;
    if (q > 0.) p = -p; else q = -q;
    r = st.e;
    st.e = st.d;
  }
  if (std::fabs(p) >= std::fabs(q * .5 * r) ||
      p <= q * (st.a - st.x) || p >= q * (st.b - st.x)) {
    st.e = (st.x < xm) ? st.b - st.x : st.a - st.x;
    st.d = c * st.e;
  } else {
    st.d = p / q;
    double u = st.x + st.d;
    if (u - st.a < t2 || st.b - u < t2) st.d = (st.x >= xm) ? -tol1 : tol1;
  }
  if (std::fabs(st.d) >= tol1) return st.x + st.d;
  return (st.d > 0.) ? st.x + tol1 : st.x - tol1;
}

// Update the brentFmin() bracket and best points after evaluating f(u)
// @param st search state
// @param u trial point
// @param fu objective at u
// @return nothing
// @noRd
static inline void brentUpdate(brentState &st, double u, double fu) {
  if (fu <= st.fx) {
    if (u < st.x) st.b = st.x; else st.a = st.x;
    st.v = st.w; st.w = st.x; st.x = u;
    st.fv = st.fw; st.fw = st.fx; st.fx = fu;
    return;
  }
  if (u < st.x) st.a = u; else st.b = u;
  if (fu <= st.fw || st.w == st.x) {
    st.v = st.w; st.fv = st.fw;
    st.w = u; st.fw = fu;
  } else if (fu <= st.fv || st.v == st.x || st.v == st.w) {
    st.v = u; st.fv = fu;
  }
}

// Port of R's Brent_fmin() (stats::optimize())
//
// Non-finite values become DBL_MAX, as optimize() does.
// @param ax lower end of the interval
// @param bx upper end of the interval
// @param f0 objective
// @param tol x tolerance
// @param nEval evaluation counter to increment, or NULL
// @return the minimizing x
// @noRd
static inline double brentFmin(double ax, double bx,
                               const std::function<double(double)> &f0,
                               double tol, int *nEval = NULL) {
  auto f = [&](double x) {
    if (nEval) (*nEval)++;
    double v = f0(x);
    return std::isfinite(v) ? v : DBL_MAX;
  };
  const double c = (3. - std::sqrt(5.)) * .5;
  // optimize() rejects tol <= 0; with tol 0 and the minimum at 0 the bracket never closes
  const double eps = std::sqrt(DBL_EPSILON), tol3 = std::max(tol, DBL_EPSILON) / 3.;
  brentState st;
  st.a = ax; st.b = bx;
  st.v = st.a + c * (st.b - st.a);
  st.w = st.v; st.x = st.v;
  st.d = 0.; st.e = 0.;
  st.fx = f(st.x);
  st.fv = st.fx; st.fw = st.fx;
  for (int it = 0; it < 1000; ++it) {
    double xm = (st.a + st.b) * .5;
    double tol1 = eps * std::fabs(st.x) + tol3;
    double t2 = tol1 * 2.;
    if (std::fabs(st.x - xm) <= t2 - (st.b - st.a) * .5) break;
    double u = brentStep(st, xm, tol1, t2);
    brentUpdate(st, u, f(u));
  }
  return st.x;
}

// Bracket a 1-D minimum by golden-ratio downhill expansion from x0
// @param f objective
// @param x0 start
// @param step first step
// @param lo lower end of the bracket (output)
// @param hi upper end of the bracket (output)
// @param nEval evaluation counter to increment, or NULL
// @param maxIt expansion limit
// @return false (leaving the widest bracket tried) if no minimum was bracketed
// @noRd
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

// Unbounded 1-D minimization: bracket from x0, then Brent on the bracket
// @param f objective
// @param x0 start
// @param tol x tolerance
// @param nEval evaluation counter to increment, or NULL
// @return the minimizing x
// @noRd
static inline double min1D(const std::function<double(double)> &f, double x0,
                           double tol, int *nEval = NULL) {
  double step = std::max(0.1 * std::fabs(x0), 0.1);
  double lo, hi;
  bracketMin1D(f, x0, step, lo, hi, nEval);
  return brentFmin(lo, hi, f, tol, nEval);
}

// C++ port of .boundedResidOpt()
//
// Infinite bounds are widened to a finite box around the start and the
// start is nudged inside it; bobyqa for n >= 2, bounded Brent for n == 1.
// @param n dimension
// @param par start
// @param lower lower bounds (may be infinite)
// @param upper upper bounds (may be infinite)
// @param fn objective
// @param rhoend bobyqa final trust radius; <= 0 selects the default
// @param maxfun bobyqa evaluation limit
// @return a minqaCRes
// @noRd
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
