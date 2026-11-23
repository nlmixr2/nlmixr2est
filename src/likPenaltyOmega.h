#ifndef __NLMIXR2EST_LIKPENALTY_OMEGA_H__
#define __NLMIXR2EST_LIKPENALTY_OMEGA_H__
// Penalized omega M-step shared by the EM-type methods (saem, imp family).
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <algorithm>

// S is the closed-form M-step update, the minimizer of the complete-data -2LL
// nsub*(log|G| + tr(G^-1 S)).  With an omega-dependent population penalty
// pen(G) (-2LL scale) S becomes the minimizer of that plus pen(G) over the
// lower-triangle entries (fr[a], fc[a]); the rest stay at S.  Damped Newton with
// an FD Hessian, G kept positive definite.  Returns false (S untouched) when
// the penalty does not depend on those entries.
template <class Pen>
static bool penOmegaMstepCore(arma::mat &S0, const std::vector<arma::uword> &fr,
                              const std::vector<arma::uword> &fc, double nsub, Pen pen) {
  const arma::uword nf = fr.size();
  if (nf == 0 || S0.n_elem == 0) return false;
  const arma::mat S = S0;
  auto build = [&](const arma::vec &x) {
    arma::mat G = S;
    for (arma::uword a = 0; a < nf; ++a) { G(fr[a], fc[a]) = x(a); G(fc[a], fr[a]) = x(a); }
    return G;
  };
  auto obj = [&](const arma::vec &x) {
    arma::mat G = build(x), L;
    if (!arma::chol(L, G)) return R_PosInf;
    arma::mat Gi = arma::inv_sympd(G);
    double v = 2.0 * nsub * arma::accu(arma::log(L.diag())) + nsub * arma::accu(Gi % S) + pen(G);
    return R_FINITE(v) ? v : R_PosInf;
  };
  arma::vec x(nf), hx(nf);
  for (arma::uword a = 0; a < nf; ++a) {
    x(a) = S(fr[a], fc[a]);
    hx(a) = 1e-5 * std::max(std::sqrt(std::fabs(S(fr[a], fr[a]) * S(fc[a], fc[a]))), 1e-8);
  }
  auto penGrad = [&](const arma::vec &xx) {
    arma::vec g(nf);
    for (arma::uword a = 0; a < nf; ++a) {
      arma::vec xp = xx, xm = xx;
      xp(a) += hx(a); xm(a) -= hx(a);
      g(a) = (pen(build(xp)) - pen(build(xm))) / (2.0 * hx(a));
    }
    return g;
  };
  if (arma::all(penGrad(x) == 0.0)) return false;
  auto grad = [&](const arma::vec &xx) {
    arma::mat Gi = arma::inv_sympd(build(xx));
    arma::mat D = nsub * (Gi - Gi * S * Gi);
    arma::vec g = penGrad(xx);
    for (arma::uword a = 0; a < nf; ++a) g(a) += (fr[a] == fc[a] ? 1.0 : 2.0) * D(fr[a], fc[a]);
    return g;
  };
  double fcur = obj(x);
  if (!R_FINITE(fcur)) return false;
  for (int it = 0; it < 50; ++it) {
    arma::vec g = grad(x);
    arma::mat H(nf, nf);
    bool hOk = true;
    for (arma::uword a = 0; a < nf; ++a) {
      double h = 10.0 * hx(a);
      arma::vec xp = x, xm = x;
      xp(a) += h; xm(a) -= h;
      if (!R_FINITE(obj(xp)) || !R_FINITE(obj(xm))) { hOk = false; break; }
      H.col(a) = (grad(xp) - grad(xm)) / (2.0 * h);
    }
    arma::vec step, ev;
    if (hOk) H = 0.5 * (H + H.t());
    if (!hOk || !arma::eig_sym(ev, H) || ev.min() <= 0.0 || !arma::solve(step, H, g) ||
        !step.is_finite()) {
      step = g * (std::max(1e-3, 0.1 * arma::norm(x)) / std::max(1e-300, arma::norm(g)));
    }
    double t = 1.0, fn = R_PosInf;
    arma::vec xn;
    for (int ls = 0; ls < 30; ++ls) {
      xn = x - t * step;
      fn = obj(xn);
      if (R_FINITE(fn) && fn <= fcur) break;
      t *= 0.5;
    }
    if (!(R_FINITE(fn) && fn <= fcur)) break;
    double dx = arma::norm(xn - x);
    x = xn; fcur = fn;
    if (dx < 1e-10 * (1.0 + arma::norm(x))) break;
  }
  S0 = build(x);
  return true;
}

#endif
