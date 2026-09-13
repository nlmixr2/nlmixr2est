#ifndef __ETADIST_ETASCALE_H__
#define __ETADIST_ETASCALE_H__
// The declared eta on ITS OWN scale: where it lives, how to map that to the
// whole real line and back, and how its log density moves with it.
//
// Everything else in etaDistFam.h works on the family's ARGUMENTS -- fitting
// them, differentiating with respect to them, inverting them.  This file is the
// other axis: the eta is the variable and the arguments are held.  It is what a
// sampler proposing the eta directly needs (a bijected random walk carries a
// log-Jacobian into the acceptance ratio) and what an inner MAP needs to use a
// non-Gaussian prior (the quadratic eta'Omega^-1 eta and its curvature are the
// GAUSSIAN SPECIAL CASE of -2 log p and its second derivative).
//
// The support is NOT lotri's `support` column.  That column says what the
// family is morally -- and for a bijector it is wrong in three places:
//
//   dunif        listed "real"      really [min, max], its own two arguments
//   pareto       listed "positive"  really [y_min, inf), y_min = a[0]
//   paretoType2  listed "nonneg"    really [mu, inf),   mu    = a[0]
//
// Taking the column at face value would put an identity bijector on a bounded
// eta, and every proposal outside the range would come back -Inf with nothing
// saying why.  So bounds are computed per family FROM THE ARGUMENTS here.

#include "etaDistFam.h"
#include <cmath>

// Support of `fam` at arguments `a`, as the interval the eta may take.
static inline void rxEtaDistBounds(int fam, const double *a,
                                   double *lo, double *hi) {
  *lo = R_NegInf; *hi = R_PosInf;
  switch (fam) {
  case RXETADIST_NORM: case RXETADIST_STDNORMAL: case RXETADIST_STUDENTT:
  case RXETADIST_CAUCHY: case RXETADIST_DBLEXP: case RXETADIST_LOGIS:
  case RXETADIST_GUMBEL:
    return;                                   // the whole line
  case RXETADIST_LNORM: case RXETADIST_CHISQ: case RXETADIST_INVCHISQ:
  case RXETADIST_SCINVCHISQ: case RXETADIST_EXP: case RXETADIST_GAMMA:
  case RXETADIST_INVGAMMA: case RXETADIST_WEIBULL: case RXETADIST_FRECHET:
  case RXETADIST_RAYLEIGH:
    *lo = 0.0; return;
  case RXETADIST_PARETO:  *lo = a[0]; return;  // y_min, not 0
  case RXETADIST_PARETO2: *lo = a[0]; return;  // mu,    not 0
  case RXETADIST_BETA: case RXETADIST_BETAPROP:
    *lo = 0.0; *hi = 1.0; return;
  case RXETADIST_UNIF:    *lo = a[0]; *hi = a[1]; return;  // NOT "real"
  default: return;
  }
}

// The four bijector shapes, chosen by which bounds are finite.  Same vocabulary
// as R's .saemPseudoEtaLine(): identity, exp, reflected exp, expit.
#define RXETADIST_BIJ_ID    0
#define RXETADIST_BIJ_LO    1
#define RXETADIST_BIJ_HI    2
#define RXETADIST_BIJ_BOTH  3

static inline int rxEtaDistBijKind(double lo, double hi) {
  bool fLo = R_finite(lo), fHi = R_finite(hi);
  if (fLo && fHi) return RXETADIST_BIJ_BOTH;
  if (fLo) return RXETADIST_BIJ_LO;
  if (fHi) return RXETADIST_BIJ_HI;
  return RXETADIST_BIJ_ID;
}

// eta -> u on the whole real line.
static inline double rxEtaDistToU(int fam, double x, const double *a) {
  double lo, hi; rxEtaDistBounds(fam, a, &lo, &hi);
  switch (rxEtaDistBijKind(lo, hi)) {
  case RXETADIST_BIJ_LO:   return std::log(x - lo);
  case RXETADIST_BIJ_HI:   return std::log(hi - x);
  case RXETADIST_BIJ_BOTH: { double p = (x - lo)/(hi - lo); return std::log(p/(1.0 - p)); }
  default:                 return x;
  }
}

// u -> eta.  The inverse of rxEtaDistToU, and the one a proposal goes through.
static inline double rxEtaDistFromU(int fam, double u, const double *a) {
  double lo, hi; rxEtaDistBounds(fam, a, &lo, &hi);
  switch (rxEtaDistBijKind(lo, hi)) {
  case RXETADIST_BIJ_LO:   return lo + std::exp(u);
  case RXETADIST_BIJ_HI:   return hi - std::exp(u);
  case RXETADIST_BIJ_BOTH: { double e = 1.0/(1.0 + std::exp(-u)); return lo + (hi - lo)*e; }
  default:                 return u;
  }
}

// log |d(eta)/d(u)| at `u` -- the term a bijected random walk must add to the
// acceptance ratio, and the reason a naive random walk on a positive eta is
// wrong rather than merely inefficient.
static inline double rxEtaDistLogJac(int fam, double u, const double *a) {
  double lo, hi; rxEtaDistBounds(fam, a, &lo, &hi);
  switch (rxEtaDistBijKind(lo, hi)) {
  case RXETADIST_BIJ_LO: case RXETADIST_BIJ_HI:
    return u;                                     // d/du (lo + exp u) = exp u
  case RXETADIST_BIJ_BOTH: {
    // log((hi-lo) * e * (1-e)) written through log1p so it survives |u| large
    double s = (u > 0) ? -u : u;                  // = -|u|
    return std::log(hi - lo) + s - 2.0*std::log1p(std::exp(s));
  }
  default:
    return 0.0;
  }
}

// A step for differencing `x` that stays inside the support.  Relative to |x|
// so it works at 1e-3 and at 1e3, and shrunk near a bound rather than clamped
// onto it -- a difference taken AT the bound of a density that diverges there
// (gamma with shape < 1) is the case this exists for.
static inline double rxEtaDistDiffStep(double x, double lo, double hi) {
  double h = 1e-5 * ((std::fabs(x) > 1.0) ? std::fabs(x) : 1.0);
  if (R_finite(lo) && x - h <= lo) h = 0.25*(x - lo);
  if (R_finite(hi) && x + h >= hi) h = 0.25*(hi - x);
  return (h > 0.0 && R_finite(h)) ? h : 0.0;
}

// d log p / d(eta), central differenced.  The log density is pure arithmetic
// with no solve behind it, so this costs a few nanoseconds -- the same trade
// rxEtaDistLoglikGrad() already makes for the argument derivatives.
static inline double rxEtaDistLogDdx(int fam, double x, const double *a) {
  double lo, hi; rxEtaDistBounds(fam, a, &lo, &hi);
  double h = rxEtaDistDiffStep(x, lo, hi);
  if (h == 0.0) return R_NaN;
  double f1 = rxEtaDistLogD(fam, x + h, a), f0 = rxEtaDistLogD(fam, x - h, a);
  if (!R_finite(f1) || !R_finite(f0)) return R_NaN;
  return (f1 - f0)/(2.0*h);
}

// d2 log p / d(eta)2.  NEGATED this is the prior curvature an inner MAP adds to
// the Hessian, of which Omega^-1 is the Gaussian case.
static inline double rxEtaDistLogD2dx(int fam, double x, const double *a) {
  double lo, hi; rxEtaDistBounds(fam, a, &lo, &hi);
  double h = rxEtaDistDiffStep(x, lo, hi);
  if (h == 0.0) return R_NaN;
  double f1 = rxEtaDistLogD(fam, x + h, a), fc = rxEtaDistLogD(fam, x, a),
    f0 = rxEtaDistLogD(fam, x - h, a);
  if (!R_finite(f1) || !R_finite(fc) || !R_finite(f0)) return R_NaN;
  return (f1 - 2.0*fc + f0)/(h*h);
}

// ---------------------------------------------------------------------------
// The JOINT density of a declared pair, on the eta scale, under a Gaussian
// copula.
//
// This is what lets the direct route carry a CORRELATED declared block, which
// it otherwise has to refuse: a Gaussian copula over non-normal marginals is
// exactly `eta = Q(phi(z))`, so refusing the block was refusing the only thing
// that gives the correlation meaning.  Written on the eta scale it is just
// another prior term:
//
//   log p(eta1, eta2) = log f1(eta1) + log f2(eta2) + log c_rho(u1, u2)
//   u_i = F_i(eta_i),  z_i = qnorm(u_i)
//   log c_rho = -0.5 log(1-rho^2)
//               - (rho^2 (z1^2 + z2^2) - 2 rho z1 z2) / (2 (1-rho^2))
//
// The point of writing it here rather than decoding in the model: the inverse
// CDF then runs once per ETA PER SUBJECT instead of once per OBSERVATION, and
// the correlation is estimated from the etas' own u values rather than from
// "combined latents", which is a biased estimator with a fixed point at the
// current rho.

// log of the bivariate Gaussian copula density at the two marginal z values.
static inline double rxEtaDistCopulaLogC(double z1, double z2, double rho) {
  if (!R_finite(z1) || !R_finite(z2)) return R_NegInf;
  double r2 = rho*rho;
  double om = 1.0 - r2;
  if (om <= 0.0) return R_NegInf;          // a degenerate copula has no density
  return -0.5*std::log(om) - (r2*(z1*z1 + z2*z2) - 2.0*rho*z1*z2)/(2.0*om);
}

// The marginal z the copula uses: z = qnorm(F(eta)).  Clamped off 0 and 1 the
// same way the decoder clamps phiU(), because an eta in the far tail otherwise
// saturates the CDF in double precision and returns +/-Inf for a value the
// model is perfectly happy with.
static inline double rxEtaDistCopulaZ(int fam, double x, const double *a) {
  double u = rxEtaDistP(fam, x, a);
  if (!R_finite(u)) return R_NaN;
  if (u < 1e-15) u = 1e-15; else if (u > 1.0 - 1e-15) u = 1.0 - 1e-15;
  return R::qnorm(u, 0.0, 1.0, 1, 0);
}

// Joint log density of a declared PAIR on the eta scale.  `fam1`/`a1` and
// `fam2`/`a2` are the two marginals; `rho` the copula correlation.
static inline double rxEtaDistPairLogD(int fam1, double x1, const double *a1,
                                       int fam2, double x2, const double *a2,
                                       double rho) {
  double l1 = rxEtaDistLogD(fam1, x1, a1);
  double l2 = rxEtaDistLogD(fam2, x2, a2);
  if (!R_finite(l1) || !R_finite(l2)) return R_NegInf;
  if (rho == 0.0) return l1 + l2;          // independent: no copula term at all
  double z1 = rxEtaDistCopulaZ(fam1, x1, a1);
  double z2 = rxEtaDistCopulaZ(fam2, x2, a2);
  double lc = rxEtaDistCopulaLogC(z1, z2, rho);
  if (!R_finite(lc)) return R_NegInf;
  return l1 + l2 + lc;
}

#endif // __ETADIST_ETASCALE_H__
