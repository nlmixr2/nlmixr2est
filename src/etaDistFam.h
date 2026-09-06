#ifndef __ETADISTFAM_H__
#define __ETADISTFAM_H__
// Declared random-effect distributions: family dispatch and the ODE-free
// maximum-likelihood M-step, shared by saem (src/saem.cpp) and imp/impmap
// (src/imp.cpp).  Kept in its own translation unit so the two estimators
// cannot drift apart in which families they support or how they fit them.
//
// Everything here is pure Rmath -- no ODE solve, no rxode2 state -- so it is
// callable from either estimator's M-step and from anywhere in the outer
// problem that needs a distribution fit to a set of etas.
#include <cmath>
#include <cstddef>
#include <vector>
// R::qnorm(), R::dgamma() and friends are Rcpp's namespaced Rmath wrappers, not
// the bare C API, so this needs Rcpp rather than <Rmath.h>.
#include <Rcpp.h>

// Family codes are the ROW NUMBER in lotri::lotriEtaDists(), so this dispatch
// and the catalog cannot drift apart.  .etaDistFamilyCode() (R/etaDistMstep.R)
// assigns them; a family this dispatch does not implement gets code 0 and falls
// back to the general R path, which evaluates the declaration's own d*() call.
// Arguments are the family's NATIVE parameters in the order the d*()/q*() pair
// takes them, EXCEPT that a rate is converted to Rmath's scale here rather
// than in R.
#define RXETADIST_NORM        1   // mean, sd
#define RXETADIST_STDNORMAL   2   // (none)
#define RXETADIST_STUDENTT    3   // nu, mu, sigma
#define RXETADIST_CAUCHY      4   // location, scale
#define RXETADIST_DBLEXP      5   // mu, sigma
#define RXETADIST_LOGIS       6   // location, scale
#define RXETADIST_GUMBEL      7   // mu, beta
#define RXETADIST_LNORM       8   // meanlog, sdlog
#define RXETADIST_CHISQ       9   // df
#define RXETADIST_INVCHISQ   10   // nu
#define RXETADIST_SCINVCHISQ 11   // nu, sigma
#define RXETADIST_EXP        12   // rate
#define RXETADIST_GAMMA      13   // shape, rate
#define RXETADIST_INVGAMMA   14   // alpha, beta
#define RXETADIST_WEIBULL    15   // shape, scale
#define RXETADIST_FRECHET    16   // alpha, sigma
#define RXETADIST_RAYLEIGH   17   // sigma
#define RXETADIST_PARETO     18   // y_min, alpha
#define RXETADIST_PARETO2    19   // mu, lambda, alpha
#define RXETADIST_BETA       20   // shape1, shape2
#define RXETADIST_BETAPROP   21   // mu, kappa
#define RXETADIST_UNIF       22   // min, max

static inline int rxEtaDistNarg(int fam) {
  switch (fam) {
  case RXETADIST_STDNORMAL:                      return 0;
  case RXETADIST_CHISQ:  case RXETADIST_INVCHISQ:
  case RXETADIST_EXP:    case RXETADIST_RAYLEIGH: return 1;
  case RXETADIST_STUDENTT: case RXETADIST_PARETO2: return 3;
  case RXETADIST_NORM:   case RXETADIST_CAUCHY: case RXETADIST_DBLEXP:
  case RXETADIST_LOGIS:  case RXETADIST_GUMBEL: case RXETADIST_LNORM:
  case RXETADIST_SCINVCHISQ: case RXETADIST_GAMMA: case RXETADIST_INVGAMMA:
  case RXETADIST_WEIBULL: case RXETADIST_FRECHET: case RXETADIST_PARETO:
  case RXETADIST_BETA:   case RXETADIST_BETAPROP: case RXETADIST_UNIF: return 2;
  default: return -1;                 // unimplemented -> R fallback
  }
}

// Native parameters constrained positive, as a bit mask.  nelder_fn is
// unbounded, so the objective optimizes log() of these.
static inline int rxEtaDistPosMask(int fam) {
  switch (fam) {
  case RXETADIST_NORM: case RXETADIST_CAUCHY: case RXETADIST_DBLEXP:
  case RXETADIST_LOGIS: case RXETADIST_GUMBEL: case RXETADIST_LNORM:
  case RXETADIST_BETAPROP:                                    return 0x2;
  case RXETADIST_STUDENTT:                                    return 0x5; // nu, sigma
  case RXETADIST_PARETO2:                                     return 0x6; // lambda, alpha
  case RXETADIST_CHISQ: case RXETADIST_INVCHISQ:
  case RXETADIST_EXP:   case RXETADIST_RAYLEIGH:              return 0x1;
  case RXETADIST_SCINVCHISQ: case RXETADIST_GAMMA:
  case RXETADIST_INVGAMMA: case RXETADIST_WEIBULL:
  case RXETADIST_FRECHET: case RXETADIST_PARETO:
  case RXETADIST_BETA:                                        return 0x3;
  default:                                                    return 0x0;
  }
}

// quantile: latent uniform -> eta.  Mirrors the catalog's own templates.
static inline double rxEtaDistQ(int fam, double u, const double *a) {
  switch (fam) {
  case RXETADIST_NORM:      return R::qnorm(u, a[0], a[1], 1, 0);
  case RXETADIST_STDNORMAL: return R::qnorm(u, 0.0, 1.0, 1, 0);
  case RXETADIST_STUDENTT:  return a[1] + a[2]*R::qt(u, a[0], 1, 0);
  case RXETADIST_CAUCHY:    return R::qcauchy(u, a[0], a[1], 1, 0);
  case RXETADIST_DBLEXP: {
    double sg = (u < 0.5) ? -1.0 : 1.0;
    return a[0] - a[1]*sg*std::log1p(-2.0*sg*(u - 0.5));
  }
  case RXETADIST_LOGIS:     return R::qlogis(u, a[0], a[1], 1, 0);
  case RXETADIST_GUMBEL:    return a[0] - a[1]*std::log(-std::log(u));
  case RXETADIST_LNORM:     return R::qlnorm(u, a[0], a[1], 1, 0);
  case RXETADIST_CHISQ:     return R::qchisq(u, a[0], 1, 0);
  // 2*gammapInv(nu/2, 1-u) IS qchisq(1-u, nu)
  case RXETADIST_INVCHISQ:  return 1.0/R::qchisq(1.0 - u, a[0], 1, 0);
  case RXETADIST_SCINVCHISQ: return a[0]*a[1]*a[1]/R::qchisq(1.0 - u, a[0], 1, 0);
  case RXETADIST_EXP:       return R::qexp(u, 1.0/a[0], 1, 0);
  case RXETADIST_GAMMA:     return R::qgamma(u, a[0], 1.0/a[1], 1, 0);
  case RXETADIST_INVGAMMA:  return a[1]/R::qgamma(1.0 - u, a[0], 1.0, 1, 0);
  case RXETADIST_WEIBULL:   return R::qweibull(u, a[0], a[1], 1, 0);
  case RXETADIST_FRECHET:   return a[1]*std::pow(-std::log(u), -1.0/a[0]);
  case RXETADIST_RAYLEIGH:  return a[0]*std::sqrt(-2.0*std::log1p(-u));
  case RXETADIST_PARETO:    return a[0]*std::pow(1.0 - u, -1.0/a[1]);
  case RXETADIST_PARETO2:   return a[0] + a[1]*(std::pow(1.0 - u, -1.0/a[2]) - 1.0);
  case RXETADIST_BETA:      return R::qbeta(u, a[0], a[1], 1, 0);
  case RXETADIST_BETAPROP:  return R::qbeta(u, a[0]*a[1], (1.0 - a[0])*a[1], 1, 0);
  case RXETADIST_UNIF:      return R::qunif(u, a[0], a[1], 1, 0);
  default:                  return NA_REAL;
  }
}

// log density at an eta value
static inline double rxEtaDistLogD(int fam, double x, const double *a) {
  switch (fam) {
  case RXETADIST_NORM:      return R::dnorm(x, a[0], a[1], 1);
  case RXETADIST_STDNORMAL: return R::dnorm(x, 0.0, 1.0, 1);
  case RXETADIST_STUDENTT:  // location-scale t: dt(z)/sigma
    return R::dt((x - a[1])/a[2], a[0], 1) - std::log(a[2]);
  case RXETADIST_CAUCHY:    return R::dcauchy(x, a[0], a[1], 1);
  case RXETADIST_DBLEXP:
    return -std::log(2.0*a[1]) - std::fabs(x - a[0])/a[1];
  case RXETADIST_LOGIS:     return R::dlogis(x, a[0], a[1], 1);
  case RXETADIST_GUMBEL: {
    double z = (x - a[0])/a[1];
    return -std::log(a[1]) - z - std::exp(-z);
  }
  case RXETADIST_LNORM:     return R::dlnorm(x, a[0], a[1], 1);
  case RXETADIST_CHISQ:     return R::dchisq(x, a[0], 1);
  case RXETADIST_INVCHISQ:  // X = 1/Y, Y~chisq(nu); |dY/dX| = 1/x^2
    return (x > 0) ? R::dchisq(1.0/x, a[0], 1) - 2.0*std::log(x) : R_NegInf;
  case RXETADIST_SCINVCHISQ: {
    if (x <= 0) return R_NegInf;
    double nu = a[0], t2 = a[1]*a[1];
    return (nu/2.0)*std::log(nu*t2/2.0) - R::lgammafn(nu/2.0)
      - (1.0 + nu/2.0)*std::log(x) - nu*t2/(2.0*x);
  }
  case RXETADIST_EXP:       return R::dexp(x, 1.0/a[0], 1);
  case RXETADIST_GAMMA:     return R::dgamma(x, a[0], 1.0/a[1], 1);
  case RXETADIST_INVGAMMA:
    return (x > 0) ? a[0]*std::log(a[1]) - R::lgammafn(a[0])
      - (a[0] + 1.0)*std::log(x) - a[1]/x : R_NegInf;
  case RXETADIST_WEIBULL:   return R::dweibull(x, a[0], a[1], 1);
  case RXETADIST_FRECHET: {
    if (x <= 0) return R_NegInf;
    double z = x/a[1];
    return std::log(a[0]/a[1]) - (1.0 + a[0])*std::log(z) - std::pow(z, -a[0]);
  }
  case RXETADIST_RAYLEIGH:
    return (x > 0) ? std::log(x) - 2.0*std::log(a[0]) - x*x/(2.0*a[0]*a[0])
      : R_NegInf;
  case RXETADIST_PARETO:
    return (x >= a[0]) ? std::log(a[1]) + a[1]*std::log(a[0])
      - (a[1] + 1.0)*std::log(x) : R_NegInf;
  case RXETADIST_PARETO2: {
    double z = (x - a[0])/a[1];
    return (z >= 0) ? std::log(a[2]/a[1]) - (a[2] + 1.0)*std::log1p(z) : R_NegInf;
  }
  case RXETADIST_BETA:      return R::dbeta(x, a[0], a[1], 1);
  case RXETADIST_BETAPROP:
    return R::dbeta(x, a[0]*a[1], (1.0 - a[0])*a[1], 1);
  case RXETADIST_UNIF:      return R::dunif(x, a[0], a[1], 1);
  default:                  return R_NegInf;
  }
}

// Fit `fam` to `vals` by maximum likelihood, starting from the native
// parameters in `a0` (length rxEtaDistNarg(fam)).  On success `a0` holds the
// fitted parameters; on failure it is left untouched and false is returned.
//
// `w`, when non-null, is a vector of nonnegative weights of the same length as
// `vals` -- the importance weights imp carries on each sampled eta.  A null `w`
// means unit weights, which is what saem's equally-weighted MCMC draws want.
bool rxEtaDistMleW(int fam, const std::vector<double> &vals,
                   const std::vector<double> *w, double *a0);
bool rxEtaDistMle(int fam, const std::vector<double> &vals, double *a0);

// Is this pooled set of latent draws worth fitting?
//
// The latent normals are standard normal BY CONSTRUCTION -- rxEtaDistExpand()
// gives them a unit, fixed omega -- so the only reason a pooled sample departs
// from that is information about the declared family, or a sampler that has not
// settled.  The M-step cannot tell those apart and acts on the difference, so
// acting on the second is a runaway: measured on Bauer's gamma model, an
// unmixed first iteration (pooled sd 2.15) drove the fitted shape 7.39 -> 1.74
// -> 0.51, each shrink widening the mapped etas and feeding the next.  Once
// settled the same fit sits at sd ~0.87.
//
// Returns false when the spread is outside [lo, hi]; `sdOut`, when non-null,
// receives the measured spread either way (for reporting).
// `wt`, when non-null, weights the draws -- imp's importance weights.  This
// matters there: imp draws from N(MAP, gamma*H^-1) with gamma >= 1, so the
// RAW sample spread is over-dispersed by design and would trip the upper bound
// on a perfectly healthy fit, silently disabling the M-step.  The weighted
// spread is the posterior's, which is what the bound is about.
bool rxEtaDistSpreadOk(const std::vector<double> &w, double lo, double hi,
                       double *sdOut, const std::vector<double> *wt = nullptr);

// Closed-form M-step for a Gaussian copula's correlation from paired latent
// draws.  `w` weights as above.
double rxEtaDistCorMleW(const std::vector<double> &w1,
                        const std::vector<double> &w2,
                        const std::vector<double> *w);
double rxEtaDistCorMle(const std::vector<double> &w1,
                       const std::vector<double> &w2);

#endif // __ETADISTFAM_H__
