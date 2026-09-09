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
#include <string>
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
// Map a family's fitted NATIVE parameters back onto the user's thetas, in C++.
// Returns false when an argument expression is outside the C++ grammar or the
// solve does not converge -- the caller then keeps the R route.
#include "etaDistExpr.h"

// ---------------------------------------------------------------------------
// The general declared-distribution M-step objective (see the design block in
// R/etaDistMstep.R).  C++ because it has to run inside the OpenMP regions the
// focei family and imp use -- the R API cannot be touched from a parallel
// region at all, so an R implementation could never be the target.
//
//   sum over (subject i, record j with evid == 0) of
//      w_ij * log p_family( eta_ij ; args_ij(theta) )
//
// Everything it needs is already thread-safe: etaDistExprParse/Eval are
// header-only pure arithmetic, and rxEtaDistGradD dispatches through
// rxode2ll's exported function pointers.
//
// The argument expressions are parsed ONCE against a name list holding the
// thetas FIRST and then the per-record symbols, so evaluation is a flat array
// lookup: vals[0..nth) are the candidate thetas and vals[nth..) that record's
// own values.  A covariate needs no special case -- it is simply a name whose
// value differs by record.
//
// `rec` is nRec x nSym in row-major order, `etaAt` is the random effect held
// at the CURRENT parameters (design Q1), and `wt` carries the per-record
// weight (design Q2: 1/n_i, so each SUBJECT contributes one unit however many
// times it was observed).
//
// Returns false when any record cannot be evaluated -- a partial sum would
// silently drop whichever subjects failed.
bool rxEtaDistLoglikObj(int fam,
                        const std::vector< std::vector<etaDistTok> > &rpn,
                        int nth, int nSym,
                        const double *theta,
                        const double *rec, const double *etaAt,
                        const double *wt, int nRec,
                        double *out);

// Objective AND its gradient with respect to the thetas.
//
// The chain rule is complete here because both halves already exist:
//
//   dL/dtheta_t = sum_r w_r * sum_k  g_k(r) * da_k/dtheta_t (r)
//
//   g_k        = d(log p)/d(a_k), returned by rxEtaDistGradD() from rxode2ll's
//                exact Stan-backed derivatives -- this is what the function
//                pointer export was FOR.
//   da_k/dth_t = derivative of the argument expression.  Central-differenced on
//                the RPN, which is pure arithmetic with no solve behind it, so
//                it costs a handful of nanoseconds per record and is accurate
//                to ~1e-8 -- ample for a quasi-Newton step.
//
// Same contract as rxEtaDistLoglikObj(): false when any record cannot be
// evaluated, since a partial sum silently drops subjects.
bool rxEtaDistLoglikGrad(int fam,
                         const std::vector< std::vector<etaDistTok> > &rpn,
                         int nth, int nSym,
                         const double *theta,
                         const double *rec, const double *etaAt,
                         const double *wt, int nRec,
                         double *out, double *grad);

// Parse a declaration's argument expressions for rxEtaDistLoglikObj().
// `names` must be the thetas followed by the per-record symbols, in the same
// order the `theta`/`rec` arrays supply them.  Returns false when any
// expression falls outside the evaluator's grammar, which is reported rather
// than guessed at.
bool rxEtaDistLoglikParse(const std::vector<std::string> &exprs,
                          const std::vector<std::string> &names,
                          std::vector< std::vector<etaDistTok> > &rpn);

bool rxEtaDistArgsToThetas(const std::vector<std::string> &exprs,
                           const std::vector<std::string> &thetaNames,
                           const double *start, const double *target,
                           double *out);

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
// Has family k's pooled latent spread STOPPED CHANGING between M-step attempts?
//
// A LEVEL test cannot do this job, and the [0.5, 1.0] band that used to be
// written at each call site is measurably wrong.  Under a wrong family a fully
// mixed chain sits far from 1 -- 1.40 on Bauer's g1 -- and that spread IS the
// information the M-step consumes; a still-burning chain passes through the
// same 1.40 on its way down from 2.6, where acting on it diverges.  The two are
// indistinguishable by value and obvious by trajectory.  Measured on g1, saem
// at 300 burn + 150 EM: the level band scored 19.7% against 9.2% for a band
// wide enough to admit the settled spread -- but that same wide band scored
// 255.1% at 60 + 30, where the chain has not settled.  Across the four Bauer
// arms the [0.5, 1.0] band admitted exactly one (g3, settled spread 0.880).
//
// `lo`/`hi` remain as a loose DIVERGENCE cap only.  `tol <= 0` disables the
// settling test and leaves the cap alone.
//
// The caller measures `lsd` -- weighted for imp, unweighted for focei -- and
// owns `prev`/`cur`, so this stays free of estimator state.  `cur` must be
// cleared per attempt and copied into `prev` ONCE per attempt, after every loop
// that consulted it (see rxEtaDistSpreadAdvance): advancing inside a loop lets
// a later loop compare an attempt against itself, which always looks settled.
static inline bool rxEtaDistSpreadSettled(int k, double lsd,
                                          std::vector<double> &prev,
                                          std::vector<double> &cur,
                                          double lo, double hi, double tol) {
  if (k < 0) return false;
  if ((int)prev.size() <= k) prev.resize((size_t)k + 1, NA_REAL);
  if ((int)cur.size() <= k) cur.resize((size_t)k + 1, NA_REAL);
  cur[(size_t)k] = lsd;
  if (!std::isfinite(lsd)) {
    // a measurement that failed is not a gap to be spanned: drop the baseline
    // so the next attempt declines for want of one rather than silently
    // comparing across two gaps or more
    prev[(size_t)k] = NA_REAL;
    return false;
  }
  if (!(lsd >= lo && lsd <= hi)) return false;
  if (!(tol > 0.0)) return true;                  // cap only
  double p = prev[(size_t)k];
  if (!std::isfinite(p) || !(p > 0.0)) return false;
  return std::fabs(lsd - p) <= tol * p;
}

// ONE advance per M-step attempt, after every loop that consulted the baseline.
static inline void rxEtaDistSpreadAdvance(std::vector<double> &prev,
                                          std::vector<double> &cur) {
  if (prev.size() < cur.size()) prev.resize(cur.size(), NA_REAL);
  for (size_t i = 0; i < cur.size(); ++i) {
    if (std::isfinite(cur[i])) prev[i] = cur[i];
  }
}

bool rxEtaDistSpreadOk(const std::vector<double> &w, double lo, double hi,
                       double *sdOut, const std::vector<double> *wt = nullptr);

// Closed-form M-step for a Gaussian copula's correlation from paired latent
// draws.  `w` weights as above.
double rxEtaDistCorPost(const std::vector<double> &z1,
                        const std::vector<double> &z2);
double rxEtaDistCorSpearman(const std::vector<double> &z1,
                            const std::vector<double> &z2);
double rxEtaDistCorMleW(const std::vector<double> &w1,
                        const std::vector<double> &w2,
                        const std::vector<double> *w);
double rxEtaDistCorMle(const std::vector<double> &w1,
                       const std::vector<double> &w2);

#endif // __ETADISTFAM_H__
