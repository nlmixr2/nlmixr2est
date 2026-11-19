#ifndef __LIKCONTRIBUTIL_H__
#define __LIKCONTRIBUTIL_H__
// Internal (not part of the contributor ABI) helpers shared by every objective
// that cycles the external likelihood-contribution registry: likInner0 and
// vaeDecoderPxzCore in inner.cpp, and the population objective in nlm.cpp.
#include <float.h>
#include <atomic>
#include "censEst.h"
#include "../inst/include/nlmixr2estLikContrib.h"

// #1051: the analytic outer gradient (foceiControl(fast=TRUE)) re-derives
// d(objective)/d(theta) from MODEL sensitivities alone, so a contributor that
// CHANGES the objective -- writes llik, or d(LL)/d(eta), which moves eta* off
// the base problem's stationary point -- makes it wrong.  A pure observer
// writes neither and stays exact.  Which one a bundle is cannot be asked of it,
// so it is observed here, the single place every objective cycles the registry,
// and read by analyticOuterGrad().  Defined in inner.cpp; reset per fit by
// foceiOuter() and whenever the registry changes.
extern std::atomic<int> _nlmixrContribSeen;     // the obs hook has run at least once
extern std::atomic<int> _nlmixrContribChanged;  // ... and something wrote back

// Exact Gaussian cotangents d(LL)/d(f) and d(LL)/d(r), honoring censoring the
// same way the base objective does: dCensNormal1 chains the uncensored score
// through (df, dr), so (df=1, dr=0) yields d(LL)/d(f) and (df=0, dr=1) yields
// d(LL)/d(r).  Uncensored records return the uncensored slope unchanged; M2 adds
// the censored adjustment, M3/M4 replace it.
static inline void nlmixrLikContribGaussCotan(double cens, double dv, double limit,
                                              double f, double r,
                                              double *dLLdf, double *dLLdr) {
  // r == 0 guard.  This is censEst.h's _safe_zero, spelled out because that
  // header #undef's the macro at its end.  DBL_EPSILON (not sqrt(DBL_EPSILON))
  // is deliberate: dCensNormal1 applies _safe_zero to r internally, so the dll
  // handed to it must use the same guard or the chain rule is inconsistent --
  // and it matches what likInner0's base objective has always used.
  double rz = (r == 0.0) ? DBL_EPSILON : r;
  double err = f - dv;
  *dLLdf = dCensNormal1(cens, dv, limit, -err / rz, f, r, 1.0, 0.0);
  *dLLdr = dCensNormal1(cens, dv, limit,
                        0.5 * err * err / (rz * rz) - 0.5 / rz, f, r, 0.0, 1.0);
}

// Cycle the registry for one observation; returns the extra log-likelihood the
// contributors added.  dLLdEta (length neta, may be NULL when neta == 0) is
// zeroed first and comes back holding the contributors' extra d(LL)/d(eta).
static inline double nlmixrLikContribObs1(int id, int k, int neta,
                                          double f, double dv, double r,
                                          double dLLdf, double dLLdr,
                                          const double *dfdEta, double *dLLdEta) {
  for (int q = 0; q < neta; ++q) dLLdEta[q] = 0.0;
  double llAdd = 0.0;
  nlmixrLikObs o;
  o.id = id; o.k = k; o.neta = neta;
  o.f = f; o.dv = dv; o.r = r; o.dLL_df = dLLdf; o.dLL_dr = dLLdr;
  o.df_deta = dfdEta; o.llik = &llAdd; o.dLL_deta = dLLdEta;
  nlmixrLikContribObs(&o);
  // Relaxed set-once-to-1 from any thread (#1051); never cleared here.
  _nlmixrContribSeen.store(1, std::memory_order_relaxed);
  if (llAdd != 0.0) {
    _nlmixrContribChanged.store(1, std::memory_order_relaxed);
  } else {
    for (int q = 0; q < neta; ++q) {
      if (dLLdEta[q] != 0.0) {
        _nlmixrContribChanged.store(1, std::memory_order_relaxed);
        break;
      }
    }
  }
  return llAdd;
}

// Per-subject brackets; nobs/eta differ per caller so they stay parameters.
static inline void nlmixrLikContribBeginSubj(int id, int neta, int nobs,
                                             const double *eta) {
  nlmixrLikSubj s;
  s.id = id; s.neta = neta; s.nobs = nobs; s.eta = eta;
  nlmixrLikContribBegin(&s);
}

static inline void nlmixrLikContribEndSubj(int id, int neta, int nobs,
                                           const double *eta) {
  nlmixrLikSubj s;
  s.id = id; s.neta = neta; s.nobs = nobs; s.eta = eta;
  nlmixrLikContribEnd(&s);
}
#endif
