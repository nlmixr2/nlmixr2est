#ifndef __LIKCONTRIBUTIL_H__
#define __LIKCONTRIBUTIL_H__
// Internal (not part of the contributor ABI) helpers shared by every objective
// that cycles the external likelihood-contribution registry: likInner0 and
// vaeDecoderPxzCore in inner.cpp, and the population objective in nlm.cpp.
#include <float.h>
#include "censEst.h"
#include "../inst/include/nlmixr2estLikContrib.h"

// Exact Gaussian cotangents d(LL)/d(f) and d(LL)/d(r), honoring censoring the
// same way the base objective does: dCensNormal1 chains the uncensored score
// through (df, dr), so (df=1, dr=0) yields d(LL)/d(f) and (df=0, dr=1) yields
// d(LL)/d(r).  Uncensored records return the uncensored slope unchanged; M2 adds
// the censored adjustment, M3/M4 replace it.
static inline void nlmixrLikContribGaussCotan(double cens, double dv, double limit,
                                              double f, double r,
                                              double *dLLdf, double *dLLdr) {
  // censEst.h #undef's its _safe_zero macro at the end of the header, so spell
  // the same zero guard out here.
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
