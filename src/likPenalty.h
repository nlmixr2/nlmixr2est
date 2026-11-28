#ifndef __NLMIXR2EST_LIKPENALTY_INTERNAL_H__
#define __NLMIXR2EST_LIKPENALTY_INTERNAL_H__
// Internal dispatch for the likelihood penalty ABI
// (inst/include/nlmixr2estLikPenalty.h).  Defined in likPenalty.cpp.
#include "../inst/include/nlmixr2estLikPenalty.h"

#ifdef __cplusplus
extern "C" {
#endif
  // number of active penalties with a subject / population hook
  int nlmixrPenNSubject(void);
  int nlmixrPenNPop(void);

  // Subject terms, all active penalties in series.  Returns the extra LL;
  // dLLdEta (may be NULL) is ADDED to.  Thread-safe (no shared writes except the
  // relaxed "changed" flag).
  double nlmixrPenSubjAll(int id, int neta, int nobs, const double *eta,
                          int ntheta, const double *theta, double *dLLdEta);

  // Population terms on the -2LL scale.  grad (may be NULL) is ADDED to; *gradOk
  // (may be NULL) is set to 0 when some penalty cannot supply a gradient.
  double nlmixrPenPopAll(int ntheta, const double *theta, int neta,
                         const double *omega, int nsub, double *grad, int *gradOk);

  // d(population terms)/d(omega) by central differences; dOmega [neta*neta] is
  // overwritten (symmetric, zero for omega-free penalties).
  void nlmixrPenPopOmegaGrad(int ntheta, const double *theta, int neta,
                             const double *omega, int nsub, double *dOmega);

  // 1 when a subject term has written a non-zero value since activation: such a
  // term moves eta* off the base problem, so the analytic outer gradient declines.
  int nlmixrPenSubjChanged(void);
#ifdef __cplusplus
}
#endif
#endif
