#ifndef __ETADIST_KERNEL_H__
#define __ETADIST_KERNEL_H__
// Metropolis kernels for a declared eta sampled ON ITS OWN SCALE.
//
// saem's existing kernels are already the saemix three: (1) an independence
// proposal from the prior, (2) a random walk with the prior's covariance, (3) a
// coordinate-wise random walk.  They are written for a GAUSSIAN prior, which is
// correct today because the expansion leaves the sampled column standard normal
// and puts the family in a decoder.  These are the same three kernels for the
// case where the sampled quantity IS the declared eta.
//
// Two things change and nothing else does:
//
//   * kernel 1 proposes from the FAMILY -- inverse-CDF through rxEtaDistQ(),
//     which every declarable family already has, so one code path covers all
//     22 and no per-family generator is needed.  The proposal is the prior, so
//     the prior terms cancel from the ratio exactly as they do now.
//
//   * kernels 2 and 3 walk on the BIJECTED scale u (etaDistEtaScale.h), so a
//     positive or bounded eta cannot be proposed outside its support.  A
//     symmetric walk in u targets p(eta(u)) |d eta/du|, so the log-Jacobian
//     enters the ratio.  Leaving it out is not an inefficiency, it is the wrong
//     stationary distribution -- biased toward wherever the map compresses.
//
// The likelihood is a callback so these are testable without a model: with a
// flat likelihood the chain must reproduce the family's own moments, which is
// the property that says the kernel is right independently of everything the
// estimator does around it.

#include "etaDistEtaScale.h"
#include <functional>

// -log p(eta) for the declared family.  The sign convention matches saem's
// U_phi, which is a NEGATIVE log density, so the existing `deltu` assembly and
// its `deltu < -log(u)` test apply unchanged.
static inline double rxEtaDistUPhi(int fam, double x, const double *a) {
  double v = rxEtaDistLogD(fam, x, a);
  return R_finite(v) ? -v : R_PosInf;
}

// One kernel-1 (independence-from-prior) step.
//
// `u01` is a uniform for the proposal, `uAcc` the acceptance uniform.  Returns
// the accepted eta.  `negLL(eta)` is -log p(y | eta); its value at the current
// point is passed in and updated on acceptance so a caller sweeping many steps
// pays one likelihood per proposal, not two.
static inline double rxEtaDistKern1(int fam, const double *a, double etaCur,
                                    double u01, double uAcc,
                                    const std::function<double(double)> &negLL,
                                    double *negLLCur) {
  double etaC = rxEtaDistQ(fam, u01, a);
  if (!R_finite(etaC)) return etaCur;
  double lc = negLL(etaC);
  if (!R_finite(lc)) return etaCur;
  // proposal IS the prior: the prior terms cancel, leaving the data alone
  double deltu = lc - *negLLCur;
  if (!R_finite(*negLLCur) || deltu < -std::log(uAcc)) {
    *negLLCur = lc;
    return etaC;
  }
  return etaCur;
}

// One kernel-2/3 (random walk on the bijected scale) step.
//
// `z` is a standard normal, `s` the step scale.  The Jacobian terms are what
// make this correct rather than merely in-support.
static inline double rxEtaDistKern2(int fam, const double *a, double etaCur,
                                    double z, double s, double uAcc,
                                    const std::function<double(double)> &negLL,
                                    double *negLLCur) {
  double uCur = rxEtaDistToU(fam, etaCur, a);
  if (!R_finite(uCur)) return etaCur;
  double uC = uCur + s*z;
  double etaC = rxEtaDistFromU(fam, uC, a);
  if (!R_finite(etaC)) return etaCur;
  double lc = negLL(etaC);
  if (!R_finite(lc)) return etaCur;
  double jc = rxEtaDistLogJac(fam, uC, a), j0 = rxEtaDistLogJac(fam, uCur, a);
  if (!R_finite(jc) || !R_finite(j0)) return etaCur;
  double pc = rxEtaDistUPhi(fam, etaC, a), p0 = rxEtaDistUPhi(fam, etaCur, a);
  if (!R_finite(pc)) return etaCur;
  // U = -log density, so the Jacobian enters with a MINUS: the target in u is
  // p(eta(u))|d eta/du|, i.e. U_u = U_eta - logJac.
  double deltu = (lc - *negLLCur) + (pc - p0) - (jc - j0);
  if (!R_finite(p0) || !R_finite(*negLLCur) || deltu < -std::log(uAcc)) {
    *negLLCur = lc;
    return etaC;
  }
  return etaCur;
}

#endif // __ETADIST_KERNEL_H__
