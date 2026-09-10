#ifndef __LOGSUMEXP_H__
#define __LOGSUMEXP_H__
#include <cmath>
#include <limits>
#include <Rmath.h>

// Weighted log-sum-exp, max-shifted:
//
//   log( sum_k w[k] * exp(v[k]) )
//
// The shift is what keeps a term that underflows from taking the whole sum with
// it, which is the entire reason a mixture marginal is computed in log space.
// `w == NULL` means unit weights.
//
// `skipNonFinite` distinguishes the two callers, and they genuinely differ:
//
//   TRUE  -- a mixture marginal.  A component that failed to solve contributes
//            nothing rather than poisoning the sum.  Returns R_NegInf when NO
//            component is finite, so the caller can charge its own penalty;
//            mixture code must not let that read as a likelihood of 1.
//   FALSE -- the importance-sampling proposal kernel, which has always let a
//            NaN propagate.  Quietly dropping one there would turn a loud
//            failure into a plausible number.
//
// When `wOut` is non-NULL it receives the normalized weights
// w[k]*exp(v[k]) / sum_j w[j]*exp(v[j]) -- the posterior responsibilities.
// They are left at 0 when the sum is not finite.
//
// static inline, in a header, because impPropLogKernelRecip() calls this once
// per importance sample; a cross-translation-unit call there is measurable.
static inline double rxLogSumExpW(const double *v, const double *w, int n,
                                  double *wOut, bool skipNonFinite) {
  double mx = -std::numeric_limits<double>::infinity();
  for (int k = 0; k < n; ++k) {
    if (wOut != NULL) wOut[k] = 0.0;
    if (skipNonFinite && !R_FINITE(v[k])) continue;
    double lp = v[k] + ((w == NULL) ? 0.0 : std::log(std::max(1e-300, w[k])));
    if (lp > mx) mx = lp;          // NaN compares false, so it never wins here
  }
  if (!R_FINITE(mx)) return R_NegInf;
  double se = 0.0;
  for (int k = 0; k < n; ++k) {
    if (skipNonFinite && !R_FINITE(v[k])) continue;
    double e = std::exp(v[k] + ((w == NULL) ? 0.0 : std::log(std::max(1e-300, w[k]))) - mx);
    se += e;                        // a NaN term reaches the result when
    if (wOut != NULL) wOut[k] = e;  // skipNonFinite is FALSE, as it always did
  }
  if (wOut != NULL && se > 0.0) {
    for (int k = 0; k < n; ++k) wOut[k] /= se;
  }
  return mx + std::log(se);
}

// Mixture marginal: weighted, and a failed component contributes nothing.
static inline double rxLogSumExpMix(const double *v, const double *w, int n,
                                    double *wOut) {
  return rxLogSumExpW(v, w, n, wOut, true);
}

// Unweighted log( sum_k exp(v[k]) ), NaN-propagating.
static inline double rxLogSumExp(const double *v, int n) {
  return rxLogSumExpW(v, NULL, n, NULL, false);
}

#endif // __LOGSUMEXP_H__
