#ifndef __NM_MCMC_RNG_H__
#define __NM_MCMC_RNG_H__
// Shared MCMC / importance-sampling RNG guard (used by saem.cpp and imp.cpp).
//
// Every rxode2 solver re-seeds the threefry engine per subject mid-solve
// (setSeedEng1(getRxSeed1() + id)), so an inner-likelihood evaluation clobbers
// whatever seed the sampling block set.  A sampling block records its seed via
// nmSetSeedEng1() and, right after the inner-likelihood call, restores it
// (nmRestoreMcmcSeed() / nmRngGuard()) so the solve's per-subject re-seed can
// never carry into the sampling draws.
//
// Requires rxode2ptr.h (for setSeedEng1) to be included first.

// One shared instance per thread across TUs (C++17 inline variable).
inline thread_local uint32_t _nmMcmcEngSeed = 0u;

// Set the engine seed AND record it as the current sampling-block seed.
static inline void nmSetSeedEng1(uint32_t s) {
  _nmMcmcEngSeed = s;
  setSeedEng1(s);
}

// Restore the recorded sampling-block seed (undo an inner solve's re-seed).
static inline void nmRestoreMcmcSeed() {
  setSeedEng1(_nmMcmcEngSeed);
}

// Evaluate an inner likelihood, then restore the sampling-block seed.
template <class F>
static inline auto nmRngGuard(F&& fn) -> decltype(fn()) {
  auto _r = fn();
  nmRestoreMcmcSeed();
  return _r;
}

#endif // __NM_MCMC_RNG_H__
