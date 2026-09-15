#ifndef __NM_SEQ_SEED_H__
#define __NM_SEQ_SEED_H__
// Sequential threefry seeds with a closed form.
//
// A sampler lays its draws out as a fixed number of seeds per iteration and seeds
// item i of a step (a subject, chain row, observation or support point) with
// seed + offset + i right before that item's draws, where offset follows from the
// iteration and step alone.  Seeds are sequential and never reused, and a fit
// stopped and resumed at an iteration draws exactly what it would have.  Distinct
// threefry keys are independent streams, so nothing is hashed.
//
// Requires Rcpp, rxode2ptr.h and nmMcmcRng.h.

// The seed `offset` past the fit's seed.
static inline uint32_t nmSeqSeed(int seed, uint64_t offset) {
  return (uint32_t)seed + (uint32_t)offset;
}

// Seed the serial (thread 0) engine with item i past `offset`.
static inline void nmSeqSeedSet(int seed, uint64_t offset, uint64_t i) {
  setRxThreadId(0);
  nmSetSeedEng1(nmSeqSeed(seed, offset + i));
}

// Continue rxode2's own per-subject solve seeds right after a sampler's block:
// the sampler owns [seed, seed + reserved) and the solves take the sequence from
// seed + reserved, so no key is shared.  Called after setup, whose rxSolve_()
// advanced the sequence by the thread count.  It also keeps getRxSeed1() off
// R's RNG.  The caller runs inside rxWithSeed(), which restores the ambient seed.
static inline void nmSeqSeedStart(int seed, uint64_t reserved) {
  Rcpp::Function rxSetSeed = Rcpp::Environment::namespace_env("rxode2")["rxSetSeed"];
  rxSetSeed((double)nmSeqSeed(seed, reserved));
}

// Whether rxode2's seed sequence is in force (rxWithSeed(rxseed=) or rxSetSeed()).
static inline bool nmSeqSeedActive() {
  Rcpp::Function rxGetSeed = Rcpp::Environment::namespace_env("rxode2")["rxGetSeed"];
  return Rcpp::as<int>(rxGetSeed()) != -1;
}

#endif // __NM_SEQ_SEED_H__
