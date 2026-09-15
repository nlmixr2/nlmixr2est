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

// Restart rxode2's own per-subject solve seeds in the other half of the 32-bit
// range, after setup: a solve never shares a seed with a sampler's draw, and the
// setup solve's thread-count-dependent advance does not carry over.  It also
// keeps getRxSeed1() off R's RNG.  The caller runs inside rxWithSeed(), which
// restores the ambient seed afterward.
static inline void nmSeqSeedStart(int seed) {
  Rcpp::Function rxSetSeed = Rcpp::Environment::namespace_env("rxode2")["rxSetSeed"];
  rxSetSeed((double)nmSeqSeed(seed, 0x80000000u));
}

#endif // __NM_SEQ_SEED_H__
