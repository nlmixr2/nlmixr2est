#ifndef __IMPQRNG_H__
#define __IMPQRNG_H__
// Scrambling for the QRPEM Sobol point set (src/imp.cpp).
//
// A raw Sobol sequence is deterministic, so the only randomization the
// importance sampler had was a Cranley-Patterson shift: u -> frac(u + s).  A CP
// shift randomizes the point set but does not break the correlation structure
// between the high-order dimensions of the sequence, which is where a Sobol set
// degrades first -- so it helps least exactly where a fit has many random
// effects.  Scrambling permutes the digits themselves and does break it.
//
// Two families, both applied to the ENGINE'S OUTPUT bits rather than to its
// direction numbers (boost does not expose those).  That is legitimate for
// both: a Sobol point is a linear (over GF(2)) image of the direction numbers,
// so a linear scramble commutes with generation, and a nested digit
// permutation is defined on the digits of the point regardless of how it
// was produced.
//
// Both are keyed by a seed derived ARITHMETICALLY from the fit's impSeed and
// the (iteration, subject, dimension) indices -- never by drawing from the
// threefry engine.  So scrambling consumes no RNG draws, leaves every other
// draw stream in the kernel untouched, and stays reproducible and independent
// of the thread count.
//
// Resolution is 32 bits.  The unscrambled path keeps its 53-bit doubles
// untouched (impSobolU0 is unchanged); a scrambled coordinate is built from
// the top 32 bits of the same 64-bit engine output, which separates far more
// points than any usable isample.

#include <cstdint>

// ---- seed mixing -----------------------------------------------------------
// A fit's scramble seed must fold in EVERY index that identifies it.  A bare
// sum collides whenever two index pairs sum alike (the convention the SAEM
// seeding comment in CLAUDE.md records as a real past bug), so mix
// multiplicatively.  This is the SplitMix64 finalizer.
static inline uint64_t impQrngMix(uint64_t z) {
  z += 0x9E3779B97F4A7C15ULL;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

// Per-(subject, iteration, dimension) scramble key.
static inline uint32_t impQrngSeed(uint32_t base, int id, int iter, int dim) {
  uint64_t z = (uint64_t)base;
  z = impQrngMix(z ^ (0x1000193ULL * (uint64_t)(uint32_t)id));
  z = impQrngMix(z ^ (0x01000193ULL * (uint64_t)(uint32_t)iter));
  z = impQrngMix(z ^ (0x85EBCA6BULL * (uint64_t)(uint32_t)dim));
  return (uint32_t)(z >> 32);
}

static inline uint32_t impQrngReverseBits(uint32_t x) {
  x = (x << 16) | (x >> 16);
  x = ((x & 0x00FF00FFu) << 8) | ((x & 0xFF00FF00u) >> 8);
  x = ((x & 0x0F0F0F0Fu) << 4) | ((x & 0xF0F0F0F0u) >> 4);
  x = ((x & 0x33333333u) << 2) | ((x & 0xCCCCCCCCu) >> 2);
  x = ((x & 0x55555555u) << 1) | ((x & 0xAAAAAAAAu) >> 1);
  return x;
}

// ---- Owen (nested uniform) scrambling --------------------------------------
// Burley (2020), "Practical Hash-based Owen Scrambling", JCGT 9(4).  The
// Laine-Karras permutation on the bit-reversed value is a hash-based
// approximation to a nested uniform scramble: each bit is flipped according to
// a function of the bits above it, which is exactly Owen's construction.
static inline uint32_t impQrngLaineKarras(uint32_t x, uint32_t seed) {
  x += seed;
  x ^= x * 0x6C50B47Cu;
  x ^= x * 0xB82F1E52u;
  x ^= x * 0xC7AFE638u;
  x ^= x * 0x8D22F6E6u;
  return x;
}

static inline uint32_t impOwenScramble(uint32_t x, uint32_t seed) {
  x = impQrngReverseBits(x);
  x = impQrngLaineKarras(x, seed);
  return impQrngReverseBits(x);
}

// ---- linear matrix scrambling ----------------------------------------------
// Matousek (1998) / Tezuka: a random non-singular lower-triangular binary
// matrix L applied to the digit vector, followed by a digital shift.  Named
// "lms" and documented as linear matrix scrambling rather than after any one
// vendor's variant -- Phoenix's "Tezuka-Faure" is not published in enough
// detail to claim we reproduce it.
//
// Row i of L is (unit diagonal | random bits strictly below it), so L is
// unit lower-triangular and therefore always invertible; output bit i is the
// parity of (v & row_i).  Bit 0 here is the MOST significant digit, matching
// the radical-inverse convention.
static inline uint32_t impLmsScramble(uint32_t x, uint32_t seed) {
  uint32_t out = 0u;
  uint32_t shift = impQrngMix(((uint64_t)seed << 1) | 1ULL) >> 32;
  for (int i = 0; i < 32; ++i) {
    // row i: diagonal bit at position i, random bits at positions > i
    uint32_t rnd = (uint32_t)(impQrngMix((uint64_t)seed * 0x2545F4914F6CDD1DULL +
                                         (uint64_t)i) >> 32);
    uint32_t row = (0x80000000u >> i);
    if (i < 31) row |= (rnd & (0x7FFFFFFFu >> i));
    uint32_t v = x & row;
    // parity of v
    v ^= v >> 16; v ^= v >> 8; v ^= v >> 4; v ^= v >> 2; v ^= v >> 1;
    if (v & 1u) out |= (0x80000000u >> i);
  }
  return out ^ shift;
}

// Scramble modes, matching impmapControl(qrScramble=).
enum impQrScrambleType { impQrScrambleNone = 0, impQrScrambleOwen = 1, impQrScrambleLms = 2 };

// One scrambled coordinate as a uniform in (0,1).  `v` is the engine's raw
// 64-bit output; only its top 32 bits are used.
static inline double impQrScrambleU(uint64_t v, int type, uint32_t seed) {
  uint32_t x = (uint32_t)(v >> 32);
  if (type == impQrScrambleOwen) x = impOwenScramble(x, seed);
  else if (type == impQrScrambleLms) x = impLmsScramble(x, seed);
  // half-step offset, matching impSobolU0's convention, so no coordinate is
  // ever exactly 0 or 1
  return std::ldexp((double)x + 0.5, -32);
}

#endif // __IMPQRNG_H__
