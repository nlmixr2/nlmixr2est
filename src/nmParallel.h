#ifndef __NM_PARALLEL_H__
#define __NM_PARALLEL_H__
// The rxode2 obligations a per-subject parallel region owes, on top of the
// loop mechanics in nmParallelCore.h.
//
// Two of them, and each was being assembled by hand at ~40 call sites:
//
//   1. POSITION -> SUBJECT ID.  rx->ordId is the solve order; a loop counts
//      positions and everything it touches -- the solve and the per-subject
//      arrays it writes -- must use the id.  Using the position for one and the
//      id for the other silently swaps subjects.
//   2. THREAD ID.  rxode2 and nlmixr2est each link their own libgomp on
//      Windows, so rxode2's omp_get_thread_num() inside a solve returns 0 for
//      every worker and collapses its per-thread buffers onto slot 0.  Our real
//      thread id has to be handed over around the solve.  It is set once per
//      thread (nmRxThreadScope, constructed inside the region) rather than once
//      per iteration -- it cannot change within a thread inside one region.
//
// `body` receives the SUBJECT ID, never the position.  A caller that needs a
// mixture-component offset adds it inside the body, where the id is in hand.
//
// Deliberately knows nothing about op_focei or nlmOp -- those are private to
// inner.cpp and nlm.cpp, which is what forced this duplication in the first
// place (see odeSwap.h for the same constraint).
//
// This does NOT take the sortIds(rx,2)/sortIds(rx,0) pairing: callers hold that
// across several regions (FdPhaseStateGuard, NpInnerParallelScope), so it stays
// theirs.
//
// Requires rxode2ptr.h (setRxThreadId, getOrdId, getRxNsim, getRxNsub).

#include "nmParallelCore.h"
#include <cstdint>  // int64_t in nmOrdId

// Hands our real thread id to rxode2 for the life of this thread's slice.
struct nmRxThreadScope {
#ifdef _OPENMP
  nmRxThreadScope() { setRxThreadId(omp_get_thread_num()); }
  ~nmRxThreadScope() { setRxThreadId(-1); }
#endif
};

// Position -> subject id.  Returns a 1-BASED id, like rxode2's getOrdId().
//
// rx->ordId is a permutation of the nsub*nsim subject-SOLVES, so reading it
// means what a loop takes it to mean only when the loop covers all of them --
// `n` is the loop's own bound, and that is the test.  A loop over nsub when
// nsim > 1 would otherwise read an arbitrary subset of 1..nsub*nsim: some
// subjects visited twice, others skipped, and the per-subject arrays indexed
// past nsub.  Fall back to the data order there; unordered is slower, mixed up
// is wrong.
static inline int nmOrdId(rx_solve *rxIn, int pos, int n) {
  // 64-bit: nsub and nsim are both int, and their product is only guaranteed to
  // fit one because rxode2's sortIds() refuses a solve where it would not --
  // a guarantee from the far side of the package boundary, so do not rely on
  // it for a signed multiply here, where overflow would be undefined.
  const int64_t nall = (int64_t)getRxNsub(rxIn) * (int64_t)getRxNsim(rxIn);
  // pos is bounded by n, so this also keeps getOrdId()'s unchecked
  // rx->ordId[pos] in range: it is only reached when n IS that length.
  return ((int64_t)n == nall) ? getOrdId(rxIn, pos) : pos + 1;
}

// Run `body(id)` once per subject.  `n` is a subject count (getRxNsub(rx) or a
// value derived from it); `par` is whether this region will actually fork.
template <class F>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    nmFlags flags, F&& body) {
  const bool mapIt = !(flags & nmNoMap);
  nmForEach<nmRxThreadScope>(n, cores, par, flags, [&](int pos) {
    body(mapIt ? nmOrdId(rxIn, pos, n) - 1 : pos);
  });
}

template <class F>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    F&& body) {
  nmForEachSubject(rxIn, n, cores, par, nmDefault, body);
}

template <class F, class E>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    nmFlags flags, F&& body, E&& onError) {
  const bool mapIt = !(flags & nmNoMap);
  nmForEach<nmRxThreadScope>(n, cores, par, flags,
                             [&](int pos) { body(mapIt ? nmOrdId(rxIn, pos, n) - 1 : pos); },
                             [&](int pos) { onError(mapIt ? nmOrdId(rxIn, pos, n) - 1 : pos); });
}

template <class F, class E>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    F&& body, E&& onError) {
  nmForEachSubject(rxIn, n, cores, par, nmDefault, body, onError);
}

#endif // __NM_PARALLEL_H__
