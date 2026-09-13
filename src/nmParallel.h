#ifndef __NM_PARALLEL_H__
#define __NM_PARALLEL_H__
// One place for the conventions a per-subject parallel region has to follow.
//
// Every such region here owes rxode2 the same three things, and each was being
// assembled by hand at ~40 call sites:
//
//   1. POSITION -> SUBJECT ID.  rx->ordId is the solve order; a loop counts
//      positions and everything it touches -- the solve and the per-subject
//      arrays it writes -- must use the id.  Using the position for one and the
//      id for the other silently swaps subjects.
//   2. THREAD ID.  rxode2 and nlmixr2est each link their own libgomp on
//      Windows, so rxode2's omp_get_thread_num() inside a solve returns 0 for
//      every worker and collapses its per-thread buffers onto slot 0.  Our real
//      thread id has to be handed over around the solve.
//   3. EXCEPTIONS.  One escaping the region is uncatchable across the OpenMP
//      boundary and calls std::terminate, taking the R session with it.
//
// The thread id is set once per thread here rather than once per iteration: it
// cannot change within a thread inside one region, so the `parallel` and the
// `for` are separated to hoist it out.
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
// Requires rxode2ptr.h (setRxThreadId, getOrdId, getRxNsim) to be included first.
#include "rxomp.h"
#include <cstdint>  // int64_t in nmOrdId

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
  // a guarantee from the other side of the package boundary, so do not rely on
  // it for a signed multiply here, where overflow would be undefined.
  const int64_t nall = (int64_t)getRxNsub(rxIn) * (int64_t)getRxNsim(rxIn);
  // pos is bounded by n, so this also keeps getOrdId()'s unchecked
  // rx->ordId[pos] in range: it is only reached when n IS that length.
  return ((int64_t)n == nall) ? getOrdId(rxIn, pos) : pos + 1;
}

// Flags.  Both exist because a region can need them, not as style knobs:
//
//   nmStatic -- the thread -> iteration assignment must be the SAME every run.
//     Needed when the body accumulates into a per-thread partial that is later
//     summed in thread order: with a dynamic schedule the partition shifts run
//     to run, the partial sums differ, and the reduction differs in its last
//     bits.  Costs load balancing, so do not reach for it otherwise.
//   nmNoMap -- do not walk rx->ordId.  Needed for the same reason: rx->ordId is
//     re-sorted by accumulated solve time, so it varies run to run, and a body
//     whose result depends on WHICH subjects a thread got must not see it.
//
// A region that needs either of these is a region whose answer depends on the
// partition, which is worth saying out loud at the call site.
// A distinct type, not a bare unsigned: with an unsigned the 6-argument
// flags call and the 6-argument error-sink call are both viable and the
// compiler picks the wrong one.
struct nmFlags {
  unsigned v;
  explicit constexpr nmFlags(unsigned x = 0u) : v(x) {}
};
static constexpr nmFlags nmDefault{0u}, nmStatic{1u}, nmNoMap{2u};
static constexpr nmFlags operator|(nmFlags a, nmFlags b) { return nmFlags(a.v | b.v); }
static constexpr bool operator&(nmFlags a, nmFlags b) { return (a.v & b.v) != 0u; }

// Run `body(id)` once per subject.  `n` is a subject count (getRxNsub(rx) or a
// value derived from it); `par` is whether this region will actually fork.
template <class F>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    nmFlags flags, F&& body) {
  const bool mapIt = !(flags & nmNoMap);
#ifdef _OPENMP
  if (par) {
    if (flags & nmStatic) {
#pragma omp parallel num_threads(cores)
      {
        setRxThreadId(omp_get_thread_num());
#pragma omp for schedule(static)
        for (int pos = 0; pos < n; ++pos) {
          body(mapIt ? nmOrdId(rxIn, pos, n) - 1 : pos);
        }
        setRxThreadId(-1);
      }
    } else {
#pragma omp parallel num_threads(cores)
      {
        setRxThreadId(omp_get_thread_num());
#pragma omp for schedule(dynamic)
        for (int pos = 0; pos < n; ++pos) {
          body(mapIt ? nmOrdId(rxIn, pos, n) - 1 : pos);
        }
        setRxThreadId(-1);
      }
    }
    return;
  }
#else
  (void)cores; (void)par;
#endif
  for (int pos = 0; pos < n; ++pos) body(pos);
}

template <class F>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    F&& body) {
  nmForEachSubject(rxIn, n, cores, par, nmDefault, body);
}

// As above, but an exception escaping `body` is caught and reported to
// `onError(id)` rather than reaching the OpenMP boundary.  `onError` must not
// throw and must not touch R.
//
// Only the parallel path catches.  Serial, an exception has a caller that can
// still see it, and swallowing it there would turn a real error into one
// subject quietly marked bad -- so it propagates, as it did before.
template <class F, class E>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    nmFlags flags, F&& body, E&& onError) {
  if (!par) {
    nmForEachSubject(rxIn, n, cores, false, flags, body);
    return;
  }
  nmForEachSubject(rxIn, n, cores, true, flags, [&](int id) {
    try {
      body(id);
    } catch (...) {
      onError(id);
    }
  });
}

template <class F, class E>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    F&& body, E&& onError) {
  nmForEachSubject(rxIn, n, cores, par, nmDefault, body, onError);
}

#endif // __NM_PARALLEL_H__
