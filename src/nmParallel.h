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
  return (n == getRxNsub(rxIn) * getRxNsim(rxIn)) ? getOrdId(rxIn, pos) : pos + 1;
}

// Run `body(id)` once per subject.  `n` is a subject count (getRxNsub(rx) or a
// value derived from it); `par` is whether this region will actually fork.
template <class F>
static inline void nmForEachSubject(rx_solve *rxIn, int n, int cores, bool par,
                                    F&& body) {
#ifdef _OPENMP
  if (par) {
#pragma omp parallel num_threads(cores)
    {
      setRxThreadId(omp_get_thread_num());
#pragma omp for schedule(dynamic)
      for (int pos = 0; pos < n; ++pos) {
        body(nmOrdId(rxIn, pos, n) - 1);
      }
      setRxThreadId(-1);
    }
    return;
  }
#else
  (void)cores; (void)par;
#endif
  for (int pos = 0; pos < n; ++pos) body(pos);
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
                                    F&& body, E&& onError) {
  if (!par) {
    nmForEachSubject(rxIn, n, cores, false, body);
    return;
  }
  nmForEachSubject(rxIn, n, cores, true, [&](int id) {
    try {
      body(id);
    } catch (...) {
      onError(id);
    }
  });
}

#endif // __NM_PARALLEL_H__
