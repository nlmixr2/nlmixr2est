#ifndef __NM_PARALLEL_CORE_H__
#define __NM_PARALLEL_CORE_H__
// The mechanics of a per-item parallel region, with NO rxode2 dependency.
//
// This is the half that every parallel region in the package has in common:
// the pragma, the choice of schedule, the decision whether to fork at all, and
// keeping an exception off the OpenMP boundary.  It is split out from
// nmParallel.h so that a translation unit with no rxode2 dependency can still
// share one implementation of all of that -- foceiGrad.cpp, npde.cpp and
// vaeEncoder.cpp do arma over pre-computed matrices and include not one rxode2
// header, and making them take rxode2's pointer API to gain a subject mapping
// that cannot apply to them would be a worse trade than the duplication it
// removes.
//
// nmParallel.h adds the rxode2 obligations on top: the rx->ordId position ->
// subject id mapping and the thread-id handoff.
//
// Requires rxomp.h.

#include "rxomp.h"

// Flags.  Both exist because a region can need them, not as style knobs:
//
//   nmStatic -- the thread -> iteration assignment must be the SAME every run.
//     Needed when the body accumulates into a per-thread partial that is later
//     summed in thread order: with a dynamic schedule the partition shifts run
//     to run, the partial sums differ, and the reduction differs in its last
//     bits.  Costs load balancing, so do not reach for it otherwise.
//   nmNoMap -- do not walk rx->ordId (see nmParallel.h).  Needed for the same
//     reason: rx->ordId is re-sorted by accumulated solve time, so it varies
//     run to run, and a body whose result depends on WHICH items a thread got
//     must not see it.  Meaningless here, where nothing maps; carried so the
//     two layers take one flag type.
//
// A region that needs either is a region whose answer depends on the
// partition, which is worth saying out loud at the call site.
//
// A distinct type, not a bare unsigned: with an unsigned the flags call and the
// error-sink call are both viable and overload resolution picks the wrong one.
struct nmFlags {
  unsigned v;
  explicit constexpr nmFlags(unsigned x = 0u) : v(x) {}
};
static constexpr nmFlags nmDefault{0u}, nmStatic{1u}, nmNoMap{2u};
static constexpr nmFlags operator|(nmFlags a, nmFlags b) { return nmFlags(a.v | b.v); }
static constexpr bool operator&(nmFlags a, nmFlags b) { return (a.v & b.v) != 0u; }

// Default per-thread scope: nothing.  nmParallel.h passes one that hands our
// real thread id to rxode2 for the life of the thread's slice of the region --
// once per thread rather than once per iteration, which is why the `parallel`
// and the `for` are separated below.
struct nmNoScope {};

// Run `body(i)` for i in [0, n).  `par` is whether this region will fork.
template <class Scope = nmNoScope, class F>
static inline void nmForEach(int n, int cores, bool par, nmFlags flags, F&& body) {
#ifdef _OPENMP
  if (par) {
    if (flags & nmStatic) {
#pragma omp parallel num_threads(cores)
      {
        Scope _nmScope;
#pragma omp for schedule(static)
        for (int i = 0; i < n; ++i) body(i);
      }
    } else {
#pragma omp parallel num_threads(cores)
      {
        Scope _nmScope;
#pragma omp for schedule(dynamic)
        for (int i = 0; i < n; ++i) body(i);
      }
    }
    return;
  }
#else
  (void)cores; (void)par; (void)flags;
#endif
  for (int i = 0; i < n; ++i) body(i);
}

template <class Scope = nmNoScope, class F>
static inline void nmForEach(int n, int cores, bool par, F&& body) {
  nmForEach<Scope>(n, cores, par, nmDefault, body);
}

// As above, but an exception escaping `body` is caught and reported to
// `onError(i)` rather than reaching the OpenMP boundary.  `onError` must not
// throw and must not touch R.
//
// Only the parallel path catches.  Serial, an exception still has a caller that
// can see it, and swallowing it there would turn a real error into one item
// quietly marked bad.
template <class Scope = nmNoScope, class F, class E>
static inline void nmForEach(int n, int cores, bool par, nmFlags flags,
                             F&& body, E&& onError) {
  if (!par) {
    nmForEach<Scope>(n, cores, false, flags, body);
    return;
  }
  nmForEach<Scope>(n, cores, true, flags, [&](int i) {
    try {
      body(i);
    } catch (...) {
      onError(i);
    }
  });
}

template <class Scope = nmNoScope, class F, class E>
static inline void nmForEach(int n, int cores, bool par, F&& body, E&& onError) {
  nmForEach<Scope>(n, cores, par, nmDefault, body, onError);
}

#endif // __NM_PARALLEL_CORE_H__
