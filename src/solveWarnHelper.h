#ifndef __NLMIXR2EST_SOLVE_WARN_HELPER_H__
#define __NLMIXR2EST_SOLVE_WARN_HELPER_H__

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Flush rxode2's aggregated ODE-solve warnings (currently intdy
   window-misses). Resolved lazily on first call so a stale rxode2 that
   predates the symbol degrades to a no-op rather than aborting the
   optimization. Call sites: end of FOCEi outer iteration print and end of
   SAEM per-iter print. */
static inline void nmFlushRxSolveWarn(int maxIds) {
  typedef void (*rxSolveWarnFlush_t)(int);
  static rxSolveWarnFlush_t fn = NULL;
  static int triedLoad = 0;
  if (!triedLoad) {
    /* R_GetCCallable errors if the symbol is missing; guard via
       R_FindSymbol's softer mechanism is not applicable here because
       R_RegisterCCallable stores into a separate table. Trust the
       version-bumped `Imports: rxode2 (>= ...)` in DESCRIPTION to
       prevent a stale install reaching this call. */
    fn = (rxSolveWarnFlush_t) R_GetCCallable("rxode2", "rxSolveWarnFlush");
    triedLoad = 1;
  }
  if (fn != NULL) fn(maxIds);
}

#endif
