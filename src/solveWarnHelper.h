#ifndef __NLMIXR2EST_SOLVE_WARN_HELPER_H__
#define __NLMIXR2EST_SOLVE_WARN_HELPER_H__

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Flush rxode2's aggregated ODE-solve warnings (intdy window-misses,
   lsoda nhnil, dop853 stiff/step-too-small). Resolved lazily on first
   call. The lookup runs inside R_ToplevelExec so a stale rxode2 that
   predates the symbol — common when CI installs the released rxode2
   from CRAN/RSPM against this branch — degrades to a no-op instead of
   crashing the fit with "function 'rxSolveWarnFlush' not provided by
   package 'rxode2'". Call sites: end of FOCEi outer iteration print
   and end of SAEM per-iter print. */

struct nmFlushLookup_ {
  void *fn; /* really a rxSolveWarnFlush_t but kept as void* in the struct */
};

static inline void nmFlushLookupBody_(void *data) {
  struct nmFlushLookup_ *ctx = (struct nmFlushLookup_ *)data;
  /* R_GetCCallable longjmps to top-level on missing symbol; capture
     under R_ToplevelExec so the longjmp turns into a return-false. */
  ctx->fn = (void *) R_GetCCallable("rxode2", "rxSolveWarnFlush");
}

static inline void nmFlushRxSolveWarn(int maxIds) {
  typedef void (*rxSolveWarnFlush_t)(int);
  static rxSolveWarnFlush_t fn = NULL;
  static int triedLoad = 0;
  if (!triedLoad) {
    struct nmFlushLookup_ ctx;
    ctx.fn = NULL;
    /* R_ToplevelExec returns TRUE iff the body completed without
       longjmping. On failure ctx.fn is whatever the body wrote before
       the error — we ignore it. */
    if (R_ToplevelExec(nmFlushLookupBody_, &ctx) && ctx.fn != NULL) {
      fn = (rxSolveWarnFlush_t) ctx.fn;
    }
    triedLoad = 1;
  }
  if (fn != NULL) fn(maxIds);
}

/* Push the real subject-id factor levels into rxode2's global solve state so
   the aggregated solve warnings flushed by nmFlushRxSolveWarn() are labelled
   with the user's ID instead of "Unknown".  FOCEi/SAEM pass a declassed
   data.frame/matrix to rxSolve_, so rxode2 never populates its factor table
   during estimation; we set it explicitly from the fit's idLvl right after
   setup.  Resolved lazily under R_ToplevelExec (same version-skew-safety as
   nmFlushRxSolveWarn): an rxode2 predating rxSetIdLvlFactors degrades to a
   no-op and the flush falls back to printing "internal #N". */

struct nmSetIdLvlLookup_ {
  void *fn; /* really a rxSetIdLvlFactors_t but kept as void* in the struct */
};

static inline void nmSetIdLvlLookupBody_(void *data) {
  struct nmSetIdLvlLookup_ *ctx = (struct nmSetIdLvlLookup_ *)data;
  ctx->fn = (void *) R_GetCCallable("rxode2", "rxSetIdLvlFactors");
}

static inline void nmSetIdLvlFactors(SEXP idLvl) {
  typedef void (*rxSetIdLvlFactors_t)(SEXP);
  static rxSetIdLvlFactors_t fn = NULL;
  static int triedLoad = 0;
  /* rxSetIdLvlFactors indexes idLvl with STRING_ELT, so only forward a
     character vector of subject-id levels; anything else is a no-op. */
  if (idLvl == R_NilValue || TYPEOF(idLvl) != STRSXP) return;
  if (!triedLoad) {
    struct nmSetIdLvlLookup_ ctx;
    ctx.fn = NULL;
    if (R_ToplevelExec(nmSetIdLvlLookupBody_, &ctx) && ctx.fn != NULL) {
      fn = (rxSetIdLvlFactors_t) ctx.fn;
    }
    triedLoad = 1;
  }
  if (fn != NULL) fn(idLvl);
}

#endif
