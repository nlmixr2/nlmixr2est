#ifndef NLMIXR2EST_NLM_PTR_H
#define NLMIXR2EST_NLM_PTR_H
// Downstream (caller) side of the nlm population-likelihood C API (#953):
// install nlmixr2est's entry points as function pointers from the small
// external-pointer table returned by _nlmixr2est_nlmPtrs()
// (nlmixr2est:::.nlmixr2estNlmPtrs()).  Include this, instantiate the
// globals macro in ONE translation unit, and call iniNlmixr2estNlm(p) at
// .onLoad.  Requires <Rinternals.h> only: every quantity crosses the DSO
// boundary as double*/int, deliberately, so no Rcpp/Armadillo/Eigen/Stan
// type ever appears in more than one shared object (same idiom as
// nlmixr2estFoceiPtr.h).
//
// Lifecycle:
//  - call iniNlmixr2estNlm(p) on EVERY .onLoad (a reloaded nlmixr2est hands
//    back new addresses).
//  - a problem must be loaded with nlmixr2est::nlmObjectiveSetup() before
//    any entry other than the version query; entries report "not loaded"
//    by return code.
//
// NONE of these entry points longjmp (Rf_error/stop) or let a C++ exception
// escape; failure is a return code, so they are safe to call from inside a
// sampler's log-density evaluation.
//
// Threading: nlmixr2NlmEval parallelizes over SUBJECTS internally (OpenMP,
// the loaded rxControl's cores).  The caller must NOT layer its own
// threading on top (run the sampler single-threaded, chains sequential) --
// the same contract as the FOCEi table, and the reason nlmixr2stan pins
// Stan to one thread.

#ifdef __cplusplus
extern "C" {
#endif

  /* ---- 0. ABI version ------------------------------------------------- */
#define NLMIXR2EST_NLM_API 1
  typedef int (*nlmixr2NlmApiVersion_t)(void);

  /* ---- 1. dimensions + capability/hazard flags ------------------------ */
  /* 0 on success, -1 if no nlm problem is loaded.  Any out pointer may be
     NULL.  flags bits:
       0x01 the gradient (sensitivity) model is loaded (solveType >= grad)
            -- REQUIRED for nlmixr2NlmEval's dTheta
       0x02 the installed scale is the identity (natural theta scale);
            without it the gradient is on the OPTIMIZER's scale -- refuse
            for samplers unless the caller owns the chain rule
       0x04 some theta's sensitivity is finite-differenced (thetaFD):
            FD noise on those columns -- refuse or accept knowingly       */
  typedef int (*nlmixr2NlmDims_t)(int *ntheta, int *nobs, int *flags);

  /* ---- 2. the evaluation ----------------------------------------------- */
  /* theta  : length ntheta, NATURAL scale (require flags & 0x02)
     value  : the summed population objective as the nlm path defines it:
              MINUS the log-likelihood, including its normalization
              constants (Gaussian endpoints carry the full 0.5 log(2 pi)
              term), i.e. log p(y | theta) = -value
     dTheta : length ntheta, d(value)/d(theta) -- same convention, so
              d/d(theta) log p = -dTheta
     Determinism: the entry resets the solver's bad-solve state before
     evaluating, and the nlm path re-solves every subject from its inits on
     every call, so the result is a pure function of theta.
     return: 0 ok; 1 the value or a gradient entry was non-finite (a
             sampler rejection); -1 not loaded; -2 gradient model not
             loaded (flags & 0x01 clear); -3 bad ntheta; -4 exception.   */
  typedef int (*nlmixr2NlmEval_t)(const double *theta, int ntheta,
                                  double *value, double *dTheta);

  extern nlmixr2NlmApiVersion_t nlmixr2NlmApiVersionP;
  extern nlmixr2NlmDims_t       nlmixr2NlmDimsP;
  extern nlmixr2NlmEval_t       nlmixr2NlmEvalP;

  static inline SEXP iniNlmixr2estNlm0(SEXP p) {
    nlmixr2NlmApiVersionP = (nlmixr2NlmApiVersion_t) R_ExternalPtrAddrFn(VECTOR_ELT(p, 0));
    nlmixr2NlmDimsP       = (nlmixr2NlmDims_t)       R_ExternalPtrAddrFn(VECTOR_ELT(p, 1));
    nlmixr2NlmEvalP       = (nlmixr2NlmEval_t)       R_ExternalPtrAddrFn(VECTOR_ELT(p, 2));
    return R_NilValue;
  }

#define iniNlmixr2estNlmGlobals                                  \
  nlmixr2NlmApiVersion_t nlmixr2NlmApiVersionP = NULL;           \
  nlmixr2NlmDims_t       nlmixr2NlmDimsP       = NULL;           \
  nlmixr2NlmEval_t       nlmixr2NlmEvalP       = NULL;           \
  SEXP iniNlmixr2estNlm(SEXP p) { return iniNlmixr2estNlm0(p); }

#ifdef __cplusplus
}
#endif
#endif
