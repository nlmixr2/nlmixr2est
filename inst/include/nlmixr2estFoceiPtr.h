#ifndef NLMIXR2EST_FOCEI_PTR_H
#define NLMIXR2EST_FOCEI_PTR_H
// Downstream (caller) side of the FOCEi conditional-likelihood C API (#937):
// install nlmixr2est's entry points as function pointers from the small
// external-pointer table returned by _nlmixr2est_foceiPtrs()
// (nlmixr2est:::.nlmixr2estFoceiPtrs()).  Include this, instantiate the
// globals macro in ONE translation unit, and call iniNlmixr2estFocei(p) at
// .onLoad.  Requires <Rinternals.h> only (SEXP, R_ExternalPtrAddrFn,
// VECTOR_ELT): every quantity crosses the DSO boundary as double*/int,
// deliberately, so no Rcpp/Armadillo/Eigen/Stan type ever appears in more
// than one shared object.
//
// Lifecycle -- the table hands out RAW function pointers:
//  - call iniNlmixr2estFocei(p) on EVERY .onLoad, not just the first.  If
//    nlmixr2est is unloaded and reloaded, previously cached pointers address
//    an unmapped DLL, so the installer below always refreshes.
//  - a problem must be loaded with nlmixr2est::foceiLikLoad() before any
//    entry other than the version query is used; entries report "not loaded"
//    by return code rather than erroring.
//
// NONE of these entry points longjmp (Rf_error/stop) or let a C++ exception
// escape, with one documented exception: nlmixr2FoceiSetTheta calls back into
// R (the Omega rebuild) and therefore contains any R error with
// R_ToplevelExec, reporting it as a return code -- but it may only be called
// from the main R thread, never from a worker or a parallel region.  All
// other entries report failure by return code only, so they are safe to call
// from inside a sampler's log-density evaluation (past live destructors and
// an in-flight autodiff tape).

#ifdef __cplusplus
extern "C" {
#endif

  /* ---- 0. ABI version ------------------------------------------------- */
  /* Compare against NLMIXR2EST_FOCEI_API at load time; a mismatch means the
     downstream binary was built against a different table layout.          */
#define NLMIXR2EST_FOCEI_API 1
  typedef int (*nlmixr2FoceiApiVersion_t)(void);

  /* ---- 1. dimensions + capability/hazard flags ------------------------ */
  /* 0 on success, -1 if no likelihood system is loaded.  Any out pointer
     may be NULL.  flags bits:
       0x01 interaction==1 ("focei")  value/gradient self-consistent
       0x02 foceType==1    ("focep")  value and gradient are gradients of
                                      DIFFERENT functions (the live
                                      conditional R enters the value but its
                                      dR/deta term is absent from lp) --
                                      REFUSE for gradient-based samplers
       0x04 fo==1          (FO)       lp is never filled -- REFUSE
       0x08 some eta uses finite-difference event sensitivities (etaFD):
                                      ~1e-6 gradient noise -- refuse or
                                      accept knowingly
       0x10 mixture model             per-id indexing is component-major
                                      (id = component*nsub + subject, #955);
                                      the batch entries ACCEPT the expanded
                                      nid = nsub * nMix layout -- query nMix
                                      via the table's "nMix" entry
       0x20 the ODE solve method is thread-safe (cores>1 honored;
                                      otherwise cores is clamped to 1)
       0x40 the theta-sensitivity model is wired
                                      (nlmixr2FoceiCondThetaGrad available) */
  typedef int (*nlmixr2FoceiDims_t)(int *nid, int *neta, int *ntheta,
                                    int *npars, int *flags);

  /* ---- 2. population parameters --------------------------------------- */
  /* theta is the loaded problem's parameter vector (length npars) on the
     scale foceiLikLoad() installed.  MAIN R THREAD ONLY, never inside a
     parallel region: the Omega rebuild calls into R (contained with
     R_ToplevelExec).  0 ok, -1 not loaded, -2 bad length, -3 an R error was
     caught.                                                               */
  typedef int (*nlmixr2FoceiSetTheta_t)(const double *theta, int npars);

  /* ---- 3. the workhorse ------------------------------------------------ */
  /* eta   : nid x neta, ROW-MAJOR, subject order == foceiLikLoad handle's
             idLvl
     value : nid, log p(y_i | eta_i) -- the CONDITIONAL data log-likelihood,
             identical to foceiLikRun(type="cond"): no eta prior, and for
             Gaussian endpoints the per-observation -0.5 log(2 pi) constant
             is omitted (constant in every parameter)
     grad  : nid x neta, ROW-MAJOR, d/d(eta_i) log p(y_i | eta_i), computed
             inside nlmixr2est as Omega^-1 eta_i - fInd->lp with the SAME
             Omega^-1 the value used, so the caller never reconstructs the
             sign convention
     cores : clamped to the solve pool's size, and to 1 unless the solve
             method is thread safe
     Determinism: every call resets the sticky solve-tolerance relaxation and
     each subject's solve state, so the result is a pure function of
     (theta, eta) -- no eta caching, one full solve per subject per call.
     return: >=0 = number of subjects whose value was non-finite (their
             value is -Inf and their gradient row zeroed -- a sampler
             rejection, not an error); <0 hard failure (-1 not loaded,
             -2 shape/mixture, -3 exception).                              */
  typedef int (*nlmixr2FoceiCondBatch_t)(const double *eta, int nid, int neta,
                                         int cores, double *value, double *grad);

  /* ---- 4. Omega conditioning knob ------------------------------------- */
  /* Overwrite the loaded problem's Omega^-1 (and its Cholesky factor and
     0.5 log|Omega^-1|) GLOBALLY from a neta x neta symmetric positive-
     definite Omega^-1 (row-major; symmetry means storage order does not
     matter).  Purely numerical: fInd->lp is stored as -G + Omega^-1 eta, so
     keeping nlmixr2est's Omega commensurate with the caller's keeps the
     subtraction in nlmixr2FoceiCondBatch's gradient well conditioned.  The
     conditional value itself does not depend on Omega.  Main R thread only,
     outside any batch (deliberately NOT the thread_local OmegaScope, which
     workers would not see).  0 ok, -1 not loaded, -2 bad dimension,
     -3 not positive definite / decomposition failure.                     */
  typedef int (*nlmixr2FoceiSetOmegaInv_t)(const double *omegaInv, int neta);

  /* ---- 5/6. d/dtheta of the conditional at fixed eta ------------------- */
  /* idx receives the 0-based ntheta indices that carry an analytic
     d(f)/d(theta) forward sensitivity (the estimated non-mu structural and
     residual-error thetas).  Returns the count -- 0 whenever the sensitivity
     model is not wired for this load (flags & 0x40 clear) -- or -1 when no
     problem is loaded; pass idx==NULL (n ignored) to query the count only.
     Columns NOT listed are left at 0 by nlmixr2FoceiCondThetaGrad and are
     the caller's to fill from the mu-reference map:
     d/dtheta_p log p(y_i|eta_i) equals grad[i,k] (entry 3) for a theta p
     that mu-references eta_k.                                             */
  typedef int (*nlmixr2FoceiThetaSensIdx_t)(int *idx, int n);

  /* dTheta: nid x ntheta, ROW-MAJOR, d/d(theta) log p(y_i | eta_i) at fixed
     eta_i, zero outside the ThetaSensIdx columns.  The derivative is on the
     NATURAL theta scale (the forward sensitivities differentiate the model's
     THETA directly), independent of the scaling foceiLikLoad() installed --
     with foceiLikLoad(scale="natural") the two coincide; with the default
     FOCEi scaling the caller owns the chain-rule factor.
     Requires flags & 0x40 (the theta-sensitivity model wired at load).
     Same determinism discipline and return convention as
     nlmixr2FoceiCondBatch, plus -4 = sensitivity model not wired.         */
  typedef int (*nlmixr2FoceiCondThetaGrad_t)(const double *eta, int nid,
                                             int neta, int cores,
                                             double *dTheta);

  /* ---- 7. mixture component count (#955, additive) --------------------- */
  /* 1 for a non-mixture load, K for a K-component mixture, -1 not loaded.
     The batch entries take nid = nsub * nMix in component-major order
     (id = component*nsub + subject); the component-conditional value and
     gradients land in the matching rows.  NULL when the loaded nlmixr2est
     predates the entry (treat as refuse-mixtures).                        */
  typedef int (*nlmixr2FoceiNMix_t)(void);

  /* ---- 8. sampler iteration-print row (additive) ----------------------- */
  /* One evaluation's FULL natural-scale parameter vector (length npars,
     the setTheta layout) + objective, printed/recorded through the scale.h
     machinery at the cadence foceiLikIterPrintStart_() installed.  Plain C,
     non-throwing; 0 ok, -1 printing inactive, -2 bad shape, -3 exception.
     NULL when the loaded nlmixr2est predates the entry (skip printing).  */
  typedef int (*nlmixr2FoceiIterPrintRow_t)(const double *par, int npars,
                                            double objf);

  /* #958 fused tier-2 entry: value + d/d(eta) + d/d(theta) from ONE
     combined-model solve per subject.  Same contracts as condBatch +
     condThetaGrad merged: value/gradEta as condBatch (row-major
     nid x neta, -Inf + zeroed rows on rejection), dTheta as condThetaGrad
     (row-major nid x ntheta; sensitivity columns filled, mu-referenced
     columns zero for the caller's identity).  >=0 bad-subject count;
     -4 theta sensitivities not wired; -5 the loaded problem is not a
     combined build (foceiLikLoad(combSens=TRUE)).  NULL when the loaded
     nlmixr2est predates the entry.  */
  /* dims flag 0x80: the loaded problem is a combined-sensitivity build
     (foceiLikLoad(combSens=TRUE)) and the fused entry below works. */
  typedef int (*nlmixr2FoceiCondBatchThetaGrad_t)(const double *eta, int nid,
                                                  int neta, int cores,
                                                  double *value,
                                                  double *gradEta,
                                                  double *dTheta);

  extern nlmixr2FoceiApiVersion_t    nlmixr2FoceiApiVersionP;
  extern nlmixr2FoceiDims_t          nlmixr2FoceiDimsP;
  extern nlmixr2FoceiSetTheta_t      nlmixr2FoceiSetThetaP;
  extern nlmixr2FoceiCondBatch_t     nlmixr2FoceiCondBatchP;
  extern nlmixr2FoceiSetOmegaInv_t   nlmixr2FoceiSetOmegaInvP;
  extern nlmixr2FoceiThetaSensIdx_t  nlmixr2FoceiThetaSensIdxP;
  extern nlmixr2FoceiCondThetaGrad_t nlmixr2FoceiCondThetaGradP;
  extern nlmixr2FoceiNMix_t          nlmixr2FoceiNMixP;
  extern nlmixr2FoceiIterPrintRow_t  nlmixr2FoceiIterPrintRowP;
  extern nlmixr2FoceiCondBatchThetaGrad_t nlmixr2FoceiCondBatchThetaGradP;

  // Always refresh: a reloaded nlmixr2est hands back new addresses, and
  // keeping the first set seen would leave the caller calling into an
  // unloaded DLL.
  static inline SEXP iniNlmixr2estFocei0(SEXP p) {
    nlmixr2FoceiApiVersionP    = (nlmixr2FoceiApiVersion_t)    R_ExternalPtrAddrFn(VECTOR_ELT(p, 0));
    nlmixr2FoceiDimsP          = (nlmixr2FoceiDims_t)          R_ExternalPtrAddrFn(VECTOR_ELT(p, 1));
    nlmixr2FoceiSetThetaP      = (nlmixr2FoceiSetTheta_t)      R_ExternalPtrAddrFn(VECTOR_ELT(p, 2));
    nlmixr2FoceiCondBatchP     = (nlmixr2FoceiCondBatch_t)     R_ExternalPtrAddrFn(VECTOR_ELT(p, 3));
    nlmixr2FoceiSetOmegaInvP   = (nlmixr2FoceiSetOmegaInv_t)   R_ExternalPtrAddrFn(VECTOR_ELT(p, 4));
    nlmixr2FoceiThetaSensIdxP  = (nlmixr2FoceiThetaSensIdx_t)  R_ExternalPtrAddrFn(VECTOR_ELT(p, 5));
    nlmixr2FoceiCondThetaGradP = (nlmixr2FoceiCondThetaGrad_t) R_ExternalPtrAddrFn(VECTOR_ELT(p, 6));
    /* additive entries: read defensively so a downstream built against this
       header still installs cleanly from an OLDER nlmixr2est's 7-entry
       table (the missing capability reads as NULL) */
    nlmixr2FoceiNMixP = (Rf_xlength(p) > 7) ?
      (nlmixr2FoceiNMix_t) R_ExternalPtrAddrFn(VECTOR_ELT(p, 7)) : NULL;
    nlmixr2FoceiIterPrintRowP = (Rf_xlength(p) > 8) ?
      (nlmixr2FoceiIterPrintRow_t) R_ExternalPtrAddrFn(VECTOR_ELT(p, 8)) : NULL;
    nlmixr2FoceiCondBatchThetaGradP = (Rf_xlength(p) > 9) ?
      (nlmixr2FoceiCondBatchThetaGrad_t) R_ExternalPtrAddrFn(VECTOR_ELT(p, 9)) : NULL;
    return R_NilValue;
  }

#define iniNlmixr2estFoceiGlobals                                       \
  nlmixr2FoceiApiVersion_t    nlmixr2FoceiApiVersionP    = NULL;        \
  nlmixr2FoceiDims_t          nlmixr2FoceiDimsP          = NULL;        \
  nlmixr2FoceiSetTheta_t      nlmixr2FoceiSetThetaP      = NULL;        \
  nlmixr2FoceiCondBatch_t     nlmixr2FoceiCondBatchP     = NULL;        \
  nlmixr2FoceiSetOmegaInv_t   nlmixr2FoceiSetOmegaInvP   = NULL;        \
  nlmixr2FoceiThetaSensIdx_t  nlmixr2FoceiThetaSensIdxP  = NULL;        \
  nlmixr2FoceiCondThetaGrad_t nlmixr2FoceiCondThetaGradP = NULL;        \
  nlmixr2FoceiNMix_t          nlmixr2FoceiNMixP          = NULL;        \
  nlmixr2FoceiIterPrintRow_t  nlmixr2FoceiIterPrintRowP  = NULL;        \
  nlmixr2FoceiCondBatchThetaGrad_t nlmixr2FoceiCondBatchThetaGradP = NULL; \
  SEXP iniNlmixr2estFocei(SEXP p) { return iniNlmixr2estFocei0(p); }

#ifdef __cplusplus
}
#endif
#endif
