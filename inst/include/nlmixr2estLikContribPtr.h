#ifndef NLMIXR2EST_LIKCONTRIB_PTR_H
#define NLMIXR2EST_LIKCONTRIB_PTR_H
// Downstream (contributor) side of the likelihood-contribution registry: install
// nlmixr2est's register/remove entry points as function pointers from the small
// external-pointer table returned by _nlmixr2est_likContribPtrs()
// (nlmixr2est:::.nlmixr2estLikContribPtrs()).  Include this, instantiate the
// globals macro in one translation unit, and call iniNlmixr2estLikContrib(p) at
// .onLoad.  Requires <Rinternals.h> (SEXP, R_ExternalPtrAddrFn, VECTOR_ELT).
#include "nlmixr2estLikContrib.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef void (*nlmixrRegisterLikContrib_t)(const nlmixrLikContrib *c);
  typedef void (*nlmixrRemoveLikContrib_t)(const nlmixrLikContrib *c);
  typedef void (*nlmixrRegisterEmLik_t)(nlmixrEmLik_fn fn);
  typedef void (*nlmixrRemoveEmLik_t)(nlmixrEmLik_fn fn);
  typedef int  (*nlmixrHasLikContrib_t)(void);
  typedef void (*nlmixrSetInnerWeightFn_t)(nlmixrInnerWeight_fn fn);

  extern nlmixrRegisterLikContrib_t nlmixrRegisterLikContribP;
  extern nlmixrRemoveLikContrib_t   nlmixrRemoveLikContribP;
  extern nlmixrRegisterEmLik_t      nlmixrRegisterEmLikP;
  extern nlmixrRemoveEmLik_t        nlmixrRemoveEmLikP;
  extern nlmixrHasLikContrib_t      nlmixrHasLikContribP;
  extern nlmixrSetInnerWeightFn_t   nlmixrSetInnerWeightFnP;

  static inline SEXP iniNlmixr2estLikContrib0(SEXP p) {
    if (nlmixrRegisterLikContribP == NULL) {
      nlmixrRegisterLikContribP = (nlmixrRegisterLikContrib_t) R_ExternalPtrAddrFn(VECTOR_ELT(p, 0));
      nlmixrRemoveLikContribP   = (nlmixrRemoveLikContrib_t)   R_ExternalPtrAddrFn(VECTOR_ELT(p, 1));
      nlmixrRegisterEmLikP      = (nlmixrRegisterEmLik_t)      R_ExternalPtrAddrFn(VECTOR_ELT(p, 2));
      nlmixrRemoveEmLikP        = (nlmixrRemoveEmLik_t)        R_ExternalPtrAddrFn(VECTOR_ELT(p, 3));
      nlmixrHasLikContribP      = (nlmixrHasLikContrib_t)      R_ExternalPtrAddrFn(VECTOR_ELT(p, 4));
      // element 5 (setInnerWeightFn) may be absent when built against an older
      // nlmixr2est; guard on the table length.
      if (Rf_length(p) > 5)
        nlmixrSetInnerWeightFnP = (nlmixrSetInnerWeightFn_t) R_ExternalPtrAddrFn(VECTOR_ELT(p, 5));
    }
    return R_NilValue;
  }

#define iniNlmixr2estLikContribGlobals                          \
  nlmixrRegisterLikContrib_t nlmixrRegisterLikContribP = NULL;  \
  nlmixrRemoveLikContrib_t   nlmixrRemoveLikContribP   = NULL;  \
  nlmixrRegisterEmLik_t      nlmixrRegisterEmLikP      = NULL;  \
  nlmixrRemoveEmLik_t        nlmixrRemoveEmLikP        = NULL;  \
  nlmixrHasLikContrib_t      nlmixrHasLikContribP      = NULL;  \
  nlmixrSetInnerWeightFn_t   nlmixrSetInnerWeightFnP   = NULL;  \
  SEXP iniNlmixr2estLikContrib(SEXP p) { return iniNlmixr2estLikContrib0(p); }

#ifdef __cplusplus
}
#endif
#endif
