#ifndef NLMIXR2EST_LIKCONTRIB_H
#define NLMIXR2EST_LIKCONTRIB_H
// External likelihood-contribution registry.
//
// A contributor package (e.g. nlmixr2nn) registers a bundle of plain C function
// pointers of the fixed signatures below; nlmixr2est cycles all registered
// bundles IN SERIES inside its per-subject objective.  No R API is used: the
// hooks fire from inside the OpenMP parallel region, so contributors must key
// any per-subject state by `id` and use no R calls.
//
// FOCEI family: begin/obs/end fire inside likInner0's per-observation loop (so
// lpInner, which copies fInd->lp, inherits any gradient contribution).  The base
// solve outputs and cotangents are READ-ONLY; the contributor MAY add an extra
// log-likelihood term and extra d(LL)/d(eta) (both optional -- leave untouched
// to only observe, e.g. to record d(LL)/d(f) for an adjoint sweep).

#ifdef __cplusplus
extern "C" {
#endif

  // per-subject bracket
  typedef struct {
    int id;             // internal subject id (0-based)
    int neta;           // number of etas
    int nobs;           // observations for this subject
    const double *eta;  // [neta] current etas
  } nlmixrLikSubj;

  // per-observation context (one call per normal/non-normal observation)
  typedef struct {
    int id;                 // subject id
    int k;                  // observation index within the subject
    int neta;               // number of etas
    double f;               // prediction
    double dv;              // observation
    double r;               // residual variance
    double dLL_df;          // d(LL)/d(f)   cotangent (the essential read)
    double dLL_dr;          // d(LL)/d(r)   cotangent
    const double *df_deta;  // [neta] base d(f)/d(eta) (NULL when neta == 0)
    double *llik;           // += optional extra LL at this observation
    double *dLL_deta;       // [neta] += optional extra d(LL)/d(eta)
  } nlmixrLikObs;

  typedef struct {
    void (*beginSubject)(const nlmixrLikSubj *s);  // optional (NULL ok)
    void (*obs)(nlmixrLikObs *o);                  // required
    void (*endSubject)(const nlmixrLikSubj *s);    // optional (NULL ok)
  } nlmixrLikContrib;

  // EM / SAEM family: additive per-subject log-likelihood estimate.
  typedef struct { int id; int neta; int nobs; const double *eta; } nlmixrEmSubj;
  typedef double (*nlmixrEmLik_fn)(const nlmixrEmSubj *s);

  // Per-subject NN weight-injection hook: called right after a subject's etas
  // are written to par_ptr and before its inner solve, so a plugin can overwrite
  // that subject's par_ptr weight block with individual weights W_i computed from
  // the eta vector (e.g. W_i = lW * exp(etaW)).  It fires on EVERY eta-set,
  // including finite-difference perturbations, so FOCEI's FD eta-sensitivity
  // captures d(f)/d(etaW) through the injected weights with no analytic
  // gradient -- the inner (individual-weight) building block.  Thread-safe: it
  // may run in the parallel per-subject region, so it must touch only subject
  // cid's par_ptr.
  typedef void (*nlmixrInnerWeight_fn)(int cid, const double *eta, int neta);

  // registry (exposed to contributor packages via R_GetCCallable-free means;
  // see nlmixr2est init).  Register/remove is idempotent per pointer.
  void nlmixrRegisterLikContrib(const nlmixrLikContrib *c);
  void nlmixrRemoveLikContrib(const nlmixrLikContrib *c);
  void nlmixrRegisterEmLik(nlmixrEmLik_fn fn);
  void nlmixrRemoveEmLik(nlmixrEmLik_fn fn);
  int  nlmixrHasLikContrib(void);   // fast check to skip hook overhead
  void nlmixrSetInnerWeightFn(nlmixrInnerWeight_fn fn);  // NULL to clear

#ifdef __cplusplus
}
#endif
#endif
