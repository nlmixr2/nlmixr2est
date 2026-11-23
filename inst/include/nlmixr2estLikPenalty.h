#ifndef NLMIXR2EST_LIKPENALTY_H
#define NLMIXR2EST_LIKPENALTY_H
// Likelihood penalty ABI: per-subject and population terms added to the
// objective of every nlmixr2est estimation method.
//
// Scales: a SUBJECT term adds log-likelihood (LL); a POPULATION term returns a
// penalty on the -2LL (objective) scale.  So a subject term of c on every one of
// N subjects moves the objective by -2*N*c, and a population term P moves it by P.
//
// Threads: `subject` runs inside the OpenMP subject loop -- no R API, no
// allocation, key any state by `id`.  `population` runs on the main thread.
//
// See vignette("likelihood-penalties", package = "nlmixr2est").

#ifdef __cplusplus
extern "C" {
#endif

#define NLMIXR_PEN_ABI 1

// flags
#define NLMIXR_PEN_IN_OBJF 1  // count the penalty in the reported objective
#define NLMIXR_PEN_NO_GRAD 2  // `population` cannot fill dObj_dtheta

  // Per-penalty context: what the spec selected (nlmixr2LikPenaltySpec()).
  typedef struct {
    int nsel;            // number of selected thetas
    const int *sel;      // [nsel] 0-based theta indices (ui$iniDf theta order)
    int npar;            // number of numeric settings
    const double *par;   // [npar] numeric settings (e.g. a lambda)
    void *data;          // opaque pointer from the spec / bundle
  } nlmixrPenCtx;

  typedef struct {
    int id;               // 0-based subject id
    int neta;             // etas (0 for population-only methods like nlm)
    int nobs;             // observations of this subject
    int ntheta;           // 0 when the method has no theta vector here
    const double *eta;    // [neta] current etas (NULL when neta == 0)
    const double *theta;  // [ntheta] natural-scale thetas (may be NULL)
    double *llik;         // += extra log-likelihood for the subject
    double *dLL_deta;     // [neta] += extra d(LL)/d(eta) (NULL when neta == 0)
  } nlmixrPenSubj;

  typedef struct {
    int ntheta;            // thetas
    int neta;              // omega dimension (0 when unavailable)
    int nsub;              // subjects
    const double *theta;   // [ntheta] natural-scale thetas
    const double *omega;   // [neta*neta] current omega, ini() eta order; every
                           // method with etas passes it (NULL when neta == 0)
    double *dObj_dtheta;   // [ntheta] += d(penalty)/d(theta); NULL if not wanted
                           // (an omega dependence is finite-differenced)
  } nlmixrPenPop;

  typedef struct {
    int abi;                // NLMIXR_PEN_ABI
    int flags;              // NLMIXR_PEN_* defaults (a spec may override IN_OBJF)
    const char *name;       // short label used in output
    void (*subject)(nlmixrPenSubj *s, const nlmixrPenCtx *ctx);       // optional
    double (*population)(nlmixrPenPop *p, const nlmixrPenCtx *ctx);  // optional
    void *data;             // default ctx->data when registered globally
  } nlmixrLikPenalty;

  // Session-wide registration (every fit).  Prefer per-fit activation through
  // nlmixr2LikPenaltySpec() on the model's ui.  Idempotent per pointer.
  void nlmixrRegisterLikPenalty(const nlmixrLikPenalty *p);
  void nlmixrRemoveLikPenalty(const nlmixrLikPenalty *p);

#ifdef __cplusplus
}
#endif
#endif
