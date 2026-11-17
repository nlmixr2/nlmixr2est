// Likelihood penalty registry: session-wide bundles (nlmixrRegisterLikPenalty)
// plus per-fit specs activated from R (.likPenaltyActivate).  Both are flattened
// into one read-only table before a fit so the hot path is a plain loop.
#define STRICT_R_HEADERS
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <atomic>
#include <vector>
#include <string>
#include <cstring>
#include "likPenalty.h"

#define NLMIXR_MAX_PEN 16

struct nlmixrPenEntry {
  const nlmixrLikPenalty *pen;
  int inObjf;
  std::vector<int> sel;
  std::vector<double> par;
  void *data;
  std::string name;
  nlmixrPenCtx ctx;
};

static const nlmixrLikPenalty *_penGlobal[NLMIXR_MAX_PEN] = {NULL};
static int _penNGlobal = 0;
static std::vector<nlmixrPenEntry> _penFit;     // per-fit specs
static std::vector<nlmixrPenEntry> _penTable;   // flattened global + per-fit
static int _penNSubj = 0, _penNPop = 0;
static std::atomic<int> _penSubjChanged{0};

static void penBindCtx(nlmixrPenEntry &e) {
  e.ctx.nsel = (int)e.sel.size();
  e.ctx.sel = e.sel.empty() ? NULL : e.sel.data();
  e.ctx.npar = (int)e.par.size();
  e.ctx.par = e.par.empty() ? NULL : e.par.data();
  e.ctx.data = e.data;
}

static void penRebuild(void) {
  _penTable.clear();
  _penTable.reserve(_penNGlobal + _penFit.size());
  for (int i = 0; i < _penNGlobal; ++i) {
    nlmixrPenEntry e;
    e.pen = _penGlobal[i];
    e.inObjf = (e.pen->flags & NLMIXR_PEN_IN_OBJF) ? 1 : 0;
    e.data = e.pen->data;
    e.name = e.pen->name == NULL ? "" : e.pen->name;
    _penTable.push_back(e);
  }
  for (size_t i = 0; i < _penFit.size(); ++i) _penTable.push_back(_penFit[i]);
  // ctx pointers into the vectors are only valid once the table stops growing
  _penNSubj = _penNPop = 0;
  for (size_t i = 0; i < _penTable.size(); ++i) {
    penBindCtx(_penTable[i]);
    if (_penTable[i].pen->subject != NULL) _penNSubj++;
    if (_penTable[i].pen->population != NULL) _penNPop++;
  }
  _penSubjChanged.store(0, std::memory_order_relaxed);
}

extern "C" void nlmixrRegisterLikPenalty(const nlmixrLikPenalty *p) {
  if (p == NULL || p->abi != NLMIXR_PEN_ABI) return;
  if (p->subject == NULL && p->population == NULL) return;
  for (int i = 0; i < _penNGlobal; ++i) if (_penGlobal[i] == p) return;
  if (_penNGlobal < NLMIXR_MAX_PEN) _penGlobal[_penNGlobal++] = p;
  penRebuild();
}

extern "C" void nlmixrRemoveLikPenalty(const nlmixrLikPenalty *p) {
  for (int i = 0; i < _penNGlobal; ++i) if (_penGlobal[i] == p) {
    for (int k = i; k < _penNGlobal - 1; ++k) _penGlobal[k] = _penGlobal[k + 1];
    _penGlobal[--_penNGlobal] = NULL;
    penRebuild();
    return;
  }
}

extern "C" int nlmixrPenNSubject(void) { return _penNSubj; }
extern "C" int nlmixrPenNPop(void) { return _penNPop; }
extern "C" int nlmixrPenSubjChanged(void) {
  return _penSubjChanged.load(std::memory_order_relaxed);
}

static inline double penSubjOne(const nlmixrPenEntry &e, int id, int neta, int nobs,
                                const double *eta, int ntheta, const double *theta,
                                double *dLLdEta) {
  double ll = 0.0;
  nlmixrPenSubj s;
  s.id = id; s.neta = neta; s.nobs = nobs; s.ntheta = ntheta;
  s.eta = neta > 0 ? eta : NULL; s.theta = ntheta > 0 ? theta : NULL;
  s.llik = &ll; s.dLL_deta = (neta > 0) ? dLLdEta : NULL;
  e.pen->subject(&s, &e.ctx);
  return ll;
}

extern "C" double nlmixrPenSubjAll(int id, int neta, int nobs, const double *eta,
                                   int ntheta, const double *theta, double *dLLdEta) {
  if (_penNSubj == 0) return 0.0;
  double ll = 0.0;
  double dSave[64];
  bool track = (dLLdEta != NULL && neta > 0 && neta <= 64);
  if (track) std::memcpy(dSave, dLLdEta, sizeof(double) * neta);
  for (size_t i = 0; i < _penTable.size(); ++i) {
    const nlmixrPenEntry &e = _penTable[i];
    if (e.pen->subject == NULL) continue;
    ll += penSubjOne(e, id, neta, nobs, eta, ntheta, theta, dLLdEta);
  }
  bool changed = (ll != 0.0);
  if (!changed && track) {
    for (int q = 0; q < neta; ++q) if (dLLdEta[q] != dSave[q]) { changed = true; break; }
  }
  if (changed) _penSubjChanged.store(1, std::memory_order_relaxed);
  return ll;
}

static inline double penPopOne(const nlmixrPenEntry &e, int ntheta, const double *theta,
                               int neta, const double *omega, int nsub, double *grad) {
  nlmixrPenPop p;
  p.ntheta = ntheta; p.neta = neta; p.nsub = nsub;
  p.theta = theta; p.omega = neta > 0 ? omega : NULL;
  p.dObj_dtheta = grad;
  return e.pen->population(&p, &e.ctx);
}

extern "C" double nlmixrPenPopAll(int ntheta, const double *theta, int neta,
                                  const double *omega, int nsub, double *grad, int *gradOk) {
  if (gradOk != NULL) *gradOk = 1;
  if (_penNPop == 0) return 0.0;
  double ret = 0.0;
  for (size_t i = 0; i < _penTable.size(); ++i) {
    const nlmixrPenEntry &e = _penTable[i];
    if (e.pen->population == NULL) continue;
    if (grad != NULL && (e.pen->flags & NLMIXR_PEN_NO_GRAD) && gradOk != NULL) *gradOk = 0;
    ret += penPopOne(e, ntheta, theta, neta, omega, nsub, grad);
  }
  return ret;
}

// --- built-in penalties -------------------------------------------------------

// thetaL2: lambda * sum((theta[sel] - center)^2) + offset on the -2LL scale.
// par = c(lambda, center, offset) (center and offset default to 0).
static double penThetaL2Pop(nlmixrPenPop *p, const nlmixrPenCtx *ctx) {
  double lambda = ctx->npar > 0 ? ctx->par[0] : 0.0;
  double center = ctx->npar > 1 ? ctx->par[1] : 0.0;
  double ret = ctx->npar > 2 ? ctx->par[2] : 0.0;
  for (int j = 0; j < ctx->nsel; ++j) {
    int t = ctx->sel[j];
    if (t < 0 || t >= p->ntheta) continue;
    double d = p->theta[t] - center;
    ret += lambda * d * d;
    if (p->dObj_dtheta != NULL) p->dObj_dtheta[t] += 2.0 * lambda * d;
  }
  return ret;
}
static const nlmixrLikPenalty _penThetaL2 = {
  NLMIXR_PEN_ABI, 0, "thetaL2", NULL, penThetaL2Pop, NULL};

// etaL2: subject LL += c - 0.5 * lambda * sum(eta^2).  par = c(lambda, c).
// Mostly a test/teaching penalty: c is a pure objective shift (-2*N*c).
static void penEtaL2Subj(nlmixrPenSubj *s, const nlmixrPenCtx *ctx) {
  double lambda = ctx->npar > 0 ? ctx->par[0] : 0.0;
  double c = ctx->npar > 1 ? ctx->par[1] : 0.0;
  double ll = c;
  for (int q = 0; q < s->neta; ++q) {
    ll -= 0.5 * lambda * s->eta[q] * s->eta[q];
    if (s->dLL_deta != NULL) s->dLL_deta[q] -= lambda * s->eta[q];
  }
  *s->llik += ll;
}
static const nlmixrLikPenalty _penEtaL2 = {
  NLMIXR_PEN_ABI, 0, "etaL2", penEtaL2Subj, NULL, NULL};

static void penPtrFinalizer(SEXP p) { (void)p; }

extern "C" SEXP _nlmixr2est_likPenaltyBuiltin(SEXP nameSxp) {
  const char *nm = CHAR(STRING_ELT(nameSxp, 0));
  const nlmixrLikPenalty *p = NULL;
  if (!strcmp(nm, "thetaL2")) p = &_penThetaL2;
  else if (!strcmp(nm, "etaL2")) p = &_penEtaL2;
  else Rf_error("unknown built-in likelihood penalty '%s'", nm);
  SEXP ret = PROTECT(R_MakeExternalPtr((void *)p, Rf_install("nlmixrLikPenalty"), R_NilValue));
  R_RegisterCFinalizer(ret, penPtrFinalizer);
  UNPROTECT(1);
  return ret;
}

// --- R interface ----------------------------------------------------------------

// ptrs: list of EXTPTRSXP (nlmixrLikPenalty*); sel: list of 0-based integer
// vectors; par: list of doubles; data: list of EXTPTRSXP or NULL; inObjf: int
// (NA -> the bundle's own flag).
extern "C" SEXP _nlmixr2est_likPenaltyActivate(SEXP ptrs, SEXP sel, SEXP par,
                                               SEXP data, SEXP inObjf) {
  _penFit.clear();
  int n = Rf_length(ptrs);
  for (int i = 0; i < n; ++i) {
    SEXP pp = VECTOR_ELT(ptrs, i);
    if (TYPEOF(pp) != EXTPTRSXP) Rf_error("likelihood penalty %d is not an external pointer", i + 1);
    const nlmixrLikPenalty *pen = (const nlmixrLikPenalty *)R_ExternalPtrAddr(pp);
    if (pen == NULL) Rf_error("likelihood penalty %d is a NULL pointer", i + 1);
    if (pen->abi != NLMIXR_PEN_ABI) {
      Rf_error("likelihood penalty %d has ABI %d, expected %d", i + 1, pen->abi, NLMIXR_PEN_ABI);
    }
    nlmixrPenEntry e;
    e.pen = pen;
    SEXP s = VECTOR_ELT(sel, i);
    for (int j = 0; j < Rf_length(s); ++j) e.sel.push_back(INTEGER(s)[j]);
    SEXP pr = VECTOR_ELT(par, i);
    for (int j = 0; j < Rf_length(pr); ++j) e.par.push_back(REAL(pr)[j]);
    SEXP d = VECTOR_ELT(data, i);
    e.data = (TYPEOF(d) == EXTPTRSXP) ? R_ExternalPtrAddr(d) : pen->data;
    int io = INTEGER(inObjf)[i];
    e.inObjf = (io == NA_INTEGER) ? ((pen->flags & NLMIXR_PEN_IN_OBJF) ? 1 : 0) : io;
    e.name = pen->name == NULL ? "" : pen->name;
    _penFit.push_back(e);
  }
  penRebuild();
  return Rf_ScalarInteger((int)_penTable.size());
}

extern "C" SEXP _nlmixr2est_likPenaltyDeactivate(void) {
  _penFit.clear();
  penRebuild();
  return R_NilValue;
}

// Evaluate every active penalty at a point (post-fit reporting).
// theta: natural-scale thetas; eta: nsub x neta matrix (or 0-col); omega: neta x
// neta (or NULL); nobs: per-subject observation counts.  Returns a list with the
// per-penalty name, inObjf, population value and summed subject LL.
extern "C" SEXP _nlmixr2est_likPenaltyEval(SEXP thetaSxp, SEXP etaSxp, SEXP omegaSxp,
                                           SEXP nobsSxp) {
  int ntheta = Rf_length(thetaSxp);
  const double *theta = ntheta > 0 ? REAL(thetaSxp) : NULL;
  int nsub = Rf_length(nobsSxp);
  int neta = 0;
  if (Rf_isMatrix(etaSxp)) neta = Rf_ncols(etaSxp);
  int nom = Rf_isMatrix(omegaSxp) ? Rf_nrows(omegaSxp) : 0;
  const double *omega = nom > 0 ? REAL(omegaSxp) : NULL;
  int n = (int)_penTable.size();
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP io = PROTECT(Rf_allocVector(LGLSXP, n));
  SEXP pop = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP sub = PROTECT(Rf_allocVector(REALSXP, n));
  std::vector<double> eta(neta > 0 ? neta : 1), dEta(neta > 0 ? neta : 1);
  for (int i = 0; i < n; ++i) {
    const nlmixrPenEntry &e = _penTable[i];
    SET_STRING_ELT(nm, i, Rf_mkChar(e.name.c_str()));
    LOGICAL(io)[i] = e.inObjf;
    REAL(pop)[i] = e.pen->population == NULL ? 0.0 :
      penPopOne(e, ntheta, theta, nom, omega, nsub, NULL);
    double ll = 0.0;
    if (e.pen->subject != NULL) {
      for (int id = 0; id < nsub; ++id) {
        for (int q = 0; q < neta; ++q) {
          eta[q] = REAL(etaSxp)[id + q * nsub];
          dEta[q] = 0.0;
        }
        ll += penSubjOne(e, id, neta, INTEGER(nobsSxp)[id], eta.data(), ntheta, theta,
                         dEta.data());
      }
    }
    REAL(sub)[i] = ll;
  }
  SEXP ret = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ret, 0, nm); SET_VECTOR_ELT(ret, 1, io);
  SET_VECTOR_ELT(ret, 2, pop); SET_VECTOR_ELT(ret, 3, sub);
  SEXP rn = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(rn, 0, Rf_mkChar("name")); SET_STRING_ELT(rn, 1, Rf_mkChar("inObjf"));
  SET_STRING_ELT(rn, 2, Rf_mkChar("population")); SET_STRING_ELT(rn, 3, Rf_mkChar("subjectLL"));
  Rf_setAttrib(ret, R_NamesSymbol, rn);
  UNPROTECT(6);
  return ret;
}

extern "C" SEXP _nlmixr2est_likPenaltyInfo(void) {
  SEXP ret = PROTECT(Rf_allocVector(INTSXP, 4));
  INTEGER(ret)[0] = (int)_penTable.size();
  INTEGER(ret)[1] = _penNSubj;
  INTEGER(ret)[2] = _penNPop;
  INTEGER(ret)[3] = _penSubjChanged.load(std::memory_order_relaxed);
  UNPROTECT(1);
  return ret;
}

// Per-subject total subject-penalty LL at the given etas (nsub x neta).
extern "C" SEXP _nlmixr2est_likPenaltySubjLL(SEXP thetaSxp, SEXP etaSxp, SEXP nobsSxp) {
  int ntheta = Rf_length(thetaSxp);
  const double *theta = ntheta > 0 ? REAL(thetaSxp) : NULL;
  int nsub = Rf_length(nobsSxp);
  int neta = Rf_isMatrix(etaSxp) ? Rf_ncols(etaSxp) : 0;
  SEXP ret = PROTECT(Rf_allocVector(REALSXP, nsub));
  std::vector<double> eta(neta > 0 ? neta : 1), dEta(neta > 0 ? neta : 1);
  for (int id = 0; id < nsub; ++id) {
    for (int q = 0; q < neta; ++q) eta[q] = REAL(etaSxp)[id + q * nsub];
    REAL(ret)[id] = nlmixrPenSubjAll(id, neta, INTEGER(nobsSxp)[id], eta.data(), ntheta,
                                     theta, NULL);
  }
  UNPROTECT(1);
  return ret;
}

// Population penalty at theta (no omega), -2LL scale.
extern "C" SEXP _nlmixr2est_likPenaltyPop(SEXP thetaSxp, SEXP nsubSxp) {
  int ntheta = Rf_length(thetaSxp);
  return Rf_ScalarReal(nlmixrPenPopAll(ntheta, ntheta > 0 ? REAL(thetaSxp) : NULL, 0, NULL,
                                       Rf_asInteger(nsubSxp), NULL, NULL));
}
