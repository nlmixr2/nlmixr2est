// [[Rcpp::plugins(openmp)]]
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "inner.h"
#include "odeSwap.h"
#include <algorithm>

#define _(String) (String)

// Registry of the peer solvers sharing one rxode2 solve pool.  Written only
// during (serial) setup; read during parallel solve regions, so nothing here
// mutates after a fit starts.
struct OdeModelReg {
  rxSolveF *fns = NULL;
  const char *name = NULL;
  std::vector<std::string> lhsNames;
  int neq = 0;
  int nlhs = 0;
  bool loaded = false;
};

static OdeModelReg _odeReg[odeSlotN];
static OdePoolPlan _odePlan;
static bool _odePlanStale = true;

// One preserved VECSXP holding each registered model, so the compiled rxode2
// objects cannot be collected while their entry points are live.  Mirrors
// storeCovSolveArgs_/releaseCovSolveArgs_ in inner.cpp.
static SEXP _odeModels = R_NilValue;

static void odeSwapModelsInit() {
  if (_odeModels != R_NilValue) return;
  SEXP v = PROTECT(Rf_allocVector(VECSXP, odeSlotN));
  R_PreserveObject(v);
  UNPROTECT(1);
  _odeModels = v;
}

static inline bool odeSlotOk(int slot) { return slot >= 0 && slot < odeSlotN; }

bool odeSwapDeclare(int slot, const char *name, SEXP obj) {
  if (!odeSlotOk(slot)) return false;
  odeSwapClear(slot);
  if (Rf_isNull(obj) || !rxode2::rxIs(RObject(obj), "rxode2")) return false;
  // Metadata only -- no rxDynLoad, no rxUpdateFuns.  See odeSwap.h.
  List mv = rxode2::rxModelVars_(RObject(obj));
  CharacterVector lhs = as<CharacterVector>(mv["lhs"]);
  OdeModelReg &m = _odeReg[slot];
  m.fns = NULL;
  m.name = name;
  m.neq = as<CharacterVector>(mv["state"]).size();
  m.nlhs = lhs.size();
  m.lhsNames.resize((size_t)lhs.size());
  for (int i = 0; i < lhs.size(); ++i) m.lhsNames[(size_t)i] = as<std::string>(lhs[i]);
  m.loaded = true;
  odeSwapModelsInit();
  SET_VECTOR_ELT(_odeModels, slot, obj);
  _odePlanStale = true;
  return true;
}

bool odeSwapRegister(int slot, const char *name, SEXP obj, rxSolveF *fns) {
  if (fns == NULL) return false;
  if (!odeSwapDeclare(slot, name, obj)) return false;
  // The entry points come from R_GetCCallable on the model's DLL, so it must be
  // loaded first.  Only do this for models we are about to solve/read.
  if (!rxode2::rxDynLoad(RObject(obj))) { odeSwapClear(slot); return false; }
  List mv = rxode2::rxModelVars_(RObject(obj));
  rxUpdateFuns(as<SEXP>(mv["trans"]), fns);
  _odeReg[slot].fns = fns;
  return true;
}

void odeSwapClear(int slot) {
  if (!odeSlotOk(slot)) return;
  OdeModelReg &m = _odeReg[slot];
  if (m.fns != NULL) rxClearFuns(m.fns);
  m.fns = NULL; m.name = NULL; m.neq = 0; m.nlhs = 0; m.loaded = false;
  m.lhsNames.clear();
  if (_odeModels != R_NilValue) SET_VECTOR_ELT(_odeModels, slot, R_NilValue);
  _odePlanStale = true;
}

void odeSwapClearAll() {
  for (int s = 0; s < odeSlotN; ++s) odeSwapClear(s);
  if (_odeModels != R_NilValue) {
    R_ReleaseObject(_odeModels);
    _odeModels = R_NilValue;
  }
  _odePlan = OdePoolPlan();
  _odePlanStale = true;
}

bool odeSwapLoaded(int slot) { return odeSlotOk(slot) && _odeReg[slot].loaded; }
int  odeSwapNeq(int slot)    { return odeSwapLoaded(slot) ? _odeReg[slot].neq  : 0; }
int  odeSwapNlhs(int slot)   { return odeSwapLoaded(slot) ? _odeReg[slot].nlhs : 0; }
const char *odeSwapName(int slot) { return odeSwapLoaded(slot) ? _odeReg[slot].name : NULL; }

SEXP odeSwapModelSEXP(int slot) {
  if (!odeSwapLoaded(slot) || _odeModels == R_NilValue) return R_NilValue;
  return VECTOR_ELT(_odeModels, slot);
}

int odeSwapLhsIndex(int slot, const char *nm) {
  if (!odeSwapLoaded(slot) || nm == NULL) return -1;
  const std::vector<std::string> &v = _odeReg[slot].lhsNames;
  for (size_t i = 0; i < v.size(); ++i) if (v[i] == nm) return (int)i;
  return -1;
}

// Pure: pick the pool model and decide whether a private lhs buffer is needed.
// Tie-break is max neq, then max nlhs (minimizing the scratch), then lowest slot
// -- deterministic, unlike the source-order "last writer wins" it replaces.
//
// A slot counts as loaded iff neq > 0 OR nlhs > 0.  neq alone is not enough: a
// solved-form (linCmt) or purely algebraic model has zero ODE states but real
// lhs outputs, and skipping it would under-report maxNlhs and silently drop the
// scratch buffer it needs.
OdePoolPlan odeSwapPlanFor(const std::vector<int> &neq, const std::vector<int> &nlhs) {
  OdePoolPlan p;
  size_t n = std::min(neq.size(), nlhs.size());
  for (size_t i = 0; i < n; ++i) {
    if (neq[i] <= 0 && nlhs[i] <= 0) continue;
    p.nLoaded++;
    if (p.poolSlot < 0 || neq[i] > p.poolNeq ||
        (neq[i] == p.poolNeq && nlhs[i] > p.poolNlhs)) {
      p.poolSlot = (int)i; p.poolNeq = neq[i]; p.poolNlhs = nlhs[i];
    }
    if (p.maxNlhsSlot < 0 || nlhs[i] > p.maxNlhs) {
      p.maxNlhs = nlhs[i]; p.maxNlhsSlot = (int)i;
    }
  }
  if (p.poolSlot < 0) return p;
  for (size_t i = 0; i < n; ++i) {
    // same loaded predicate as above; a loaded 0-state model still needs the
    // override to compact the stride (rxEffNeq accepts 0)
    if ((neq[i] > 0 || nlhs[i] > 0) && neq[i] < p.poolNeq) { p.overrideNeeded = true; break; }
  }
  // The widest lhs does not belong to the pool model: rxode2's per-thread slice
  // is only poolNlhs wide, so reads of that model need our own buffer.
  if (p.maxNlhs > p.poolNlhs) p.scratchNlhs = p.maxNlhs;
  return p;
}

const OdePoolPlan &odeSwapPlan() {
  if (_odePlanStale) {
    std::vector<int> neq((size_t)odeSlotN, 0), nlhs((size_t)odeSlotN, 0);
    for (int s = 0; s < odeSlotN; ++s) {
      if (!_odeReg[s].loaded) continue;
      neq[(size_t)s] = _odeReg[s].neq;
      nlhs[(size_t)s] = _odeReg[s].nlhs;
    }
    _odePlan = odeSwapPlanFor(neq, nlhs);
    _odePlanStale = false;
  }
  return _odePlan;
}

SEXP odeSwapPoolModelSEXP() {
  const OdePoolPlan &p = odeSwapPlan();
  return (p.poolSlot < 0) ? R_NilValue : odeSwapModelSEXP(p.poolSlot);
}

int odeSwapCanPool(int slot) {
  if (!odeSwapLoaded(slot) || _odeReg[slot].fns == NULL ||
      _odeReg[slot].fns->calc_lhs == NULL) return odeDenyNotLoaded;
  if (_odeReg[slot].nlhs <= 0) return odeDenyNlhsUnknown;
  const OdePoolPlan &p = odeSwapPlan();
  // rxEffNeq only ever compacts: an override above op->neq silently falls back,
  // so a model larger than the pool cannot be solved in it.
  if (p.poolSlot < 0 || _odeReg[slot].neq > p.poolNeq) return odeDenyPoolNotSized;
  return odeDenyNone;
}

// ---- R introspection ----------------------------------------------------
//
// Exposes the pool decision so tests can assert the MECHANISM (which model sized
// the pool, whether the scratch was needed, whether an override leaked), not just
// that the numbers came out the same.

//[[Rcpp::export]]
List odeSwapPlanFor_(IntegerVector neq, IntegerVector nlhs) {
  std::vector<int> a(neq.begin(), neq.end()), b(nlhs.begin(), nlhs.end());
  OdePoolPlan p = odeSwapPlanFor(a, b);
  return List::create(_["poolSlot"] = p.poolSlot, _["poolNeq"] = p.poolNeq,
                      _["poolNlhs"] = p.poolNlhs, _["maxNlhs"] = p.maxNlhs,
                      _["maxNlhsSlot"] = p.maxNlhsSlot,
                      _["scratchNlhs"] = p.scratchNlhs,
                      _["needsScratch"] = (p.scratchNlhs > 0),
                      _["nLoaded"] = p.nLoaded,
                      _["overrideNeeded"] = p.overrideNeeded);
}

//[[Rcpp::export]]
List odeSwapInfo_() {
  const OdePoolPlan &p = odeSwapPlan();
  CharacterVector nm(odeSlotN);
  IntegerVector neq(odeSlotN), nlhs(odeSlotN), deny(odeSlotN);
  LogicalVector loaded(odeSlotN), sizesPool(odeSlotN);
  for (int s = 0; s < odeSlotN; ++s) {
    nm[s] = _odeReg[s].name == NULL ? NA_STRING : Rf_mkChar(_odeReg[s].name);
    neq[s] = _odeReg[s].neq;
    nlhs[s] = _odeReg[s].nlhs;
    loaded[s] = _odeReg[s].loaded;
    sizesPool[s] = (s == p.poolSlot);
    deny[s] = odeSwapCanPool(s);
  }
  List models = List::create(_["slot"] = seq_len(odeSlotN) - 1, _["name"] = nm,
                             _["neq"] = neq, _["nlhs"] = nlhs,
                             _["loaded"] = loaded, _["sizesPool"] = sizesPool,
                             _["deny"] = deny);
  models.attr("class") = "data.frame";
  models.attr("row.names") = IntegerVector::create(NA_INTEGER, -odeSlotN);
  // op->neq / op->nlhs are only meaningful once a solve pool exists.
  int opNeq = NA_INTEGER, opNlhs = NA_INTEGER;
  rx_solve *rxl = getRxSolve_();
  if (rxl != NULL) {
    rx_solving_options *op = getSolvingOptions(rxl);
    if (op != NULL) { opNeq = getOpNeq(op); opNlhs = getOpNlhs(op); }
  }
  return List::create(
    _["models"] = models,
    _["poolSlot"] = p.poolSlot,
    _["poolName"] = (p.poolSlot < 0 || _odeReg[p.poolSlot].name == NULL) ?
      CharacterVector::create(NA_STRING) : CharacterVector::create(_odeReg[p.poolSlot].name),
    _["poolNeq"] = p.poolNeq, _["poolNlhs"] = p.poolNlhs,
    _["maxNlhs"] = p.maxNlhs, _["maxNlhsSlot"] = p.maxNlhsSlot,
    _["scratchNlhs"] = p.scratchNlhs,
    _["needsScratch"] = (p.scratchNlhs > 0),
    _["overrideNeeded"] = p.overrideNeeded,
    _["nLoaded"] = p.nLoaded,
    _["opNeq"] = opNeq, _["opNlhs"] = opNlhs);
}
