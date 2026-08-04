// [[Rcpp::plugins(openmp)]]
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "inner.h"
#include <set>
#include "odeSwap.h"
#include <algorithm>
#include <atomic>

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
  // Endpoint (CMT) rebasing -- see odeSwapCmtRebase().  rxode2 compiles each
  // model's endpoint switch against the USER compartment numbering and emits
  //   #define _CMT ((fabs(CMT)<=nPhys) ? CMT : CMT - nSens)
  // so a model normalizes the raw CMT with its OWN sensitivity-compartment
  // count.  Peers sharing one pooled event table therefore need the raw value
  // rebased by the difference; nSens and the CMT parameter slot are recorded
  // here so that can be done without going back to R.
  int nSens = 0;      // length(rxModelVars(obj)$sens)
  int cmtPar = -1;    // index of "CMT" in $params, or -1 if the model has none
  // event-sensitivity shape (see OdeSwapEsBatch); esActive == 0 for a model
  // with no jump sensitivities, e.g. rxPred
  int esActive = 0;   // does this model carry event ("jump") sensitivities?
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
  // Matches rxode2's op->neq, which counts "state" only -- NOT stateExtra
  // (measured: state=8, stateExtra=1, op->neq=8).  Adding stateExtra pushes neq
  // above op->neq, and rxEffNeq() then ignores the override entirely, silently
  // disabling every compaction.
  m.neq = as<CharacterVector>(mv["state"]).size();
  m.nlhs = lhs.size();
  m.lhsNames.resize((size_t)lhs.size());
  for (int i = 0; i < lhs.size(); ++i) m.lhsNames[(size_t)i] = as<std::string>(lhs[i]);
  // CMT rebasing inputs (see the struct comment and odeSwapCmtRebase)
  m.nSens = mv.containsElementNamed("sens") ?
    as<CharacterVector>(mv["sens"]).size() : 0;
  m.cmtPar = -1;
  if (mv.containsElementNamed("params")) {
    CharacterVector pars = as<CharacterVector>(mv["params"]);
    for (int i = 0; i < pars.size(); ++i) {
      if (as<std::string>(pars[i]) == "CMT") { m.cmtPar = i; break; }
    }
  }
  // Record only WHETHER this model carries event ("jump") sensitivities.  The
  // shape itself is installed through rxode2's R entry point, which derives the
  // dims from the model, so the registry does not need to duplicate them --
  // and must not try: eventSensInfo$map's fields are not all integer
  // (map2$q is character), so parsing them here throws.
  m.esActive = 0;
  try {
    RObject _o(obj);
    SEXP _v = R_NilValue;
    if (_o.hasAttribute("eventSensInfo")) {
      _v = _o.attr("eventSensInfo");
    } else if (TYPEOF(obj) == ENVSXP &&
               R_existsVarInFrame(obj, Rf_install("eventSensInfo"))) {
      // check the binding first: Rf_eval on a missing symbol raises an R error,
      // i.e. a longjmp past C++ destructors that try/catch cannot catch
      _v = Rf_eval(Rf_install("eventSensInfo"), obj);
    }
    if (!Rf_isNull(_v) && TYPEOF(_v) == VECSXP) {
      List _e(_v);
      bool _fd = _e.containsElementNamed("mode") &&
        as<std::string>(_e["mode"]) == "fd";
      if (!_fd && _e.containsElementNamed("map")) m.esActive = 1;
    }
  } catch (...) {
    m.esActive = 0;
  }
  m.loaded = true;
  odeSwapModelsInit();
  SET_VECTOR_ELT(_odeModels, slot, obj);
  _odePlanStale = true;
  return true;
}

bool odeSwapRegister(int slot, const char *name, SEXP obj, rxSolveF *fns) {
  if (fns == NULL) return false;
  if (!odeSwapDeclare(slot, name, obj)) return false;
  // Deliberately NO rxDynLoad here.  The code this replaced called only
  // rxUpdateFuns for each peer model, and adding a load broke the modeled-dosing
  // (f/lag) analytic gradient: rxDynLoad rebinds rxode2's event-sensitivity
  // globals, and the f/lag jump sensitivities ARE event sensitivities, so
  // re-loading a peer mid-setup corrupts them.  R hands these models in already
  // loaded (foceiSetup_ rxDynLoads the inner model itself); if one were not,
  // rxUpdateFuns would fail loudly rather than silently compute a wrong result.
  List mv = rxode2::rxModelVars_(RObject(obj));
  rxUpdateFuns(as<SEXP>(mv["trans"]), fns);
  _odeReg[slot].fns = fns;
  return true;
}

bool odeSwapHasEs(int slot) {
  return odeSlotOk(slot) && _odeReg[slot].loaded && _odeReg[slot].esActive != 0;
}

// Which ROLE's ES shape is installed (odeEsUnknown = we have not recorded one),
// and which SLOT installed it.  Both are needed and they are NOT interchangeable:
// the ROLE decides whether a solve may compact, while restoring has to re-install
// a concrete slot's model.  Conflating them installed slot 1 (pred) when restoring
// role 1 (inner).
static int _odeEsSlot = odeEsUnknown;
static int _odeEsSlotIdx = -1;

int odeSwapEsModelForSlot(int slot) {
  switch (slot) {
  case odeSlotPred:      return odeEsPred;
  case odeSlotInner:     return odeEsInner;
  case odeSlotHess2:     return odeEsHess2;
  case odeSlotThetaSens:
  case odeSlotOuter:
  case odeSlotOuterNode:
  case odeSlotOuterCov:  return odeEsOuter;
  default:               return odeEsUnknown;
  }
}

int  odeSwapEsInstalledModel()          { return _odeEsSlot; }
void odeSwapEsNoteInstalled(int esModel) { _odeEsSlot = esModel; _odeEsSlotIdx = -1; }

// Install a slot's ES shape.
//
// The INSTALL still goes through R: rxode2EventSensLoadFull() wants all six dims and
// they are derived from model$eventSensInfo, an R-level structure that
// rxEventSensLoadModel() already knows how to read.  Acceptable because this is a
// batch boundary -- once per model batch, outside any parallel region.  It must never
// be called from inside an OpenMP loop (Rcpp::Function is not safe there).
//
// The RESTORE is what moved onto the new C API (rxode2 5.1.6, PR #1175): see
// OdeSwapEsBatch, which snapshots the live shape with rxode2EventSensShapeSave().
//
// Per-slot shape buffers were tried and REMOVED.  They would have taken R off this
// path entirely, but a saved shape holds pointers into the model's shared library and
// rxUnloadAll()/rxUnload() is driven from R without going through odeSwapClear(), so a
// cache entry can outlive the library it points into and the next restore would install
// dangling function pointers.  The short-lived snapshot below cannot: it is taken and
// consumed inside one batch scope.
static bool odeSwapEsInstall(int slot) {
  if (!odeSlotOk(slot) || !_odeReg[slot].loaded) return false;
  SEXP _m = odeSwapModelSEXP(slot);
  if (_m == R_NilValue) return false;
  try {
    Environment _rx = Environment::namespace_env("rxode2");
    Function _f = as<Function>(_rx["rxEventSensLoadModel"]);
    _f(RObject(_m));
  } catch (...) {
    // leave the shape as it was; the caller falls back rather than mis-solving
    return false;
  }
  return true;
}

static void odeSwapEsDeactivate() {
#ifdef NLMIXR2EST_HAS_ESSHAPE
  if (rxode2EventSensDeactivate != NULL) { rxode2EventSensDeactivate(); return; }
#endif
  {
  try {                                   // COMPATIBILITY (rxode2 5.1.5): via R
    Environment _rx = Environment::namespace_env("rxode2");
    Function _f = as<Function>(_rx["rxEventSensDeactivate"]);
    _f();
  } catch (...) {
  }
  }
}

OdeSwapEsBatch::OdeSwapEsBatch(int slot)
  : prevSlot_(_odeEsSlot), prevSlotIdx_(_odeEsSlotIdx), armed_(false) {
  // Snapshot whatever shape is live, whoever installed it.  This is what lets the
  // destructor put back a shape the registry never saw -- focei.R's fit-wide
  // rxEventSensLoadModel(model$inner) -- without having to recognize it.
  int want = odeSwapEsModelForSlot(slot);
  // Compare the SLOT, not the role: thetaSens/outer/outerNode/outerCov all share the
  // odeEsOuter role but are DIFFERENT compiled models with different ES shapes, so a
  // role match would skip installing the one we are about to solve.
  if (slot == _odeEsSlotIdx) return;                       // already this model
  saveLive();                        // BEFORE anything is installed over it
  bool esOk = true;
  if (odeSwapHasEs(slot)) {
    esOk = odeSwapEsInstall(slot);
    // Record the role ONLY when the install actually happened.  Recording a role
    // we failed to install would tell OdeSwapScope it may compact against a shape
    // that is not there.
    if (!esOk) return;
  } else {
    // A model with no jump sensitivities must not be solved under another
    // model's shape either: the injection would fire with dims that belong to
    // a different model.  Deactivate for the batch.
    odeSwapEsDeactivate();
  }
  _odeEsSlot = want;
  _odeEsSlotIdx = slot;
  armed_ = true;
}

// Capture the live shape into prevShape_ (called before anything is installed).
bool OdeSwapEsBatch::saveLive() {
#ifdef NLMIXR2EST_HAS_ESSHAPE
  if (rxode2EventSensShapeSize == NULL || rxode2EventSensShapeSave == NULL) return false;
  int sz = rxode2EventSensShapeSize();
  if (sz <= 0) return false;
  prevShape_.resize((size_t)sz);
  if (!rxode2EventSensShapeSave(prevShape_.data(), sz)) { prevShape_.clear(); return false; }
  return true;
#else
  // rxode2 5.1.5: no way to snapshot the shape, so the destructor falls back to
  // reconstructing which slot to reinstall (see below).  COMPATIBILITY -- remove with
  // the follow-up issue.
  return false;
#endif
}

OdeSwapEsBatch::~OdeSwapEsBatch() {
  if (!armed_) return;
  // Restore the exact shape that was live on entry.  The snapshot covers the case
  // the old code had to special-case -- a shape installed OUTSIDE the registry, by
  // focei.R's fit-wide rxEventSensLoadModel(model$inner) -- because it puts back
  // what was there rather than reconstructing who put it there.  Leaving the
  // batch's shape live would mis-specify every inner solve after the first
  // gradient call of an iterating fit.
  bool restored = false;
#ifdef NLMIXR2EST_HAS_ESSHAPE
  restored = rxode2EventSensShapeRestore != NULL && !prevShape_.empty() &&
    rxode2EventSensShapeRestore(prevShape_.data(), (int)prevShape_.size());
#endif
  if (restored) {
    // put back verbatim -- covers a shape installed outside the registry
  } else if (prevSlotIdx_ >= 0) {
    odeSwapEsInstall(prevSlotIdx_);          // a SLOT, never a role
  } else if (prevSlot_ == odeEsInner && odeSwapHasEs(odeSlotInner)) {
    odeSwapEsInstall(odeSlotInner);
  } else {
    odeSwapEsDeactivate();
  }
  _odeEsSlot = prevSlot_;
  _odeEsSlotIdx = prevSlotIdx_;
}

void odeSwapClear(int slot) {
  if (!odeSlotOk(slot)) return;
  OdeModelReg &m = _odeReg[slot];
  if (m.fns != NULL) rxClearFuns(m.fns);
  m.fns = NULL; m.name = NULL; m.neq = 0; m.nlhs = 0; m.loaded = false;
  m.nSens = 0; m.cmtPar = -1;
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
int  odeSwapNSens(int slot)  { return odeSwapLoaded(slot) ? _odeReg[slot].nSens : 0; }
int  odeSwapCmtPar(int slot) { return odeSwapLoaded(slot) ? _odeReg[slot].cmtPar : -1; }

// How much to subtract from the pooled table's raw CMT before `slot`'s calc_lhs
// reads it.  See the OdeModelReg comment: rxode2 bakes `- nSens` into each
// model's own _CMT macro, but the translated event table belongs to whichever
// model sized the pool, so a peer de-normalizes by the wrong offset.
//
// UPSTREAM DEFICIENCY (rxode2): the shared table should carry an offset-free
// endpoint id, or _CMT should normalize against the table's model rather than
// the reading model's.  Until rxode2 offers that, this rebase is the fix.
// Returns 0 when the slot IS the pool, when either model has no CMT parameter,
// or when the counts already agree -- so it is a no-op on every single-model fit.
int odeSwapCmtDelta(int slot) {
  if (!odeSwapLoaded(slot)) return 0;
  const OdePoolPlan &p = odeSwapPlan();
  if (p.poolSlot < 0 || p.poolSlot == slot) return 0;
  if (!odeSwapLoaded(p.poolSlot)) return 0;
  return _odeReg[p.poolSlot].nSens - _odeReg[slot].nSens;
}

// Re-base the endpoint rows of one subject's CMT covariate column.  See odeSwap.h.
//
// Goes through rxode2's accessor pair getIndCmt()/setIndCmt() (function-pointer
// table), so nothing here is ABI-linked to op->cmtCov / ind->cov_ptr any more --
// setIndCmt() landed in rxode2 5.1.6 (PR #1175) as the writer this needed.  Both
// no-op for a model with no CMT covariate and for an individual with no covariate
// array, so the callers need no separate guard for those.
//
// NA rows: getIndCmt() reports a missing CMT as 1 -- it cannot distinguish NA from
// compartment 1 -- but 1 is a physical compartment for any real model, so the
// |cmt| > nPhys test below skips those rows exactly as reading the raw NA did.
// That reasoning needs nPhys >= 1, which the constructor now requires.
//
// Only this subject's rows are touched, so it stays safe inside the per-subject
// parallel regions.  Sign handling mirrors the _CMT macro: a negative CMT offsets
// the other way.
// COMPATIBILITY (remove with the follow-up issue): rxode2 5.1.6 added setIndCmt(),
// the writer half of the accessor pair.  Use it when the header has it; on 5.1.5 fall
// back to the ABI-linked field access this used to do.  Reading is always through
// getIndCmt(), which both versions have -- only the WRITE needed the new entry point.
static inline int odeCmtGet(rx_solving_options *op, rx_solving_options_ind *ind, int kk) {
  return getIndCmt(op, ind, kk);
}
static inline void odeCmtSet(rx_solving_options *op, rx_solving_options_ind *ind,
                             int kk, int cmt) {
#ifdef NLMIXR2EST_HAS_SETINDCMT
  // The #ifdef only says the header DECLARED it -- it does not say the installed
  // rxode2 FILLED that slot.  nlmixr2est compiled against 5.1.6 can be loaded against
  // 5.1.5 (DESCRIPTION allows it), where the pointer table is shorter and this stays
  // NULL from iniRxodePtrs0(); calling it would segfault.  Compile-time decides what
  // is compilable, runtime decides what is used.
  if (setIndCmt != NULL) { setIndCmt(op, ind, kk, cmt); return; }
#endif
  {
    if (op == NULL || op->cmtCov < 0 || ind == NULL || ind->cov_ptr == NULL) return;
    int n = getIndNallTimes(ind);
    if (kk < 0 || kk >= n) return;
    ind->cov_ptr[(size_t)n * (size_t)op->cmtCov + (size_t)kk] = (double) cmt;
  }
}

void OdeSwapCmtScope::shift(int by) {
  if (by == 0 || _op == NULL || _ind == NULL) return;
  int n = getIndNallTimes(_ind);
  if (n <= 0) return;
  _saved.clear();
  for (int kk = 0; kk < n; ++kk) {
    int v = odeCmtGet(_op, _ind, kk);
    // rxode2 >= 5.1.6 reports a missing CMT as NA_INTEGER; 5.1.5 reports it as 1, which
    // the |v| > _nPhys test below already skips (1 is physical for any real model, which
    // is why the constructor requires _nPhys >= 1).  Handle the new value explicitly: it
    // is INT_MIN, so `v < -_nPhys` would otherwise ACCEPT it and write back garbage.
    if (v == NA_INTEGER) continue;
    // dose / physical compartments are left alone -- the macro does not offset them
    if (!(v > _nPhys || v < -_nPhys)) continue;
    _saved.push_back(std::make_pair(kk, v));      // remember the pool-basis value
    odeCmtSet(_op, _ind, kk, (v < 0) ? v - by : v + by);
  }
}

// Put back exactly what shift() overwrote.
//
// Restoring the SAVED value rather than shifting again is what makes this exact.  Two
// ways a recomputed reverse goes wrong:
//   * re-testing |v| > _nPhys on the SHIFTED value: a pool-basis CMT just above _nPhys
//     can land at or below it after shift(-delta), so the reverse pass would classify it
//     as physical and skip it, leaving the column in the peer's basis for every later
//     reader;
//   * re-deriving the direction from the shifted value: a positive CMT that shifts below
//     zero would take the negative branch and move further away instead of back.
// Neither can arise when the original value is simply written back.
void OdeSwapCmtScope::unshift() {
  if (_op == NULL || _ind == NULL) return;
  for (size_t i = 0; i < _saved.size(); ++i) {
    odeCmtSet(_op, _ind, _saved[i].first, _saved[i].second);
  }
  _saved.clear();
}

OdeSwapCmtScope::OdeSwapCmtScope(int slot, rx_solving_options *op,
                                 rx_solving_options_ind *ind) {
  int d = odeSwapCmtDelta(slot);
  if (d == 0 || op == NULL || ind == NULL) return;   // no cmtCov -> accessors no-op
  const OdePoolPlan &p = odeSwapPlan();
  // physical compartments = states that are not sensitivities; identical for every
  // peer (they share one user model), so the pool's own counts define it
  int nPhys = odeSwapNeq(p.poolSlot) - odeSwapNSens(p.poolSlot);
  // >= 1, not >= 0: shift() relies on compartment 1 being physical so that an NA
  // CMT (which getIndCmt reports as 1) is skipped rather than shifted.
  if (nPhys < 1) return;
  _op = op; _ind = ind; _delta = d; _nPhys = nPhys;
  shift(-_delta);                 // pool basis -> this peer's basis
}

OdeSwapCmtScope::~OdeSwapCmtScope() {
  if (_delta == 0) return;
  unshift();                      // back to the pool basis for the next reader
  _delta = 0;
}

int  odeSwapNeq(int slot)    { return odeSwapLoaded(slot) ? _odeReg[slot].neq  : 0; }
int  odeSwapNlhs(int slot)   { return odeSwapLoaded(slot) ? _odeReg[slot].nlhs : 0; }
const char *odeSwapName(int slot) { return odeSwapLoaded(slot) ? _odeReg[slot].name : NULL; }

SEXP odeSwapModelSEXP(int slot) {
  if (!odeSwapLoaded(slot) || _odeModels == R_NilValue) return R_NilValue;
  return VECTOR_ELT(_odeModels, slot);
}

rxSolveF *odeSwapFns(int slot) {
  return odeSwapLoaded(slot) ? _odeReg[slot].fns : NULL;
}

void odeSwapSolveInd(int slot, int rxId) {
  rxSolveF *f = odeSwapFns(slot);
  if (f == NULL || f->dydt == NULL) return;
  ind_solve(getRxSolve_(), rxId, f->dydt_liblsoda, f->dydt_lsoda_dum,
            f->jdum_lsoda, f->dydt, f->update_inis, f->global_jt);
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
  // The registry records the pool as PLANNED; verify the LIVE solve still
  // matches.  Post-fit table generation (and any other re-solve) replaces the
  // global solve with one sized for the inner/pred model, while the registry
  // still describes the fit's pool -- solving this model there overflows the
  // state and lhs buffers (measured: post-fit analytic gradient after
  // calcTables=TRUE segfaulted exactly this way).
  rx_solve *_rx = getRxSolve_();
  if (_rx == NULL) return odeDenyPoolNotSized;
  rx_solving_options *_op = getSolvingOptions(_rx);
  if (_op == NULL ||
      getOpNeq(_op) < _odeReg[slot].neq ||
      getOpNlhs(_op) < _odeReg[slot].nlhs) return odeDenyPoolNotSized;
  return odeDenyNone;
}

// ---- per-individual solve scope -----------------------------------------

static std::atomic<long> _odeOverrideArmedN(0);
static std::atomic<long> _odeLhsWidthMismatchN(0);
static std::atomic<long> _odeOverrideNeutralizedN(0);
static std::atomic<long> _odeScratchUsedN(0);
static std::atomic<long> _odeScratchResizeN(0);
// cumulative, NOT reset by clearAll/unpin -- post-fit `pinned` is always FALSE
// because teardown unpins, so this is the only way a test can see it happened
static std::atomic<long> _odePinnedN(0);
// Counts COMPLETED pooled outer-gradient solves.  Cumulative and never reset, so
// a test can tell a working pooled path from a silent fallback to rxode2::rxSolve
// -- the two are numerically equivalent, which is how a dead path went unnoticed.
static std::atomic<long> _odePooledSolveN(0);
static std::atomic<long> _odePinCalledN(0);
static std::atomic<int>  _odePinDeny(0);

long odeSwapOverrideArmedN() { return _odeOverrideArmedN.load(std::memory_order_relaxed); }
long odeSwapLhsWidthMismatchN() { return _odeLhsWidthMismatchN.load(std::memory_order_relaxed); }
long odeSwapOverrideNeutralizedN() { return _odeOverrideNeutralizedN.load(std::memory_order_relaxed); }

bool odeSwapCheckLhsWidth(int slot, rxSolveF *fns, rx_solve *rx, rx_solving_options *op) {
  if (fns == NULL || fns->calc_lhs == NULL || rx == NULL || op == NULL) return false;
  int want = odeSwapNlhs(slot);
  int room = getOpNlhs(op);
  if (want <= 0 || room < want) return false;
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, 0);   // base subject 0
  if (ind == NULL) return false;
  double *st = getIndSolve(ind);
  if (st == NULL) return false;
  // Only WHICH slots get written matters, not the values -- so use a sentinel no
  // model produces.  NaN would compare unequal and so read as "written".
  static const double _sent = -9.87654321e37;
  std::vector<double> probe((size_t)room, _sent);
  // Probe at this subject's own first time, not t=0: an lhs sitting inside a
  // time-dependent branch would go unassigned at an arbitrary time and look missing.
  fns->calc_lhs(0, getTime(getIndIx(ind, 0), ind), st, probe.data());
  // EVERY declared slot must be written, not merely the last one -- a writer that
  // assigns lhs[0] and lhs[want-1] while skipping the middle would otherwise pass.
  int missing = 0;
  for (int k = 0; k < want; ++k) if (probe[(size_t)k] == _sent) missing++;
  if (missing > 0) {
    // A model whose lhs are all inside conditional branches could in principle
    // report missing here and lose the pooled route.  That is the safe direction:
    // the caller then takes the rxode2::rxSolve reference path, which is correct,
    // just slower -- and the counter plus the warning make it visible rather than
    // silent.  odeSwapCanPool() computes deny reasons on demand.
    if (_odeLhsWidthMismatchN.fetch_add(1, std::memory_order_relaxed) == 0) {
      Rf_warning("analytic gradient: pooled solve disabled (model/code mismatch)");
    }
    return false;
  }
  return true;
}
long odeSwapScratchUsedN()   { return _odeScratchUsedN.load(std::memory_order_relaxed); }
long odeSwapScratchResizeN() { return _odeScratchResizeN.load(std::memory_order_relaxed); }
long odeSwapPinnedN()        { return _odePinnedN.load(std::memory_order_relaxed); }
long odeSwapPooledSolveN()   { return _odePooledSolveN.load(std::memory_order_relaxed); }
void odeSwapNotePooledSolve() { _odePooledSolveN.fetch_add(1, std::memory_order_relaxed); }
long odeSwapPinCalledN()     { return _odePinCalledN.load(std::memory_order_relaxed); }
int  odeSwapPinDeny()        { return _odePinDeny.load(std::memory_order_relaxed); }

// A private buffer is needed exactly when this model writes more lhs values than
// rxode2's per-thread slice holds.  That slice is op->nlhs wide, which is the
// POOL model's width -- so a peer narrower than or equal to the pool reads
// through rxode2's own slice, and only a wider one needs our buffer.
//
// Whether that happens is a property of the model set, not of which slot it is:
// est="impmap" gives each residual-error theta a d(V)/d(theta) lhs column but no
// sensitivity state, so with more residual parameters than etas the
// theta-sensitivity model is narrower in neq yet wider in lhs than the inner
// model and the buffer is required.  See test-odeswap-fit.R.
bool odeSwapWantsScratch(int slot, rx_solving_options *op) {
  if (op == NULL) return true;              // no pool to measure against; be safe
  return odeSwapNlhs(slot) > getOpNlhs(op);
}

// One buffer per thread, grown on demand and reused across subjects, replacing
// the three per-call std::vector allocations.
static std::vector<double> &odeSwapScratch(int want) {
  static thread_local std::vector<double> buf;
  if ((int)buf.size() < want) {
    if (!buf.empty()) _odeScratchResizeN.fetch_add(1, std::memory_order_relaxed);
    buf.resize((size_t)want);
  }
  return buf;
}

OdeSwapScope::OdeSwapScope(int slot, rx_solving_options_ind *ind, rx_solving_options *op)
  : ind_(ind), op_(op), slot_(slot), saved_(-1), neq_(0), wide_(false), armed_(false) {
  neq_ = odeSwapNeq(slot);
  wide_ = odeSwapWantsScratch(slot, op);
  // Only arm for a registered slot: odeSwapNeq() returns 0 for an unloaded one,
  // and compacting the stride to 0 would corrupt the solve rather than fail.
  // AND only when the event path is bound to THIS model's role.  handle_evid
  // callocs its scratch from the EFFECTIVE neq and then calls the INSTALLED model's
  // dydt, so an override taken against a different installed model overflows that
  // scratch (valgrind-confirmed).  Unknown (odeEsUnknown) never compacts.  pred is
  // exempt and has to be: it carries no event sensitivities of its own, so the jump
  // injection does not fire, and the finite-difference fallback REQUIRES compaction.
  int _esCur = odeSwapEsInstalledModel();
  int _esWant = odeSwapEsModelForSlot(slot);
  bool _pathMatches = (_esWant == odeEsPred) ||
    (_esCur != odeEsUnknown && _esCur == _esWant);
  if (ind_ != NULL && odeSwapLoaded(slot) && _pathMatches) {
    saved_ = getIndNeqOverride(ind_);
    setIndNeqOverride(ind_, neq_);
    armed_ = true;
    _odeOverrideArmedN.fetch_add(1, std::memory_order_relaxed);
  } else if (ind_ != NULL) {
    // Declining must mean "run at the POOL's full width" -- not "inherit whatever
    // override happens to be pinned".  A focei fast fit pins the inner override for
    // the whole fit, so inheriting it would integrate only the wider model's leading
    // states and leave every later one stale.
    // ANY pinned override, not just a narrower one: a WIDER one left in place would
    // integrate more states than this model has and run its dydt out of bounds.
    int _cur = getIndNeqOverride(ind_);
    if (_cur >= 0) {
      saved_ = _cur;
      setIndNeqOverride(ind_, -1);       // -1 = no override -> full op->neq
      armed_ = true;
      _odeOverrideNeutralizedN.fetch_add(1, std::memory_order_relaxed);
    }
  }
}

OdeSwapScope::~OdeSwapScope() {
  if (armed_ && ind_ != NULL) setIndNeqOverride(ind_, saved_);
}

double *OdeSwapScope::lhs() const {
  if (!wide_) return getIndLhs(ind_);
  // Size against the live op->nlhs too: the registry can only under-predict if a
  // model was registered after the pool was built, and silently truncating there
  // would be the exact overflow this exists to prevent.
  int want = odeSwapNlhs(slot_);
  int opN = (op_ != NULL) ? getOpNlhs(op_) : 0;
  if (opN > want) want = opN;
  _odeScratchUsedN.fetch_add(1, std::memory_order_relaxed);
  return odeSwapScratch(want).data();
}

// ---- persistent pin ------------------------------------------------------

static int   _odePinnedSlot = -1;
static rx_solve *_odePinnedRx = NULL;
static int   _odePinnedNsub = 0;

bool odeSwapPinned()     { return _odePinnedSlot >= 0; }
int  odeSwapPinnedSlot() { return _odePinnedSlot; }

static void odeSwapPinWalk(int neq) {
  rx_solve *rxl = getRxSolve_();
  if (rxl == NULL) return;
  int nsub = (int)getRxNsub(rxl);
  // Base subjects only: getRxId(id) is id % nsub in inner.cpp, i.e. the identity
  // over [0, nsub), and mixture pseudo-subjects share these same inds.
  for (int id = 0; id < nsub; ++id) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rxl, id);
    if (ind != NULL) setIndNeqOverride(ind, neq);
  }
  _odePinnedRx = rxl;
  _odePinnedNsub = nsub;
}

void odeSwapPinAll(int slot) {
  _odePinCalledN.fetch_add(1, std::memory_order_relaxed);
  if (!odeSwapLoaded(slot)) { _odePinDeny.store(1, std::memory_order_relaxed); return; }
  const OdePoolPlan &p = odeSwapPlan();
  int neq = odeSwapNeq(slot);
  // Only meaningful when a larger peer sized the pool; otherwise the stride
  // already matches and pinning would just be state to unwind.
  if (p.poolSlot < 0) { _odePinDeny.store(2, std::memory_order_relaxed); return; }
  if (neq >= p.poolNeq) { _odePinDeny.store(3, std::memory_order_relaxed); return; }
  odeSwapPinWalk(odeSwapNeq(slot));
  _odePinnedSlot = slot;
  _odePinnedN.fetch_add(1, std::memory_order_relaxed);
}

void odeSwapRepin() {
  if (_odePinnedSlot < 0) return;
  odeSwapPinWalk(odeSwapNeq(_odePinnedSlot));
}

void odeSwapUnpinAll() {
  if (_odePinnedSlot < 0) return;
  // Only walk when this is still the solve structure we pinned: after
  // rxSolveFree() the inds are gone, and after a rebuild they are different
  // objects that were never pinned.
  if (_odePinnedRx != NULL && getRxSolve_() == _odePinnedRx) {
    int nsub = (int)getRxNsub(_odePinnedRx);
    if (nsub > _odePinnedNsub) nsub = _odePinnedNsub;
    for (int id = 0; id < nsub; ++id) {
      rx_solving_options_ind *ind = getSolvingOptionsInd(_odePinnedRx, id);
      if (ind != NULL) setIndNeqOverride(ind, -1);
    }
  }
  _odePinnedSlot = -1; _odePinnedRx = NULL; _odePinnedNsub = 0;
}

// ---- override-aware solve-buffer scanning -------------------------------

int odeSwapIndSolveSpan(rx_solving_options *op, rx_solving_options_ind *ind) {
  if (op == NULL || ind == NULL) return 0;
  int full = getOpNeq(op);
  int eff = odeSwapEffNeq(ind, op);
  int nlin = getOpNlin(op);
  // rxode2 lays getAdvan() out from the FULL neq when numLin > 0, so an override
  // does not compact that part of the buffer.  Scan the full span there rather
  // than risk under-scanning a region the solve did write.
  if (nlin > 0 && eff != full) eff = full;
  int n = eff + nlin;
  if (n <= 0) return 0;
  return n * getIndNallTimes(ind);
}

bool odeSwapIndBadSolveSlot(rx_solving_options *op, rx_solving_options_ind *ind,
                            int slot) {
  int n = odeSwapIndSolveSpan(op, ind);
  int own = odeSwapNeq(slot);
  // Only clamp DOWN, and only when this model is genuinely narrower than the
  // span: never widen, or a compacted solve would be under-scanned.
  if (own > 0 && own < n && getOpNlin(op) <= 0) n = own;
  if (n <= 0) return false;
  double *solve = getIndSolve(ind);
  if (solve == NULL) return false;
  for (int i = 0; i < n; ++i) {
    if (ISNA(solve[i]) || std::isnan(solve[i]) || std::isinf(solve[i])) return true;
  }
  return false;
}

bool odeSwapIndBadSolve(rx_solving_options *op, rx_solving_options_ind *ind) {
  int n = odeSwapIndSolveSpan(op, ind);
  if (n <= 0) return false;
  double *solve = getIndSolve(ind);
  if (solve == NULL) return false;
  for (int i = 0; i < n; ++i) {
    if (ISNA(solve[i]) || std::isnan(solve[i]) || std::isinf(solve[i])) return true;
  }
  return false;
}

// ---- R introspection ----------------------------------------------------
//
// Exposes the pool decision so tests can assert the MECHANISM (which model sized
// the pool, whether the scratch was needed, whether an override leaked), not just
// that the numbers came out the same.

// Drives odeSwapRetryCore -- the REAL loop, not a copy -- with stub side effects,
// so the retry logic is testable without a pathological ODE.  `nFail` solves
// report bad, then they succeed.
//[[Rcpp::export]]
List odeSwapRetryTest_(int nFail, int maxOdeRecalc, int stickyRecalcN,
                       double odeRecalcFactor, int relaxMode, int sticky0,
                       bool restoreTolOnSuccess) {
  OdeRetryOpts o;
  o.maxOdeRecalc = maxOdeRecalc;
  o.stickyRecalcN = stickyRecalcN;
  o.odeRecalcFactor = odeRecalcFactor;
  o.relaxMode = relaxMode;
  o.restoreTolOnSuccess = restoreTolOnSuccess;

  int solves = 0, indRelax = 0, globalRelax = 0, onRetry = 0, onSticky = 0;
  double tol = 1.0;
  int sticky = sticky0;

  struct StubHooks {
    int *r; int *s;
    void onRetry() { (*r)++; }
    void onSticky() { (*s)++; }
  } h; h.r = &onRetry; h.s = &onSticky;

  int retries = odeSwapRetryCore(
    sticky,
    [&]{ solves++; },
    [&]{ return solves <= nFail; },              // the first nFail solves "fail"
    [&](int mode) { if (mode == odeRelaxInd) { indRelax++; tol *= o.odeRecalcFactor; }
                    else globalRelax++; },
    [&]{ return tol; },
    [&](double x){ tol = x; },
    o, h);

  return List::create(_["retries"] = retries, _["solves"] = solves,
                      _["tolFactor"] = tol, _["stickyRecalcN2"] = sticky,
                      _["indRelax"] = indRelax, _["globalRelax"] = globalRelax,
                      _["onRetry"] = onRetry, _["onSticky"] = onSticky);
}

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

//' Record which model role rxode2's event path is bound to (R-side installs).
//' Roles: -1 unknown, 0 pred, 1 inner, 2 outer, 3 hess2.
//' @param slot role id
//' @return NULL
//' @noRd
//[[Rcpp::export]]
RObject odeSwapEsNoteInstalled_(int slot) {
  odeSwapEsNoteInstalled(slot);   // a ROLE, not a slot index
  return R_NilValue;
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
  IntegerVector activeOv(0);
  if (rxl != NULL) {
    rx_solving_options *op = getSolvingOptions(rxl);
    if (op != NULL) { opNeq = getOpNeq(op); opNlhs = getOpNlhs(op); }
    int nsub = (int)getRxNsub(rxl);
    if (nsub > 0) {
      activeOv = IntegerVector(nsub);
      for (int id = 0; id < nsub; ++id) {
        rx_solving_options_ind *ind = getSolvingOptionsInd(rxl, id);
        activeOv[id] = (ind == NULL) ? NA_INTEGER : getIndNeqOverride(ind);
      }
    }
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
    _["opNeq"] = opNeq, _["opNlhs"] = opNlhs,
    _["pinned"] = odeSwapPinned(),
    _["pinnedSlot"] = odeSwapPinnedSlot(),
    // -1 per subject means "no override in force"; a non -1 left behind after a
    // fit is the leak this component exists to make impossible.
    _["activeOverride"] = activeOv,
    _["overrideArmedN"] = (double)odeSwapOverrideArmedN(),
    _["overrideNeutralizedN"] = (double)odeSwapOverrideNeutralizedN(),
    _["lhsWidthMismatchN"] = (double)odeSwapLhsWidthMismatchN(),
    _["scratchUsedN"] = (double)odeSwapScratchUsedN(),
    _["scratchResizeN"] = (double)odeSwapScratchResizeN(),
    _["pinnedN"] = (double)odeSwapPinnedN(),
    _["pooledSolveN"] = (double)odeSwapPooledSolveN(),
    _["pinCalledN"] = (double)odeSwapPinCalledN(),
    _["pinDeny"] = (double)odeSwapPinDeny());
}
