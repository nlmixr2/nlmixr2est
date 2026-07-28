// [[Rcpp::plugins(openmp)]]
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "inner.h"
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
  return odeDenyNone;
}

// ---- per-individual solve scope -----------------------------------------

static std::atomic<long> _odeOverrideArmedN(0);
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
long odeSwapScratchUsedN()   { return _odeScratchUsedN.load(std::memory_order_relaxed); }
long odeSwapScratchResizeN() { return _odeScratchResizeN.load(std::memory_order_relaxed); }
long odeSwapPinnedN()        { return _odePinnedN.load(std::memory_order_relaxed); }
long odeSwapPooledSolveN()   { return _odePooledSolveN.load(std::memory_order_relaxed); }
void odeSwapNotePooledSolve() { _odePooledSolveN.fetch_add(1, std::memory_order_relaxed); }
long odeSwapPinCalledN()     { return _odePinCalledN.load(std::memory_order_relaxed); }
int  odeSwapPinDeny()        { return _odePinDeny.load(std::memory_order_relaxed); }
void odeSwapResetCounters() {
  _odeOverrideArmedN.store(0, std::memory_order_relaxed);
  _odeScratchUsedN.store(0, std::memory_order_relaxed);
  _odeScratchResizeN.store(0, std::memory_order_relaxed);
}

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
  if (ind_ != NULL && odeSwapLoaded(slot)) {
    saved_ = getIndNeqOverride(ind_);
    setIndNeqOverride(ind_, neq_);
    armed_ = true;
    _odeOverrideArmedN.fetch_add(1, std::memory_order_relaxed);
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
    _["scratchUsedN"] = (double)odeSwapScratchUsedN(),
    _["scratchResizeN"] = (double)odeSwapScratchResizeN(),
    _["pinnedN"] = (double)odeSwapPinnedN(),
    _["pooledSolveN"] = (double)odeSwapPooledSolveN(),
    _["pinCalledN"] = (double)odeSwapPinCalledN(),
    _["pinDeny"] = (double)odeSwapPinDeny());
}
