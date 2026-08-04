#ifndef __ODESWAP_H__
#define __ODESWAP_H__
// Shared-pool ODE model swapping.
//
// The FOCEi-family inner problem runs several DIFFERENT compiled rxode2 models
// inside ONE rxode2 solve pool:
//
//   1. the pool is sized (rxSolve_) for the registered model with the largest
//      neq -- see odeSwapPlan();
//   2. a single individual is solved as another registered model by swapping
//      ind->neqOverride (rxode2 honors it via rxEffNeq); allocations stay at
//      op->neq and only the per-event solve stride compacts, so this is safe
//      per individual and never touches the shared op->neq;
//   3. results are read back with calc_lhs.
//
// Step 3 is the subtle one: rxode2's per-thread ind->lhs slice is exactly
// op->nlhs wide and rxode2 has NO lhs override, so reading a model whose lhs is
// wider than the pool model's would run off the end of this thread's slice.
// That happens whenever the widest-lhs model is not the largest-neq model.  It
// is reachable from an ordinary model: est="impmap" gives each residual-error
// theta a d(V)/d(theta) lhs column but NO sensitivity state, so once a model has
// more residual parameters than etas the theta-sensitivity model is narrower in
// neq yet wider in lhs than the inner model (measured: 2 etas with
// add+prop+boxCox gives inner 6/6 against thetaSens 4/8).
// odeSwapPlan() detects it and reports scratchNlhs; OdeSwapScope::lhs() then
// hands back a private buffer instead of rxode2's slice, so a call site can
// never pick the wrong one.
//
// Assumes armahead.h (hence rxode2ptr.h) and inner.h were included first.  This
// header must NOT reference op_focei/nlmOp -- those are private to inner.cpp and
// nlm.cpp respectively, which is exactly what forced the duplication this
// component removes.
#include "inner.h"
#include <vector>
#include <utility>
#include <string>

// Registered peer solvers.  Add a slot here and register it -- nothing else in
// the pool-sizing or buffer logic needs to change.
enum OdeSwapSlot {
  odeSlotInner = 0,   // rxInner:     focei inner (eta sensitivities) / nlm thetaGrad
  odeSlotPred,        // rxPred:      predNoLhs FD model / nlm predOnly
  odeSlotThetaSens,   // rxThetaSens: impmap / advi d(f)/d(theta)
  odeSlotHess2,       // rxHess2:     fast=TRUE ll() d2(logLik)/deta2
  odeSlotOuter,       // rxVaeOuter:     augmented outer-gradient model (order 2)
  odeSlotOuterNode,   // rxOuterNode:    same directions at order 1, for AGQ nodes
  odeSlotOuterCov,    // rxOuterCov:     covariance model over its own direction set
  odeSlotN
};

// The analytic path compiles up to THREE augmented models for one fit -- the
// order-2 gradient model, an order-1 model for AGQ quadrature nodes (fast=TRUE
// with nAGQ > 1), and a covariance model over a different direction set.  They
// are distinct cache keys (digest | dirs | order | flags), hence distinct
// compiles, and with nAGQ > 1 the first two are solved during the SAME
// optimization.  Each therefore needs its own rxSolveF and its own slot: they all
// stay registered together for the life of the fit, and solving one is a matter
// of calling ITS entry points, not of re-registering a shared slot.  The pool is
// still built exactly once, sized by odeSwapPlan() for whichever of them has the
// most ODE states.

// Why a pooled solve was refused, so a fallback is loud rather than silent.
enum OdeSwapDeny {
  odeDenyNone = 0,
  odeDenyNotLoaded,      // slot never registered / calc_lhs NULL
  odeDenyPoolNotSized,   // this model needs more states than the pool has
  odeDenyNlhsUnknown     // lhs width unknown, so no buffer can be sized
};

// ---- registration -------------------------------------------------------

// Record a model's neq/nlhs/lhs names WITHOUT loading its DLL or binding its
// entry points.  That is all odeSwapPlan() needs, and it must stay metadata-only:
// rxDynLoad-ing a sensitivity model before rxSolve_ has built the pool rebinds
// rxode2's event-sensitivity globals to it, and the inner solve then frees an
// ES buffer sized for the wrong neq ("free(): invalid next size").
bool odeSwapDeclare(int slot, const char *name, SEXP obj);

// odeSwapDeclare plus rxDynLoad + rxUpdateFuns, for when the entry points are
// actually going to be called.  Returns false (leaving the slot unloaded) when
// obj is not a usable rxode2 model.
bool odeSwapRegister(int slot, const char *name, SEXP obj, rxSolveF *fns);
void odeSwapClear(int slot);
void odeSwapClearAll();          // clears every slot and releases preserved SEXPs

bool odeSwapLoaded(int slot);
// ---- event-sensitivity (jump) shape ---------------------------------------
// rxode2 keeps the ES shape in PROCESS GLOBALS (_rxEsActive/_rxEsNState/
// _rxEsNParam/_rxEsNParam2), installed through the registered C-callable
// rxode2EventSensLoad().  It therefore describes exactly ONE model at a time:
// rxPred has no event sensitivities, rxInner has them, and the augmented outer
// model has them with a DIFFERENT shape.  Solving one model under another's
// shape mis-specifies its jumps.
//
// Because the state is global it cannot be swapped per individual or inside an
// OpenMP region.  Solves are batched by model: install the batch's shape, solve
// every individual for that model, then restore.  OdeSwapEsBatch is that
// boundary and MUST be constructed outside any parallel region -- unlike
// OdeSwapScope, which is per-ind and safe inside one.
//
// The shape belongs to a MODEL ROLE, not to a registry slot: several slots share
// one shape (thetaSens and the augmented outer models), while hess2 has its own.
// handle_evid sizes its scratch from the effective neq but calls the INSTALLED
// model's dydt, so a solve may only be compacted when the two agree -- and the
// installer can be R-side (focei.R loads the inner model's shape for the whole
// fit), which the registry cannot see on its own.  Hence odeSwapEsNoteInstalled().
enum OdeEsModel {
  odeEsUnknown = -1,   // not recorded -- never compact
  odeEsPred    =  0,   // no event sensitivities of its own
  odeEsInner   =  1,
  odeEsOuter   =  2,   // augmented outer models (and thetaSens)
  odeEsHess2   =  3    // its OWN shape, NOT inner's
};
// Which ROLE a slot's event-sensitivity shape belongs to.
int  odeSwapEsModelForSlot(int slot);
// Which ROLE is installed right now (odeEsUnknown when we have not recorded it).
int  odeSwapEsInstalledModel();
// Record a role installed outside the registry (focei.R's fit-wide load).
void odeSwapEsNoteInstalled(int esModel);

struct OdeSwapEsBatch {
  explicit OdeSwapEsBatch(int slot);
  ~OdeSwapEsBatch();
  bool armed() const { return armed_; }
private:
  bool saveLive();
  int prevSlot_;      // previous ROLE
  int prevSlotIdx_;   // previous SLOT -- fallback when the snapshot is unavailable
  std::vector<char> prevShape_;   // the shape live on entry (rxode2EventSensShapeSave)
  bool armed_;
  OdeSwapEsBatch(const OdeSwapEsBatch &);
  OdeSwapEsBatch &operator=(const OdeSwapEsBatch &);
};
// Does this slot carry event sensitivities at all?  (false for rxPred.)
bool odeSwapHasEs(int slot);

int  odeSwapNeq(int slot);       // 0 when unloaded; matches rxode2's op->neq
int  odeSwapNlhs(int slot);      // 0 when unloaded
int  odeSwapNSens(int slot);     // length($sens): sensitivity compartments
int  odeSwapCmtPar(int slot);    // index of "CMT" in $params, -1 when absent

// Endpoint (CMT) rebasing for a pooled solve.
//
// rxode2 compiles a multi-endpoint model's endpoint switch in USER compartment
// numbering and normalizes the raw CMT with the model's OWN sensitivity count:
//     #define _CMT ((fabs(CMT)<=nPhys) ? CMT : CMT - nSens)
// Each model is therefore self-consistent on its own translated event table
// (inner: raw 5 -> 3; outer: raw 63 -> 3).  But a pooled fit translates ONCE,
// against whichever model sized the pool, so a peer subtracts the wrong offset
// (inner reading the outer model's table: 63 - 2 = 61, matching no endpoint) and
// its rx_pred_/rx_r_/d(f)/d(eta)/rx_yj_ all silently evaluate to 0.
//
// odeSwapCmtDelta(slot) is how much the pooled table's raw CMT is off for that
// slot: 0 unless a LARGER peer sized the pool and the two disagree.
//
// It is NOT applied anywhere yet.  Everything it needs is already available --
// nSens comes straight from rxModelVars(obj)$sens, recorded at declare time -- so
// this is purely ours to finish; nothing is required from rxode2.  What is not
// simple is WHERE to apply it: rxode2 delivers CMT as a COVARIATE (op->cmtCov,
// rxData.cpp) and refreshes it from the covariate arrays at the top of calc_lhs,
// so writing ind->par_ptr beforehand is overwritten before the model reads it
// (measured).  The interception point has to be the covariate column, for the
// duration of a pooled solve.  Until that is built, multi-endpoint models simply
// do not pool -- see the outerPoolOk gate in R/focei.R.
//
// NOTE the _CMT normalization itself is CORRECT and must not be "fixed" upstream:
// it is what makes a model solved against its OWN translated table read the right
// endpoint, which is every standalone solve (npde, cwres, tables, plain rxSolve).
// What is unsupported is running several peers of different sensitivity depth over
// ONE table -- that is the pool's own doing, so the rebase belongs here and applies
// to pooled solves only.
int  odeSwapCmtDelta(int slot);

// RAII: re-base one subject's CMT covariate column for the duration of a pooled
// solve/read with `slot`'s model, and restore it on every exit path.
//
// Only ENDPOINT rows are touched -- those whose |CMT| exceeds the physical
// compartment count, which is the same `baseSize` rxode2's _CMT macro tests
// against and is identical across peers (they share one user model).  Dose rows
// (CMT <= nPhys) must NOT be shifted: the macro leaves them alone, so rebasing
// them would send a dose into a different compartment.
//
// A no-op unless a larger peer sized the pool (odeSwapCmtDelta == 0), so plain
// single-model fits and every standalone solve pay nothing.
struct OdeSwapCmtScope {
  OdeSwapCmtScope(int slot, rx_solving_options *op, rx_solving_options_ind *ind);
  ~OdeSwapCmtScope();
  OdeSwapCmtScope(const OdeSwapCmtScope &) = delete;
  OdeSwapCmtScope &operator=(const OdeSwapCmtScope &) = delete;
  bool active() const { return _delta != 0; }
private:
  rx_solving_options *_op = NULL;
  rx_solving_options_ind *_ind = NULL;
  int _delta = 0;
  int _nPhys = 0;
  void shift(int by);
  void unshift();                                  // restore what shift() overwrote
  std::vector<std::pair<int,int> > _saved;         // (row, pool-basis CMT) shift() wrote
};
const char *odeSwapName(int slot);
SEXP odeSwapModelSEXP(int slot); // R_NilValue when unloaded

// 0-based index of an lhs output in this model, or -1 when absent.  Replaces the
// hand-rolled "scan mvts$lhs for this name" loops.
int odeSwapLhsIndex(int slot, const char *nm);

// This slot's compiled entry points, or NULL when unloaded.  Solving "as model X"
// means calling X's own dydt/calc_lhs through here.  Models never displace one
// another -- each keeps its own rxSolveF for the life of the fit.
rxSolveF *odeSwapFns(int slot);

// Integrate ONE individual with `slot`'s entry points, in the shared pool.
// Generalizes the per-model `#define <x>Ode(id) ind_solve(...)` macros.
void odeSwapSolveInd(int slot, int rxId);


// ---- the pool decision --------------------------------------------------

struct OdePoolPlan {
  int poolSlot = -1;        // sizes the pool: max neq, ties -> max nlhs, then lowest slot
  int poolNeq = 0;          // becomes op->neq once rxSolve_ has run
  int poolNlhs = 0;         // becomes op->nlhs
  int maxNlhs = 0;
  int maxNlhsSlot = -1;
  int scratchNlhs = 0;      // maxNlhs when maxNlhs > poolNlhs, else 0 (rxode2's slice suffices)
  int nLoaded = 0;
  bool overrideNeeded = false;  // some loaded slot has neq < poolNeq
};

const OdePoolPlan &odeSwapPlan();
SEXP odeSwapPoolModelSEXP();

// Pure form, so the tie-break and the scratch inversion are testable without a fit.
OdePoolPlan odeSwapPlanFor(const std::vector<int> &neq, const std::vector<int> &nlhs);

int odeSwapCanPool(int slot);   // OdeSwapDeny reason code

// ---- per-individual solve scope ----------------------------------------

// True when reading `slot` needs a private lhs buffer rather than rxode2's
// per-thread ind->lhs slice.  One place decides this, so a call site cannot pick
// the wrong buffer.
bool odeSwapWantsScratch(int slot, rx_solving_options *op);

// Arms ind->neqOverride for `slot` for as long as this object lives, and hands
// back the lhs buffer to read that model through.
//
// Solve one individual as a different model:
//   OdeSwapScope sc(odeSlotThetaSens, ind, op);
//   ... solve ...
//   iniSubjectE(rxId, 1, ind, op, rx, rxThetaSens.update_inis);
//   double *lhs = sc.lhs();          // AFTER iniSubjectE -- see below
//   rxThetaSens.calc_lhs(rxId, t, getOpIndSolve(op, ind, j), lhs);
//
// lhs() is deliberately lazy: iniSubjectE(..., inLhs=1, ...) re-points ind->lhs
// at the CALLING thread's slice, and callers construct the scope well before
// they call it.  Caching the pointer in the constructor would hand back another
// thread's slice under the parallel inner loop.
struct OdeSwapScope {
  OdeSwapScope(int slot, rx_solving_options_ind *ind, rx_solving_options *op);
  ~OdeSwapScope();
  double *lhs() const;
  int neq() const { return neq_; }
  int slot() const { return slot_; }
  bool usesScratch() const { return wide_; }
  void disarm() { armed_ = false; }
private:
  OdeSwapScope(const OdeSwapScope &);
  OdeSwapScope &operator=(const OdeSwapScope &);
  rx_solving_options_ind *ind_;
  rx_solving_options *op_;
  int slot_, saved_, neq_;
  bool wide_, armed_;
};

// ---- override-aware solve-buffer scanning -------------------------------

// rxode2's rxEffNeq is compiled out for us (rxode2.h guards it behind
// __RXODE2PTR_H__), so carry the same rule: the override applies only when it is
// within [0, op->neq], otherwise the full state count does.
static inline int odeSwapEffNeq(rx_solving_options_ind *ind, rx_solving_options *op) {
  int o = getIndNeqOverride(ind), n = getOpNeq(op);
  return (o >= 0 && o <= n) ? o : n;
}

// Doubles this individual's last solve actually wrote.  Scanning further reads
// slots this solve never touched -- STALE values a previous, wider solve left in
// the shared buffer, not uninitialised memory (rxode2 allocates ind->solve once
// for op->neq and reuses it, so valgrind cannot see this; a run on origin/main
// reports zero errors).  A stale NaN/Inf there is reported as a failed solve
// that did not happen, which loosens ODE tolerances needlessly.
int odeSwapIndSolveSpan(rx_solving_options *op, rx_solving_options_ind *ind);

// NaN/Inf anywhere in that span.  Per-individual, so it answers "did THIS subject
// fail" without reading op->badSolve, which another thread can flip mid-loop.
bool odeSwapIndBadSolve(rx_solving_options *op, rx_solving_options_ind *ind);

// Same, but bounded by what the model in `slot` actually writes.  Use this when
// the pool is sized for a LARGER peer and no override is armed: the span is then
// op->neq, so the scan runs past this model's own states into slots it never
// wrote.  That is how pool SIZE leaks into results -- a stale NaN in the tail
// reads as a failed solve, loosens tolerances, and shifts the EBEs, so the same
// fit gives different answers depending on which peer happened to size the pool.
// Bounding the scan makes the pool size irrelevant without changing any stride
// (the solve and every read keep using op->neq, so all indexing is untouched).
bool odeSwapIndBadSolveSlot(rx_solving_options *op, rx_solving_options_ind *ind,
                            int slot);

// ---- bad-solve retry with tolerance relaxation --------------------------

// How a retry loosens tolerances.  Carried per call site rather than unified:
// the global form mutates rxode2's shared grtol2/gatol2 arrays and races under a
// parallel loop, while the per-individual form is thread-safe but not
// bit-identical to it.  Which one a site uses is a behavioural choice, so this
// helper preserves it rather than picking one.
enum OdeRelaxMode { odeRelaxGlobal = 0, odeRelaxInd = 1 };

struct OdeRetryOpts {
  int maxOdeRecalc = 0;
  int stickyRecalcN = 0;
  int relaxMode = odeRelaxGlobal;
  double odeRecalcFactor = 1.0;
  // FOCEi un-sticks a subject that recovered within budget; nlm deliberately
  // leaves it loosened ("stiff subjects retain loosened tolerance").
  bool restoreTolOnSuccess = true;
  // FOCEi clears op->badSolve before each re-solve; nlm does not.  op->badSolve
  // is racy under cores>1 and this helper does not read it (odeSwapIndBadSolve is
  // per-individual), so this only keeps the flag tidy for diagnostics.
  bool resetBadSolveEachRetry = true;
};

// The retry LOGIC, with every side effect injected.  Splitting it this way is
// what makes the loop testable: a unit test supplies stub solve/bad/relax/tol
// functions and drives this exact code, rather than a copy of it.
//
//   solve()      perform one solve
//   bad()        did that solve fail?
//   relax(mode)  loosen tolerances (global or per-individual)
//   tolGet()     current per-individual tolerance factor
//   tolSet(x)    restore it
template <typename SolveFn, typename BadFn, typename RelaxFn,
          typename TolGet, typename TolSet, typename Hooks>
int odeSwapRetryCore(int &stickyRecalcN2, SolveFn solve, BadFn bad, RelaxFn relax,
                     TolGet tolGet, TolSet tolSet, const OdeRetryOpts &o, Hooks &h) {
  double prevTol = tolGet();
  solve();
  int j = 0;
  while (stickyRecalcN2 <= o.stickyRecalcN && bad() && j < o.maxOdeRecalc) {
    stickyRecalcN2++;
    h.onRetry();
    relax(o.relaxMode);
    solve();
    j++;
  }
  if (j != 0) {
    if (stickyRecalcN2 <= o.stickyRecalcN) {
      if (o.restoreTolOnSuccess) tolSet(prevTol);
    } else {
      h.onSticky();
    }
  }
  return j;
}

// Solve, then while the solve is bad and budget remains, loosen and re-solve.
// Returns the retry count.  `stickyRecalcN2` is the caller's per-subject counter
// (op_focei's atomic or nlmOp's plain int -- hence the template).  Hooks supplies
// onRetry() ("tolerances were reduced") and onSticky() ("budget exhausted, the
// loosening is now permanent").
template <typename SolveFn, typename Hooks>
int odeSwapSolveRetry(rx_solving_options *op, rx_solving_options_ind *ind,
                      int &stickyRecalcN2, SolveFn solveFn,
                      const OdeRetryOpts &o, Hooks &h) {
  return odeSwapRetryCore(
    stickyRecalcN2,
    [&]{ solveFn(); },
    [&]{ return odeSwapIndBadSolve(op, ind); },
    [&](int mode) {
      if (mode == odeRelaxInd) {
        setIndTolFactor(ind, getIndTolFactor(ind) * o.odeRecalcFactor);
      } else {
        atolRtolFactor_(o.odeRecalcFactor);
      }
      setIndSolve(ind, -1);
      if (o.resetBadSolveEachRetry) resetOpBadSolve(op);
    },
    [&]{ return getIndTolFactor(ind); },
    [&](double x){ setIndTolFactor(ind, x); },
    o, h);
}

// ---- persistent pin ------------------------------------------------------

// Pin EVERY subject's effective state count to `slot`'s, for methods whose inner
// solves run compacted against a pool sized for a larger peer.  Unlike
// OdeSwapScope this outlives a single solve, so it must be balanced by
// odeSwapUnpinAll().
//
// No-op unless the slot is genuinely smaller than the pool.  Records the
// rx_solve* it pinned, so unpinning cannot walk a different (rebuilt) solve
// structure and cannot depend on caller state that may already have been reset.
void odeSwapPinAll(int slot);
void odeSwapUnpinAll();     // idempotent; safe after the pool was freed or rebuilt
void odeSwapRepin();        // re-apply after rxSolve_ rebuilt the solve structure
bool odeSwapPinned();
int  odeSwapPinnedSlot();

// Usage counters, so tests can assert the mechanism ran rather than infer it.
long odeSwapOverrideArmedN();
// A bound calc_lhs whose written width disagreed with the registry's nlhs.
long odeSwapLhsWidthMismatchN();
// Declined compactions that had to neutralize a narrower pinned override.
long odeSwapOverrideNeutralizedN();
// Verify the function pointers bound for `slot` really write the registry's nlhs.
// rxUpdateFuns resolves symbols with R_GetCCallable(lib, name) -- BY NAME -- so the
// same model can re-resolve to another dll's symbol later in a session.  Returns
// false on a mismatch, and the caller must then refuse to pool.
bool odeSwapCheckLhsWidth(int slot, rxSolveF *fns, rx_solve *rx, rx_solving_options *op);
// Number of times a bound calc_lhs was found NOT to match its registry width.
long odeSwapLhsWidthMismatchN();
// Verify the function pointers bound for `slot` really write the registry's nlhs.
// R_GetCCallable() resolves by symbol NAME, so a peer can silently re-resolve to a
// different model's dll within one session.  Returns false on a mismatch (caller
// must refuse to pool); the probe runs once per call and costs one calc_lhs.
bool odeSwapCheckLhsWidth(int slot, rxSolveF *fns, rx_solve *rx, rx_solving_options *op);
long odeSwapScratchUsedN();
long odeSwapScratchResizeN();
long odeSwapPinnedN();
long odeSwapPooledSolveN();    // completed pooled outer solves; survives teardown
void odeSwapNotePooledSolve();
long odeSwapPinCalledN();
int  odeSwapPinDeny();   // cumulative pins; survives teardown, unlike odeSwapPinned()

#endif // __ODESWAP_H__
