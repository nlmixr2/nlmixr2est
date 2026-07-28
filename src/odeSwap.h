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
int  odeSwapNeq(int slot);       // 0 when unloaded
int  odeSwapNlhs(int slot);      // 0 when unloaded
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

// Solve, then while the solve is bad and budget remains, loosen and re-solve.
// Returns the retry count.  `stickyRecalcN2` is the caller's per-subject counter
// (op_focei's atomic or nlmOp's plain int -- hence the template).  Hooks supplies
// onRetry() ("tolerances were reduced") and onSticky() ("budget exhausted, the
// loosening is now permanent").
template <typename SolveFn, typename Hooks>
int odeSwapSolveRetry(rx_solving_options *op, rx_solving_options_ind *ind,
                      int &stickyRecalcN2, SolveFn solveFn,
                      const OdeRetryOpts &o, Hooks &h) {
  double prevTol = getIndTolFactor(ind);
  solveFn();
  int j = 0;
  while (stickyRecalcN2 <= o.stickyRecalcN && odeSwapIndBadSolve(op, ind) &&
         j < o.maxOdeRecalc) {
    stickyRecalcN2++;
    h.onRetry();
    if (o.relaxMode == odeRelaxInd) {
      setIndTolFactor(ind, getIndTolFactor(ind) * o.odeRecalcFactor);
    } else {
      atolRtolFactor_(o.odeRecalcFactor);
    }
    setIndSolve(ind, -1);
    if (o.resetBadSolveEachRetry) resetOpBadSolve(op);
    solveFn();
    j++;
  }
  if (j != 0) {
    if (stickyRecalcN2 <= o.stickyRecalcN) {
      if (o.restoreTolOnSuccess) setIndTolFactor(ind, prevTol);
    } else {
      h.onSticky();
    }
  }
  return j;
}

// ---- persistent pin ------------------------------------------------------

// Pin EVERY subject's effective state count to `slot`'s, for methods whose inner
// solves run compacted against a pool sized for a larger peer.  Unlike
// OdeSwapScope this outlives a single solve, so it must be balanced by
// odeSwapUnpinAll().
//
// No-op unless the slot is genuinely smaller than the pool.  Records the
// rx_solve* it pinned, so unpinning cannot walk a different (rebuilt) solve
// structure and cannot depend on caller state that may already have been reset --
// the old clear early-returned on op_focei.innerNeq <= 0, so zeroing that first
// silently skipped it.
void odeSwapPinAll(int slot);
void odeSwapUnpinAll();     // idempotent; safe after the pool was freed or rebuilt
void odeSwapRepin();        // re-apply after rxSolve_ rebuilt the solve structure
bool odeSwapPinned();
int  odeSwapPinnedSlot();

// Usage counters, so tests can assert the mechanism ran rather than infer it.
long odeSwapOverrideArmedN();
long odeSwapScratchUsedN();
long odeSwapScratchResizeN();
void odeSwapResetCounters();

#endif // __ODESWAP_H__
