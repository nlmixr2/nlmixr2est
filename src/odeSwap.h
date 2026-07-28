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
  odeSlotOuter,       // rxVaeOuter:  augmented outer-gradient model
  odeSlotN
};

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
// slots the solve never touched, which is both a use of uninitialised memory and
// a way to false-trigger the ODE retry loop.
int odeSwapIndSolveSpan(rx_solving_options *op, rx_solving_options_ind *ind);

// NaN/Inf anywhere in that span.  Per-individual, so it answers "did THIS subject
// fail" without reading op->badSolve, which another thread can flip mid-loop.
bool odeSwapIndBadSolve(rx_solving_options *op, rx_solving_options_ind *ind);

// Usage counters, so tests can assert the mechanism ran rather than infer it.
long odeSwapOverrideArmedN();
long odeSwapScratchUsedN();
long odeSwapScratchResizeN();
void odeSwapResetCounters();

#endif // __ODESWAP_H__
