# Phase C / Phase 9: per-individual stride + lhs sizing (UNFINISHED)

Status: NOT implemented.  Branch `feat/dry-ode-swap` is clean at 1144cd2ac
(pushed).  Three attempts were made and reverted; what each proved is below so
the next attempt does not repeat them.

## The bug being fixed

`test-focei-fast-grad.R` -> "fast=TRUE fit matches the finite-difference fit"
fails identically on this branch and on origin/main (44 pass / 3 fail / 0 err;
the 3 are assertions in that one test).  It passes in isolation and fails once
block 6 ("analytic outer gradient matches FD for a multiple-endpoint model")
has run.  Reproducer: `subRun.sh 7` fails, `subRun.sh 6` passes.

Measured cause: the shared pool is sized for the INNER model (op->neq = 8)
while the augmented outer model needs 26 states.  rxEffNeq() only ever
COMPACTS, so the augmented solve is truncated to its first-order states.  The
truncation boundary is exact:

    last finite : rx_f1_ETA_3_          (first order; 2 base states x (1+3 etas) = 8)
    first NaN   : rx_f2_ETA_1__ETA_1_   (second order)

25 of 57 columns come back all-NaN, including rx_rvarf_, so R and Rsig are NaN
for every observation, the gradient is NaN, .foceiAnalyticGradCore returns NULL
at foceiGradAnalytic.R:884, and the optimizer silently uses finite differences.

Ruled out by direct measurement: scope downgrade (fast=TRUE does reach the
fit), augmented-model build (amOK=TRUE, cached), column mapping (identical
prefix/nlhs=29/rvarfIdx=11 in the passing isolated run), non-finite params
(nonfin=0), kernel exception (returns cleanly), the all-or-nothing ok==0L gate
(never hit), and a leaked neqOverride (nArmed=0 everywhere).

## The three constraints any fix must satisfy

1. The augmented outer model MUST size the pool (declare odeSlotOuter /
   odeSlotOuterNode in foceiFitCpp_), or the solve is truncated as above.

2. A PERSISTENT pin cannot be used.  odeSwapPinAll() writes
   ind->neqOverride on every subject and leaves it armed; once the pool model
   is the LARGE one, solving it on the same ind runs the narrow stride while
   its dydt writes the full width -> heap corruption.  Attempt 1 (declare +
   odeSwapPinAll) crashed focei_fast AND impmap.

3. A scope must NOT be introduced at a read-only site.  getOpIndSolve(op, ind, j)
   indexes by the LIVE override, so a read must run under the same stride the
   solve was written at.  Attempt 3 armed a scope around a read-back loop and
   crashed impmap -- which does not even use the outer model, proving the
   scoping itself was at fault.

   => The unit of work is SOLVE + ITS READS inside ONE scope, per individual.

Attempt 2 (declare + per-solve OdeSwapScope in likInner0 only) also crashed:
scoping one site is not enough while other sites still solve unscoped.

## Sites that still need converting (whole solve+read regions)

`op_focei.innerNeq` must also move AFTER foceiSetup_ -- foceiSetup_ calls
rxOptionsFreeFocei(), which placement-news op_focei and resets it to 0, so
every "innerNeq > 0" gate below it is dead code today, on this branch and on
main.

- inner.cpp likInner0            innerOde + all calc_lhs reads
- inner.cpp getPopR (~1440)      innerOde is UNSCOPED today, plus a raw-slice read
- inner.cpp ~1262 read loop      reads a solve done elsewhere: must share ITS scope
- inner.cpp calcEtaHessian       rxHess2 (already scoped at 2104 -- pattern to copy)
- inner.cpp impThetaSensCollect  thetaSens (scoped at 9570)
- inner.cpp vaeOuterSolve_       rxVaeOuter
- nlm.cpp                        predOde / thetaGrad

lhs buffers: `getIndLhs(ind)` is rxode2's per-thread slice and is only
op->nlhs (the POOL model's) wide.  Every read must take its buffer from
`guard.lhs()`, which returns that slice only when the model fits it and a
private scratch otherwise.  Remaining raw uses: inner.cpp 1282 and 1451.

## Verification loop (fast signal first)

    Rscript odeSwapBaseline.R <worktree> out.rds   # impmap + focei_fast fail in ~1 min
    subRun.sh 7   # currently FAIL -- the bug above
    subRun.sh 6   # currently ok  -- control
    # then the full baseline for bit-identity on the other 6 cases

Test harness note: run test files with the small helpers sourced but WITHOUT
helper-zzz-fits.R (its ~10 cached fits pin ~21GB and make a 59s file take
hours).  `.testSeed` lives in helper-quiet.R, so load_helpers=FALSE alone makes
tests error spuriously.  See scratchpad/runFile.R.

## Attempts 4 and 5 (fresh-context session): pool sizing is necessary but not sufficient

**Attempt 4 -- declare the augmented model ONLY (no pin, no new scopes).**
Crash-free.  impmap intact (192.487097125745); every baseline case matches main
except focei_fast, which moves 2e-10 relative.  This confirms the key insight
that correctness needs the stride to be UNIFORM between a solve and its reads,
not minimal: leaving every solve at the pool width is uniform by construction,
exactly as main is uniform at 8 today.  So compaction (the pin / per-solve
scopes) is a PERFORMANCE optimization, not a correctness prerequisite -- which
is the opposite of what attempts 1-3 assumed.

But it does NOT fix the bug: subRun.sh 7 went 3 fails -> 4.  Reason: the
augmented model is solved by a SEPARATE rxode2::rxSolve (measured: "PATH
rxSolve" on all 23 calls), so sizing the FIT's pool does not touch it.  The
first-order-only integration happens inside that separate solve.

**Attempt 5 -- also register the outer model so the pooled path activates.**
`odeSwapRegister(odeSlotOuter, "outer", model["outer"], &rxVaeOuter)` after
foceiSetup_, plus op_focei.vaeOuterNeq/Nlhs.  The pooled path DID engage --
`pool=outer poolNeq=26 pooledSolveN=4`, baseline crash-free, impmap intact,
focei_fast 123.663223389668 (vs main 123.665906512565) -- but the full
test-focei-fast-grad.R SEGFAULTS (exit 139).

The baseline does not cover the multiple-endpoint model, and that is the prime
suspect: vaeOuterSolve_ was written for est="vae"'s single-endpoint M-step and
reads through offsets from am$cols.  A second modeled endpoint changes the
direction set and the lhs layout, so the pooled reader overruns.

**Next step:** make vaeOuterSolve_ multi-endpoint safe, or gate the pooled
branch on the augmented model shape it supports (single endpoint, no censoring)
and let the others keep the rxSolve route.  The registry already records
neq/nlhs per slot, so the gate is a check, not new machinery.  Verify with
odeSwapBaseline.R (fast) THEN the full test-focei-fast-grad.R (it is the file
that segfaults; the baseline alone will not catch this).

## The third swapped quantity: event-sensitivity (jump) shape -- GLOBAL

A swap has to adapt THREE things per solve, not two:

1. neq        -- ind->neqOverride (per individual, per solve)
2. lhs width  -- OdeSwapScope::lhs(); getIndLhs(ind) is only op->nlhs (POOL) wide
3. **event-sensitivity shape** -- rxode2 GLOBALS, and this is the missing one

rxode2 keeps the jump-sensitivity shape in globals (src/par_solve.cpp):

    _rxEsActive, _rxEsNState, _rxEsNParam, _rxEsNParam2, _rxEsNParam3,
    _rxEsUseCalcJac

set via the C-callable `rxode2EventSensLoad(trans, active, nState, nParam,
nParam2)` (registered in rxode2 src/init.c) or from R by
`rxEventSensLoadModel(model)` / `rxEventSensDeactivate()`.  rxode2's own docs
name this exact use case: "For downstream packages (e.g. nlmixr2est's FOCEi)
that solve a sensitivity model through a direct C++ ind_solve() loop, bypassing
rxSolve()".

rxPred has NO event sensitivities; rxInner HAS them; rxOuter has them too but a
DIFFERENT shape (different sens directions -> different nParam/nParam2).

nlmixr2est loads them once per fit, from the INNER model only
(R/focei.R:~2417, `rxEventSensLoadModel(.ret$model$inner)`).  With
foceiControl(fast=TRUE) the augmented outer model is solved through the SAME
ind_solve() loop under the inner model's ES dims, so its jumps are mis-shaped.

Do NOT "fix" this by loading the largest model instead (tried, reverted): the
injection is compartment-count guarded, so smaller models SKIP it -- loading
outer strips the INNER problem's own jump sensitivities and mis-specifies it.
One global load cannot serve both models.

### Consequence for the architecture

Because the ES shape is GLOBAL it cannot be switched per individual inside an
OpenMP region, and it cannot be interleaved.  Solves must be BATCHED BY MODEL:
load model X's ES, solve every individual for X, then swap the ES and solve the
batch for Y.  This is the same "one stride at a time, no swapping in between"
rule that governs neq and the lhs width, applied to a quantity that happens to
be process-global rather than per-ind.

Any implementation therefore needs, per registry slot: the model's ES dims
captured once, plus a scope/batch boundary that installs them outside the
parallel region and restores the previous shape after.

### Still unexplained

With the pool declared+registered for the outer model (attempt 5) and the
scope-entry rxDynLoad removed, the full test-focei-fast-grad.R STILL aborts.
rxDynLoad is NOT an ES installer -- it loads the DLL if needed and calls
rxAssignPtr(obj), i.e. repoints rxode2's global current-model pointer, which is
why calling it mid-solve corrupts the run.  The remaining abort is consistent
with the ES shape being wrong for whichever model is solved second, but that has
not been confirmed by measurement yet.

## DESIGN (as specified): batch solves by model, load each model's ES per batch

The ES shape is process-global, so it cannot be swapped per individual or inside
an OpenMP region.  The solve loop is therefore organised as BATCHES BY MODEL,
each batch bracketed by that model's ES load:

    rxEventSensLoadModel(outer)          # outer ES active
      solve every individual for OUTER   # parallel over individuals is fine:
                                         # ES globals are constant for the batch
      FLAG each individual whose outer solve is bad -- do NOT fall back inline,
      and do NOT abandon the whole gradient (drop the all-or-nothing contract)

    rxEventSensLoadModel(inner)          # swap ES: inner ES active
      for each FLAGGED individual only:
        finite-difference d(llik)/d(theta) using the INNER problem + inner ES
      restore

Why the fallback cannot be inline: the FD fallback re-optimises the subject
through the inner problem (innerOptId), which needs the INNER model's event
sensitivities.  Applying it while the outer ES is loaded mis-specifies it -- the
same failure as loading one model's ES for both.  focei (fast=TRUE) and vae
(nonMuTheta="grad") both need inner and outer to coexist across a fit, so the
only correct ordering is batch-outer -> flag -> batch-inner-FD.

### Work items

1. Registry: capture each slot's ES dims once (active/nState/nParam/nParam2/
   nParam3/useCalcJac) so a batch boundary can install them without an R call.
   rxode2 exposes `rxode2EventSensLoad` as a C-callable and
   `rxEventSensLoadModel()` / `rxEventSensDeactivate()` from R.

2. Batch boundary API in the core, e.g. odeSwapEsBatch(slot) RAII: installs that
   slot's ES on entry, restores the previous shape on exit.  MUST be constructed
   OUTSIDE any parallel region -- it is process-global state, unlike
   OdeSwapScope which is per-ind and may be used inside one.

3. Drop the all-or-nothing contract in BOTH places:
     - src/inner.cpp vaeOuterSolve_: `if (!E.ok) return R_NilValue`
       -> return per-subject ok flags with the E structures
     - R/foceiGradAnalytic.R: `any(.res$ok == 0L)` -> keep the good subjects,
       collect the flagged ones
   (The gradient is a SUM of per-subject terms, so replacing one term is
   mathematically coherent.)

4. Phase 8D2 FD fallback, run in the inner-ES batch over the flagged subjects
   only: EtaRestoreGuard + innerOptId to re-optimise and restore the subject's
   eta, per-subject Gill steps in fInd->outerThetaHf (storage already added in
   2992c46e3), shi21 central differences over np = nth + nsg + nom with the
   perturbed Omega^-1/log|Omega| precomputed before the loop.

5. Warn when the fallback engages so it reaches $runInfo: under 75 chars, no
   method-name prefix (CLAUDE.md).

6. focei.R:~2417 keeps loading the INNER model's ES for the fit as a whole; the
   outer ES is loaded only for the duration of the outer batch.

### Verification

odeSwapBaseline.R first (impmap + focei_fast fail fast), then the FULL
test-focei-fast-grad.R -- the baseline passes while that file segfaults, so the
baseline alone is not a sufficient gate.  Then subRun.sh 7 (currently FAIL) vs
subRun.sh 6 (ok).

### Scope limit: the pred fallback stays INLINE

The batch-and-defer rule applies ONLY to the outer problem.  rxPred implements
no sensitivities, so event sensitivities are out of scope for it and its solve
is unaffected by whichever ES shape is loaded.  The existing inline fallback in
likInner0 -- `fInd->doFD` -> `OdeSwapScope(odeSlotPred, ind, op)` -> predOde(),
with the guard held so the later reads of ind->solve see the same predNeq
stride -- is therefore already correct and must NOT be restructured into a
batch.  Same for the nlm pred fallback.

Only a failed OUTER (augmented) solve has to be flagged and deferred to the
inner-ES batch, because only that fallback re-enters the inner problem.

## AGREED SHAPE for installing the ES per batch (confirmed)

Install the ES shape from the C++ side, at fit setup, BEFORE the solve pool goes
live -- the same point foceiFitCpp_ already registers the peer models -- keeping
the per-slot shape in the registry and re-installing only when the ACTIVE MODEL
CHANGES BETWEEN PHASES.  Never around an individual solve, and never while the
pool solve is live.

### Why not the two placements already tried and reverted

- Inside vaeOuterSolve_ (either at the batch, or hoisted to function entry):
  CRASHES the vae grad fit.  The install goes through rxode2's R entry point,
  which reads the model's vars and can repoint rxode2's global model pointer;
  doing that while the fit's global solve is live kills the run.  Position
  within the function is irrelevant -- the solve is live for the whole call.
- In R around the gradient call (.foceiAnalyticGradFocei): out of scope for this
  refactor; the swap belongs in the C++ core.

### Constraint that forces it

The install must be an R round trip: rxode2EventSensLoad / rxode2EventSensSetActive
are registered with R_RegisterCCallable but are NOT part of rxode2's linked
function-pointer API, so calling them directly is an ABI hazard that crashes on
CRAN updates.  (Adding them to the linked API is a legitimate follow-up in
rxode2, but is explicitly NOT part of this work.)  An R round trip is only safe
where the pool solve is not live -- i.e. at setup / phase boundaries.

### Already in place (0de6841f5)

- registry captures each model's ES shape at declare time from eventSensInfo$map
- OdeSwapEsBatch RAII installs/restores a slot's shape via
  rxode2::rxEventSensLoadModel(), documented as batch-boundary-only
- odeSwapHasEs(slot) -- false for rxPred, which is why its fallback stays inline

What remains is to drive it from the fit's PHASE boundaries in foceiFitCpp_
rather than from inside the solve entry points.

## RESOLVED (f65b73a5a): the pooled fast gradient works, with two exclusions

Implemented: the augmented outer model sizes the pool (declared before
foceiSetup_) and is registered after it, so the gradient runs through
vaeOuterSolve_.  All three swapped quantities are handled -- state count (pool
sized for the augmented model, nothing compacts), ES jump shape
(OdeSwapEsBatch, per batch, outside the parallel region), and solve tolerance
(OdeSolveTolGuard).  f/lag went from segfault to passing; focei_fast is within
2e-10 of main; test-focei-fast-grad.R is at exact main parity.

### The hidden crash that invalidated attempts 1-5

odeSwapDeclare's ES capture parsed eventSensInfo$map dims with
as<IntegerVector>(map2$q) -- but map2$q is CHARACTER, so it threw inside
declare.  Combined with an unguarded Rf_eval on a missing env binding (an R
longjmp that C++ try/catch cannot catch), every earlier attempt to declare the
augmented model crashed for THIS reason, not because of the pool design.  The
capture is now a boolean; rxEventSensLoadModel() derives the dims itself.

### Why multiple endpoints are still excluded -- and where the cause actually is

NOT the gradient route.  Measured, in order:

1. Diffing both routes' per-subject E structures on the same fit:
   R/aR/AR/Rsig* IDENTICAL, only f/a/A differ, by 1e-5..1e-4.  So the symbolic
   assembly and the dvid-conditional reads are fine.
2. emax is the one component whose per-subject contributions nearly cancel (FD
   reference 63 against components of 600-34000), so a 1e-4 perturbation IS
   enough to move it by hundreds.  An earlier note here claiming the diff was
   "far too small to explain" the gap was wrong.
3. Zero subjects take the per-subject relax/retry path (relaxAfter=0,
   tolFactor=1 for all 32) -- so it is not degraded/relaxed solves.
4. DECISIVE: with pooling forced on, the SAME wrong gradient appears whether
   the pooled solve succeeds (pooledSolveN=1) or fails and falls back to
   rxSolve (pooledSolveN=0).  Two different gradient routes, one identical
   wrong answer.

=> Declaring the augmented model resizes the FIT's pool (62 states here).  The
fit's own inner solves then run inside that larger pool, the EBEs shift, and the
gradient follows.  The gradient route is not the variable.

This is the same defect as the very first one found in this work, from the other
side: an inner model solved inside a larger pool without per-individual stride
compaction.  Single-endpoint models tolerate it (focei_fast: 2e-10); emax has no
margin, so it does not.

**The real dependency is therefore the per-individual stride compaction
(Phases 9/C), not an endpoint gate.**  Compaction requires every solve AND its
reads to sit inside one OdeSwapScope -- see the three constraints above, in
particular that a scope must not be introduced at a read site, because
getOpIndSolve() indexes by the live override.  Once inner solves inside a large
pool are numerically identical to inner solves in their own pool, the endpoint
exclusion can be removed and re-measured with subRun.sh + the multi-endpoint
test.

### Also worth knowing

- rxSetSolveAtolRtol (OdeSolveTolGuard) is honoured by a FRESH rxSolve but NOT
  by an ind_solve() inside an existing pool, which reads per-subject
  tolerances.  Tightening those per subject (setIndTolFactor) makes the
  augmented solve fail outright -- tried and reverted.
- setIndSolve(ind, -1) before the augmented solve changes nothing measurable.
- delay() models cannot be pooled at all: op->delayState/delayCol are built at
  solve setup from the POOL model's stateProp and are part of the solve
  structure, not a swappable global.

## Multi-endpoint: the cause is now ISOLATED by measurement (not inferred)

An earlier note here said the cause was "the fit-wide effect of sizing the pool",
reasoned from "same wrong gradient whether the pooled solve succeeds or falls
back".  That was an inference: three things change together when the augmented
model is declared+registered (pool resize, ES batch, tolerance guard), and the
fallback case still ran the latter two.

Isolating experiment: DECLARE the augmented model but do NOT register it, so the
pool resizes to 62 while odeSwapCanPool returns odeDenyNotLoaded -- the pooled
solve, the ES batch and the tolerance guard are all definitively off (verified:
canPoolOuter=1, pooledSolveN=0).  The gradient is STILL wrong, byte-identical
(emax 569.5 vs an FD reference of 63).

=> The POOL RESIZE ALONE causes it.  The pooled solve, the ES batch and the
tolerance guard are exonerated.

### The concrete mechanism

test-focei-fast-grad.R builds its FD reference with ofvAt(), which refits using
`foceiControl(print=0, covMethod="", maxOuterIterations=0, maxInnerIterations=100)`
-- crucially WITHOUT fast=TRUE.  Those refits therefore have no augmented model
and an INNER-sized pool, while the analytic fit's pool is 62.  So the analytic
gradient is evaluated at EBEs from a 62-state-pool fit and compared against a
finite-difference reference assembled from inner-pool fits.  The two are not at
the same EBEs.

Single-endpoint models tolerate this: the EBE shift is ~2e-10 (focei_fast
matches main to that).  emax does not: its per-subject contributions cancel to 63
against components of 600-34000, so it has no margin.

### Why this makes stride compaction load bearing

The fix is to make an inner solve inside a LARGER pool numerically identical to
an inner solve in its own pool.  That is exactly per-individual stride
compaction (ind->neqOverride), and it is why the endpoint exclusion cannot be
lifted by any gate or tolerance change.  Compaction in turn requires every solve
AND its reads to share one stride -- so likInner0's stride has to span its
external reader (the read-only function near inner.cpp:1262, which reads a solve
likInner0 performed and must NOT introduce a scope of its own).  That is a
scope-lifetime change across a function boundary, not another site conversion.

## SCOPE (explicit): one rxSolve vector per fit, never rebuilt

IN scope -- this is the refactor:
  - sizing the pool's neq / nlhs for the LARGEST model, once, at setup
  - per-individual neqOverride to compact a smaller model within that pool
  - per-solve lhs buffers (rxode2's slice when it fits, a private scratch when not)

OUT of scope -- do not do this:
  - resizing or REBUILDING the rxSolve vector (the global solve structure).
    One rxSolve_ at setup; models are swapped inside it thereafter.

This indicts a route previously treated as the safe fallback:
`.foceiAnalyticSolveAll`'s `rxode2::rxSolve(am$augMod, ...)` call REBUILDS the
global solve on every gradient evaluation (its own comment says so), which is
exactly the out-of-scope operation.  The pooled path is not an optimisation over
it -- the pooled path is the supported design, and the rxSolve route is the thing
to remove, not to fall back to.

Consequence for the multi-endpoint exclusion: gating multi-endpoint back onto the
rxSolve route is NOT an acceptable resolution, because that route rebuilds the
solve vector.  The pooled path has to be made correct for it.

## Next lead: reset to the BEST eta before the outer problem

The inner problem optimises eta.  The solve state left behind is whatever eta the
optimiser tried LAST, which is not necessarily the BEST eta it found.  Anything
downstream that reuses the solve -- including the outer problem / the analytic
gradient -- then works from the wrong eta.

So before handing off to the outer problem the model must be reset to the best
eta: if it is already that eta the cost is minimal (the ODE solve is cached), and
otherwise the eta must be updated and re-solved.

This is consistent with what has been measured on the multi-endpoint case, where
four other explanations were tested and rejected:
  - per-subject relax/retry: zero subjects relaxed (relaxAfter=0, tolFactor=1)
  - solve tolerance: OdeSolveTolGuard active (1e-10 over 1e-7/1e-4), no change
  - stale-tail bad-solve scan: bounding the scan to the solved model's own neq
    (odeSwapIndBadSolveSlot) changed nothing
  - the pooled gradient route itself: declare-without-register isolates the pool
    size from the pooled solve, and the wrong gradient persists

What has NOT been tested is whether the fit's inner optimum is left on the last
eta rather than the best one, which would make the EBEs the gradient is evaluated
at inconsistent with the solve state -- and would explain why the discrepancy
lands on the one gradient component whose per-subject terms nearly cancel.
