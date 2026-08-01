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

## Best-eta + solve reset: applied at the outer solve; the fit-side half remains

rxode2 copies the initial conditions into ind->solve and carries that same vector
forward, so the BUFFER never needs clearing.  What must be reset is the solver's
POSITION state (setIndSolve(ind, -1)) -- "this subject starts a new solve".  Two
things then have to hold together before a solve:

  1. the eta installed is the BEST eta the inner problem found, not whatever eta
     its optimiser tried last; and
  2. the ODE system is reset to the start of the solve.

When the eta is unchanged this costs nothing (the solve is cached); otherwise the
eta must be updated and re-solved.

DONE: vaeOuterSolve_ now does both, in order (install eta -> setIndSolve(-1) ->
iniSubjectE -> solve).  It was the only solve site in inner.cpp that never reset
the position -- likInner0, shi21ThetaGeneral and getPopR all did.

STILL OPEN: this does NOT fix the multi-endpoint gradient, and the reason is
positional.  The declare-without-register isolation already showed the wrong
gradient appears with the pooled solve entirely off, i.e. the discrepancy is
present in the FIT's EBEs before the outer problem is ever entered.  So the
best-eta reset has to be applied at the INNER problem's EXIT -- "reset to the best
eta before it goes to the outer problem" -- in innerOpt/innerOptId, not at the
outer solve's entry.  That touches every fit and needs its own pass with a full
baseline, rather than being bolted onto this work.

## RETRACTED: "an uncompacted inner solve in a big pool destroys the fit"

An earlier revision of this file claimed a decisive measurement that pool size
alone wrecks the fit (multi-endpoint objf -10208 with poolNeq=62 vs 53695 with
poolNeq=4), and concluded that the damage scales with the pool-to-model ratio.
BOTH the experiment and the conclusion were wrong.

The experiment did not isolate pool size.  Forcing outerPoolOk also DECLARES and
REGISTERS the augmented model, routes the gradient through vaeOuterSolve_, opens
the ES batch and applies the tolerance guard.  It compared two different
estimation paths, so the objf difference cannot be attributed to the resize --
and with two endpoints it was an untested combination, so that number most
likely reflects a broken configuration rather than a property of pool sizing.

### What the passing baselines actually prove

Sizing the pool for a LARGER model and clamping neq/lhs as designed does NOT move
the objective.  This is demonstrated, not inferred:

  - impmap_scratch: poolNeq=6 (thetaSens) with a 4-state inner model, clamped via
    ind->neqOverride -> objf 192.487097125745, EXACTLY origin/main
  - focei_fast: poolNeq=26 (augmented outer) with an 8-state inner model
    -> within 2e-10 of origin/main
  - nlm_grad, focei_c1/c2, focei_full, posthoc: unchanged throughout

So a larger pool with correct clamping is numerically neutral.  Any claim that
pool size per se perturbs a fit contradicts these results.

### Status of the multi-endpoint failure: OPEN

No established finding about pool sizing.  Six gradient-side hypotheses have been
tested and rejected (per-subject relaxation, solve tolerance, the stale-tail
bad-solve scan, the pooled gradient route in isolation, the solver-position reset,
best-eta staging).  The cause is not yet identified, and the exclusion stays in
place until it is.

## CORRECTION: both exclusions are IN scope

This file previously framed the multi-endpoint and delay exclusions as settled and
out of scope.  Both framings are withdrawn.

### Multiple endpoints -- IN scope

Accuracy is not in question: evaluating the analytic gradient twice on the SAME
fit (identical thetas, identical best etas, no inner re-optimisation), once
pooled and once through rxode2::rxSolve, gives BIT-IDENTICAL gradients for a
two-endpoint model -- all 8 components, relative error 0.  test-odeswap.R pins
this via .foceiAnalyticGradViaRxSolve().

The FD comparison that originally motivated the exclusion is an invalid
instrument for the purpose: test-focei-fast-grad.R's ofvAt() refits WITHOUT
fast=TRUE, so it references unpooled fits with re-optimised etas against a pooled
analytic value, and its flat h=1e-3 divides by 2e-3, amplifying
inner-optimisation noise ~500x on a component whose per-subject terms cancel
to ~63.

### delay() -- IN scope

The stated justification was that DDE pins method="dop853"/dense, which a shared
pool fixes at setup and cannot change per solve.  That is wrong: focei ALREADY
forces the DDE configuration at the FIT level (R/focei.R, the hasDelay block
around line 1711 -- method 0, stiff2 13, dense TRUE), so a DDE fit's pool is
already built dense-dop853 and there is nothing to change per solve.

The only untested DDE-specific concern left is the delay-history column map
(op->delayState/delayCol), built at solve setup from the pool model's stateProp.
Whether the augmented model's map is compatible with the inner solves has NOT
been measured -- and the same route-agreement probe used for endpoints answers it
directly, so it should be measured rather than assumed.

### What both exclusions actually rest on now

Only the heap corruption seen when they were enabled -- and that corruption has
since been traced by valgrind to getPopR's OdeSwapScope (a83a22f96), an unrelated
and now-reverted cause: handle_evid callocs a scratch from the EFFECTIVE neq
(the override, 8 doubles) while the dydt it invokes belongs to whichever model
rxode2's globals point at (the augmented model, 26 states), overflowing three
lines after the allocation.

=> Re-test BOTH exclusions once the getPopR revert is confirmed.  If the full
file is clean with them enabled, they come out.  Use the route-agreement probe
(same fit, both routes) as the accuracy gate, not the FD comparison.


## Open, pre-existing: the augmented model can return columns nobody wrote

`test-focei-fast-grad.R`'s "fast=TRUE fit matches the finite-difference fit" fails 3
assertions (analytic gradient never consumed).  Present on `main` with the same
harness, so Phase 0's un-skipping EXPOSED it rather than causing it; the fit still
converges (objf/fixef assertions pass) -- what breaks is that `fast=TRUE` silently
buys nothing.

Characterized:
- Cumulative: blocks 0-3 alone pass, 4-6 alone pass, 0-6 together fail.
- `.foceiAnalyticGradCore` bails at the `foceiGradAllFR_` finiteness check, because
  the augmented solve returns `f`/`f1` finite and correct while EVERY later declared
  lhs (`rx_f2_*`, `rx_rvarf_`, `rx_rvar1/2_*`, `rx_rsig*`) is unwritten.  thetas,
  states, `Oi`, `dOiCube`, `tr28` are all finite -- it is not a numerical blow-up.
- Measured: the bound `calc_lhs` writes 4 of the model's 29 declared lhs, and the
  SAME model (same prefix `rx_6afe1706...`, same modelVars) resolved to two different
  `calc_lhs` addresses within one session.  So the first pointer also dangles.
- NOT the pooled path: plain `rxode2::rxSolve` reproduces it (29 columns present,
  3300 non-finite), with `pooledSolveN` at 0 for that fit.

Ruled out by measurement: rxode2 model cache (`rxClean()` before the fit changes
nothing), stale global solve atol/rtol (both guards are RAII), a wrong compaction
stride (`override=26 opNeq=26 opNlhs=29 esOk=1 canPool=0`), stale `rxInv`/omega
derivatives, prefix collision, the model being unloaded (`rxIsLoaded` is TRUE on all
20 calls), and stale on-disk artifacts (rxTempDir is per-session).

NOT reduced to a standalone rxode2 reproducer: compiling and solving 96 unrelated
models in one session leaves the first model perfect (`tools/` probe), so plain model
churn is not the trigger.  No upstream issue filed -- it needs a reproducer first.

Do NOT "detect and rebuild": a finiteness probe is unreliable here, because the
failure is UNWRITTEN memory, which can read back finite.  Tried and reverted.

Landed defence instead: `odeSwapCheckLhsWidth()` probes the bound `calc_lhs` width
once per pooled outer call and refuses to pool on a mismatch (counter
`lhsWidthMismatchN`; measured firing 15/15 on the failing fit), so the shared pool can
never consume a mismatched binding.

Review fixes on the guard itself (antigravity, 5 findings, all real):
- the probe counted the HIGHEST written slot, so a writer hitting lhs[0] and
  lhs[want-1] while skipping the middle passed; now every slot in [0, want) must be
  written;
- it probed at t=0 with the raw state, which would mis-flag an lhs inside a
  time-dependent branch; now probes at the subject's own first time, and a refusal
  warns once instead of silently dropping the pooled route;
- OdeSwapEsBatch's early return compared the ROLE, so switching between two slots
  that share odeEsOuter (thetaSens/outer/outerNode/outerCov) skipped installing the
  model about to be solved; it compares the SLOT now, while the compaction gate keeps
  using the role;
- declining to compact neutralized only a NARROWER pinned override; a wider one
  would integrate more states than the model has, so any pinned override is cleared;
- .foceiAnalyticEnsureLoaded() was removed outright: re-loading a dll in R does not
  refresh the C++ rxSolveF pointers, so it introduced a stale-pointer segfault path
  in exchange for a reload that never fired (rxIsLoaded was TRUE on all 20 calls).


## RESOLVED: the augmented model was executing another build's code

Root cause (nlmixr2/rxode2#1171): rxode2 names a model's generated .c/.so from
`.rxPre(model, modName)` = `rx_<parsed_md5>_<arch>_` -- the PARSED model text alone --
while the emitted C also depends on the event-sensitivity code, which is generated
afterwards and injected.  Two builds of one model text whose event-sensitivity code
differs therefore land on ONE .so; the later build wins, and because entry points
resolve by name (R_GetCCallable) a model object bound to the earlier one silently
executes the replacement.

Traced in one session, same prefix, same directory:

    rx_6afe1706..._.so  193856 bytes   <- eventSens="jump" build
    rx_6afe1706..._.so  167792 bytes   <- replacement, emitted from foceiFitCpp_

after which that model's calc_lhs wrote 4 of its 29 declared lhs, so f2/rvar/rvar1/
rvar2/rsig were read from slots nobody wrote and every analytic outer gradient came
back non-finite.  The same failure mode is what nlmixr2Est.R's "not provided by
package" cache-reset recovery papers over.

Fix here (independent of the rxode2 fix): the FOCEi-family generated models are built
with a role-tagged, content-hashed `modName` (role + eventSens + parsed md5), in a
build directory of our own under R's session tempdir().  BOTH parts are needed and
neither alone is sufficient -- verified:

- naming alone, artifacts in rxTempDir():        28/3 and 100/3 (still fails)
- build dir alone, no naming:                    not sufficient either
- naming + build dir outside rxTempDir():        31/0/0 and 103/0/0

What is NOT the mechanism, each ruled out by measurement: the rxode2 model cache
(rxClean() before the fit changes nothing), stale global solve tolerances, a wrong
compaction stride, stale rxInv/omega derivatives, a prefix collision between DIFFERENT
models, the model being unloaded (rxIsLoaded TRUE on all 20 calls), stale on-disk
artifacts, and the deferred thunk lacking event-sensitivity info (both colliding
builds carry it).

Known costs and follow-ups:
- tempdir() is per-session, so these artifacts are rebuilt once per session and are
  NOT shared with an opted-in rxCreateCache().  Accepted deliberately while #1171 is
  open; revert to rxTempDir() once it is fixed.
- The convention covers every model nlmixr2est generates (focei family, nlm, nls,
  saem, nlme, pruning).  rxPipeline is left alone: it re-materializes the USER's model,
  not a generated one.
- test-nlm.R's "matExp event sensitivities build HdTheta" errors, and does so
  identically at 7e679993a with every change stashed -- pre-existing, not this work.


## Phase 8D2 design, as corrected by the user (supersedes the "inline, no second pass" note)

The earlier note in this plan said the per-individual FD should run INLINE in the
augmented solve's parallel loop, with no second pass.  That is wrong and cannot work:
the augmented solve runs inside `OdeSwapEsBatch(odeSlotOuter)`, i.e. under the OUTER
model's event-sensitivity shape, while the FD needs the INNER problem.  The ES shape is
a process global and can only change at a batch boundary, so the two cannot interleave.

Corrected design:

1. **Do NOT compute the FD during the analytic gradient's ODE solve.**  A subject whose
   augmented solve fails is FLAGGED and the loop continues; failure stops being fatal to
   the whole gradient (today `inner.cpp:11090` returns R_NilValue on the first failure,
   and `R/foceiGradAnalytic.R:798` discards on `any(ok == 0L)` -- both must change to
   carry per-subject flags).

2. **Re-establish the inner problem for the flagged individual.**  A separate phase,
   outside the outer ES batch, sets the inner model's ODE/ES shape back up for that
   subject before anything is evaluated.

3. **Use `innerOpt1()` for the difference**, so the subject is re-optimized exactly the
   way the inner problem does it, rather than approximating it.

4. **Shi CENTRAL differences with an optimized step size for the fit**, and those shi
   differences keep their OWN saved `h`: they are differencing a DIFFERENT problem than
   the inner problem's `h` is tuned for, so the two step-size stores must not be shared.
   `op_focei.gouterThetaHf` / `fInd->outerThetaHf` already exists for exactly this and
   is allocated and freed but never written -- it is the store to use.

Still to confirm before coding: the eta save/restore around `innerOpt1()`
(`EtaRestoreGuard`), and whether the omega directions need the same treatment or can be
precomputed once outside the per-subject phase.


## Phase 8F: DDE is IN scope -- do not add a deny code for it

Recorded because it kept coming back.  8F's deny codes are for linCmt() and
pool-not-sized ONLY.  DDE is NOT excluded from the pooled route.

The old justification was that a delay model pins method="dop853"/dense per solve and a
shared pool cannot change that.  It does not hold: focei already forces the DDE
configuration at the FIT level (R/focei.R hasDelay block -- method 0, stiff2 13, dense
TRUE), so a delay fit's pool is built that way from the start.  The delay-history column
map (op->delayState/delayCol) is then built at setup from a pool model that already
carries the delays.

Removed accordingly: the `.dde` gate on the pooled branch of .foceiAnalyticSolveAll and
the hasDelay conjunct in focei.R's outerPoolOk.  test-dde-focei.R covers it.


### 8D2 progress and the open correctness problem

Landed and building:
- `vaeOuterSolve_` no longer discards the whole gradient on the first failed subject: it
  returns R_NilValue for that subject's entry and an integer `ok` attribute.  R records
  the flagged ids in `.foceiOuterFlagged$ids`.  Behaviour is unchanged for now -- a
  flagged subject still falls through to the rxSolve route.
- `shi21LikTheta()` (src/inner.cpp): the function shi21Central() differences.  Installs
  theta into the subject's par_ptr, re-optimizes with `innerOpt1(id, 0)`, returns
  `fInd->lik[0]` as a 1-vector.
- `foceiOuterFdInd_(ids0)`: its own phase, opens `OdeSwapEsBatch(odeSlotInner)` (the FD
  needs the INNER shape and the augmented loop runs under the OUTER one, so they cannot
  interleave), guards each subject's eta with `EtaRestoreGuard`, and shi CENTRAL
  differences each theta with the step size cached in `fInd->outerThetaHf[j]` -- the
  separate store, not the inner problem's `etahf`.
- `.foceiOuterFdForFlagged()` on the R side; init.c updated by hand (arity 1).
- The FD phase saves and restores everything `innerOpt1()` moves: the subject's eta
  (EtaRestoreGuard), its n1qn1 warm-start Hessian (zm/mode/uzm -- otherwise the fit
  would warm-start a later inner problem from a Hessian belonging to a theta it never
  visited), and the shared `op_focei.fullTheta`.

Smoke test: 12 subjects x 4 thetas, all finite, plausible magnitudes.

**OPEN -- the values do not yet agree with the analytic gradient.**  `LikInner2()` sets
`fInd->lik[likId] = -2*lik`, so `lik[0]` is the individual's -2LL contribution and
`colSums(fd)` should reproduce the analytic outer gradient's theta components with ratio
1.  Measured on theo_sd (12 subjects):

    sum of per-subject FD:  -0.605  -10.227   22.302   35.014
    analytic gradient:       0.388   -3.241  -47.584  -69.691
    ratio:                  -0.642    0.317   -2.134   -1.990

Non-constant, so this is not a scale factor.

RULED OUT: that `shi21LikTheta()` was not syncing `op_focei.fullTheta`.  It now writes
theta through to fullTheta as well as par_ptr, and the FD values are BIT-IDENTICAL --
so theta reaches the individual likelihood entirely through par_ptr and fullTheta plays
no part in it.  (The write is kept anyway: it costs nothing and keeps the two views
consistent while the subject is perturbed.)

FIXED SINCE (two real defects, both confirmed by the numbers moving):

1. **Scale.** `fInd->lik[0]` on the innerOpt1 path holds the individual's +LL, not -2LL
   (line 2424 assigns `lik`, not the `-2*lik` of line 2428).  The outer objective is
   -2LL, so `shi21LikTheta()` now returns `-2 * lik[0]`.
2. **Eta path-dependence.** Each shi evaluation re-optimized from wherever the PREVIOUS
   perturbation left the eta, so the + and - legs of a central difference were taken
   about different points.  Every evaluation now starts from one pinned reference eta
   (`_fdRefEta`, captured per subject before its loop).  This is what was wrecking the
   eta-bearing thetas: tka went from ratio -2.68 to 0.963 on this change alone.

Comparison at a common theta (maxOuterIterations=0, theo_sd, 12 subjects):

    theta    analytic     sum of FD    ratio
    tka       -1.5611      -1.6210     0.963
    tcl       -1.1936      22.3964    -0.053   <- STILL WRONG
    tv        82.6138      83.2853     0.992
    add.sd   -79.1820     -79.4173     0.997

STILL OPEN: `tcl`.  It is a real disagreement (22.4 against -1.19), not near-zero
ratio noise, and it is not simply "has an eta" -- tka also has one and now agrees.  Next
step is to difference tcl by hand, per subject, against `innerOpt1()` at +/- h to see
whether the discrepancy is in the FD or in the analytic term for that direction.

SUPERSEDED SUSPECT -- the comparison itself.  `.foceiGradAnalyticCalc(fit)`
evaluates at the FIT's final theta, while `foceiOuterFdInd_()` starts from whatever
`op_focei.fullTheta` holds after the fit returns, which is the optimizer's LAST trial
point, not necessarily the optimum.  Two gradients evaluated at different thetas will
differ by an arbitrary, per-direction amount, which is exactly the pattern above.  Fix
the harness before touching the FD again: pin both to the same theta (set fullTheta to
the fit's estimates, or take the analytic gradient at the FD's starting point), then
require ratio 1.

Do NOT wire the substitution until that comparison passes: a wrong per-subject term
would silently bias the total gradient, which is precisely the class of bug this
session spent its time on.


## NEXT: refactor the shi step search to be per inner problem (user-specified)

The per-individual d(llik)/d(theta) machinery is committed (3485a22b8) and unwired.
What must change before it can be wired in.

### The problem

The shi `h` search currently runs from the THETA loop in `foceiOuterFdInd_`: for each
subject it calls `shi21Central(shi21LikTheta, ...)`, which drives `innerOpt1()` per
perturbation.  The step is therefore chosen without the subject's inner problem being
settled, unlike `etahf`, which is optimized INSIDE the inner problem where that
subject's context is already established.  Two symptoms trace to this:

* `cl`'s steps look wrong systematically rather than sporadically;
* the finite difference is THREAD-COUNT DEPENDENT.  Proven not to be the FD loop:
  forcing that loop serial while leaving rxode2 at 4 threads reproduces the 4-thread
  answer bit for bit.  It is the inner solve running in a pool sized for the OUTER
  model under a neq override, which the relocation should remove.

### What to build

A NEW per-individual shi optimization, modelled on the `etahf` block in
`calcEtaHessian` -- same `shi21Central` machinery, same optimize-once/cache/reuse
discipline -- but SEPARATE from it.  The individual's contribution to the OVERALL
likelihood is a different function from the inner problem's conditional objective, so
it cannot share the eta routine or its step store.

* driven by `shi21LikTheta()`;
* caches into `fInd->outerThetaHf[j]` (allocated, freed, and per subject already);
* gated to the individual: the search runs where that subject's inner problem is
  established, not from an outer theta loop.

### State that must be cached across a differencing call

`innerOpt1()` moves more than the eta.  Save and restore, per subject:

* the eta (`EtaRestoreGuard`);
* the n1qn1 warm start (`zm`, `mode`, `uzm`) -- otherwise a later inner problem warm
  starts from a Hessian belonging to a theta the fit never visited;
* `op_focei.fullTheta`;
* `fInd->lik[0..2]` -- innerOpt1() WRITES it, so differencing otherwise replaces the
  subject's likelihood with one evaluated at a perturbed theta.  Added in this commit.
  The same hazard applies wherever else `innerOpt1()` is called for a side purpose
  rather than to advance the fit; those call sites are worth auditing.

### Retained unchanged

The outlier handling is orthogonal to where `h` is optimized and measured to work:
modified z-score (Iglewicz-Hoaglin, 3.5) on raw slopes against a median/MAD centre and
scale, MAD floored at sqrt(eps), analytic subjects in the reference population and
never recomputed, an outlier treated as a diagnostic ABOUT THE PARAMETER, and flagged
parameters pushed through Chartrand TV with discrepancy-principle alpha.

### Gate

`expect_identical` across 1 vs 4 threads -- the same bar Phase 10 sets for imp.  Then
wire the substitution into the gradient (C only, inside the existing gradient call) and
add the `$runInfo` warning.

### DONE -- what the refactor actually found

Built as specified: `calcOuterThetaHf()` is a new per-individual routine modelled on the
`etahf` block, driven by `shi21LikTheta()`, caching into `fInd->outerThetaHf[j]`, with
the optimize-once/reuse discipline the eta path has (search only when the cached step is
0; reuse it with a plain central difference afterwards).  Reset in `foceiOuterFinal()`
alongside `getahf`/`getahr`, so the final objective cannot inherit a step chosen at an
earlier theta.

Relocating the search then exposed TWO REAL DEFECTS, both silent, and both explaining
more than the relocation itself did:

1. **`likInner0()` short-circuits on `fInd->oldEta` alone -- THETA IS NOT IN THAT
   CHECK.**  Pinning the reference eta (which the previous commit did) makes an eta match
   the common case, so the cached `fInd->llik` from the PREVIOUS theta was returned and
   the difference across that leg was silently zero.  `setIndSolve(ind,-1)` does not help:
   the short-circuit returns before reaching the solve.  Fixed by poisoning `oldEta` with
   NA_REAL before each evaluation -- any comparison against NaN is unequal, so the
   recompute is forced exactly rather than merely made unlikely (the -42 sentinel used
   elsewhere is only "unlikely if normal").
2. **`innerOpt1()` starts from `fInd->eta`, not from `par_ptr`.**  `EtaRestoreGuard`
   only covers `par_ptr`, so the pinned reference eta was being installed in the wrong
   store: the optimizer still began wherever the previous perturbation ended.  Fixed by
   installing the reference into both, and by a new `FdInnerStateGuard` that saves and
   restores the whole inner-problem state -- `eta`, `oldEta`, `setup`, `zm/mode/uzm`,
   `lik[0..2]` -- replacing three ad-hoc partial save/restore blocks.

**CORRECTION to the root-cause attribution above.**  The thread-count dependence was NOT
"the inner solve running in a pool sized for the OUTER model under a neq override".  It
was this carried inner-optimizer state: the per-subject path depended on what the
previous evaluation left behind, and thread scheduling changed which evaluation that
was.  Nothing about the pool or the override was touched, and the dependence is gone.

Gate, theo_sd, 12 subjects (`$CLAUDE_JOB_DIR/tmp/fdgate.R`):

    1 vs 4 threads:  max |g_t1 - g_t4| = 0, steps identical      PASS
    repeat calls at the cached step: identical                   PASS

`tcl`, the direction that was wrong by a factor of ~19, is now the second most accurate.

### STILL OPEN -- the comparison harness, not the FD

Ratios against the analytic reference recorded earlier in this file:

    tka 2.018   tcl 1.016   tv 0.980   add.sd 0.996

Three directions land near 1; `tka` sits near 2.

### RESOLVED -- the whole residual error was steps that hit the shi bound

The harness suspicion above was WRONG and is retracted.  Pinning both sides in one
process (`$CLAUDE_JOB_DIR/tmp/fdref.R`) showed `fit$theta` and the FD's starting theta
are the SAME under `maxOuterIterations=0`, so there was never a theta mismatch.  The
`tka` factor of 2 was a real FD defect.

Localizing it (`fdper.R`: per-subject FD against a hand central difference at a step
inside that parameter's usable window) was decisive -- every subject agreed to ~1e-2 or
better except ONE:

    tka  subject 12:  hand  0.6427   FD -0.8965   err -1.539  (total err -1.546)
    tv   subject 12:  hand 90.1914   FD 88.4517   err -1.740  (total err -1.777)

Subject 12 carried essentially the entire error in both directions, and it is exactly
the subject whose step search terminated ON `shi21hMin` (1e-4) for those two parameters.

A step sweep of the total objective shows why, and shows the usable window is per
parameter and spans two orders of magnitude:

    ratio to analytic     h=0.3   0.03    0.01   1e-3    1e-4
    tka                   1.026  1.038   1.035  1.029   2.082
    tcl                   1.014  1.011   1.024  0.850  -0.391
    tv                    1.972  1.010   1.002  0.991   0.980
    add.sd                3.190  1.016   1.002  1.000   1.001

`tka`/`tcl` (eta-bearing) need a LARGE step and fall apart below ~1e-3; `tv`/`add.sd`
need a small one and fall apart above ~0.03.  This is a profile likelihood, so its noise
floor is the inner optimizer's convergence, not machine epsilon -- `hMin=1e-4` is inside
that noise for the eta-bearing directions.

**Fix: treat termination on the bound as the failure signal it is.**  shi returns hMin
(or hMax) exactly when its ratio test never landed in [rl, ru] and the clamp stopped it
-- the search gave up.  Those subject/parameter cells are re-differenced at the median
step of the subjects whose search DID converge for the same parameter, and the cache is
overwritten so the degenerate step is not handed back on the next call.

Detecting it on the BOUND rather than by a threshold is what makes this exact, and it is
why the earlier median-step heuristic was the right instinct with the wrong trigger.
The modified-z pass cannot see this case and testing the slope harder would not help:
subject 12's wrong `tka` slope (-0.897) is unremarkable among slopes spanning -5.0 to
+2.7.  It is not an outlier in slope, only in step -- 1e-4 against a median of 0.255.

Result (theo_sd, same process, same theta):

    ratio FD/analytic    tka 1.031   tcl 1.016   tv 1.002   add.sd 0.996   (tka was 2.018)

and the FD now matches the HAND difference to ~0.1% (tka -1.6090 vs -1.6069, tv 82.7715
vs 82.7717).  The remaining 1-3% against the analytic value is the truncation error of a
central difference on this objective -- the hand difference shows the same gap at every
step in the window -- not a defect.  Both gates still pass: bit-identical across 1 and 4
threads, and reproducible across repeat calls at a cached step.

### Audit: other innerOpt1() call sites (user-requested)

Clean, and one correction to the framing above.

* `didInnerResetPointFail` (3344) and `innerOptId` (3405) call innerOpt1() to ADVANCE
  the fit -- its writes to lik/eta are the intended output, not a side effect.  No hazard.
* The established outer FD theta gradient (7650/7665/7696) calls `innerOpt1(gid, 1)` and
  `innerOpt1(gid, 2)`, deliberately using the lower/upper likelihood slots so `lik[0]`
  is never touched -- that is what `lik[1]`/`lik[2]` are for.  No hazard, and it is the
  cleaner convention: the new path could have used `likId=2` instead of guarding
  `lik[0]`.  Not changed now because it would move results; worth doing if this code is
  ever revisited.

**CORRECTION to the likInner0 finding.**  That path ALREADY guards the eta-only
short-circuit: inner.cpp:7628 poisons `op_focei.goldEta` with the -42 sentinel before
every parameter's perturbation, for exactly this reason.  So the short-circuit is a known
hazard the codebase handles, NOT a latent package bug -- what was missing is that the new
per-individual path did not guard it.  The NA_REAL poisoning is the same mechanism made
exact rather than "unlikely if normal".  Nothing that ships today is affected by it.

### Antigravity review (commit range c0a6a0b11..eaa4f0302) -- triaged

CONFIRMED and fixed:

* **Global eta Welford accumulators.**  `innerOpt1()` (inner.cpp:3007) updates
  `op_focei.n/etaM/etaS` whenever `_innerParallel == 0`.  The TV pass ran after
  `_innerParallel.store(0)`, and `foceiIndLik_()` never set the flag at all, so both
  folded etas evaluated at PERTURBED thetas into the population statistics that drive
  the standardized-eta reset thresholds.  Fixed by a new `FdPhaseStateGuard` held for
  the whole phase: it sets `_innerParallel`, and restores `n/etaM/etaS` plus the
  `didEtaReset`/`didHessianReset`/`didEtaNudge` diagnostic flags, which are set
  regardless of `_innerParallel` and are reported to the user.
* **Per-subject state still leaking.**  `FdInnerStateGuard` now also covers
  `saveEta` (LikInner2 writes it under likId==0 -- our likId -- and foceiFinalize
  reports it), `llikObs` (per-observation conditional log-likelihoods, handed to R as
  `e["llikObs"]`), `etahf`/`etahr`/`etahh` (the INNER problem's own step caches, which
  calcEtaHessian would fill at a perturbed theta), and `stickyRecalcN2` (once it passes
  `stickyRecalcN` the subject's solver tolerance stays loosened for the whole fit).
* **Clamp detector was exact equality.**  shi's `rcur == -1` branch does `h = h*0.5/3.0`
  and `continue`, bypassing the floor clamp, so a sub-hMin return is reachable and `==`
  would miss it.  Now `<= hMin` / `>= hMax`.

REJECTED, with reasons:

* `parErrorNoEta` / `parWarnBadHess` "left set by a failed perturbation" -- these are
  written by `innerOptId()` (3403/3404/3412), not by `innerOpt1()`.  This path calls
  `innerOpt1()` directly and never reaches them.
* "Repaired steps re-enter the median set on later calls, skewing it" -- by design: the
  repair overwrites the cache precisely so the degenerate step is not reused, and the
  replacement is itself the converged median.  Not a defect.
* "If every subject clamps, nothing is repaired" -- deliberate and documented: with no
  converged step for that parameter there is nothing to repair TO, and inventing one
  would be worse than leaving the value flagged.

NOTED, not changed: shi can also return exactly hMin having ACCEPTED it (ratio test in
range at the floor), which the detector cannot distinguish from giving up there.  Both
are treated as suspect.  That conflation is deliberate -- a step at or below the floor is
untrustworthy for this objective either way -- but it is a conflation, not a precise
test.  Distinguishing them would mean returning a termination reason from shi21Central,
which the eta and nlm paths also call.

## Phase 8E: pull the analytic gradient assembly into C++ (user-specified)

### Why

The FOCEI analytic gradient currently makes a C++ -> R -> C++ round trip:

    vaeOuterSolve_        (C++)  builds per-subject E structs, wraps them to R lists
    .foceiAnalyticGradCore (R)   stacks them into aB/AB/aRB/ARB/RsigB/RsigDirB/... 
    foceiGradAllFR_       (C++)  takes those 25 arguments and computes the gradient

Two problems, both the user's stated reason for the change:

1. **Pool state.**  Returning to R between the solve and the assembly gives R the
   chance to disturb the shared pool -- exactly the class of failure this whole branch
   exists to remove.  The solve and everything that reads it should be one C++ region.
2. **Cost.**  `AB` and `ARB` are `totObs x ndir x ndir`.  They are materialized twice
   (wrapped out of C++, then read back in) for no reason: the data never leaves the
   process and the kernel that consumes them is C++ anyway.

The FD substitution for flagged subjects belongs in the same place -- R currently
zero-fills a flagged subject's E purely so the R-side stacking keeps working
(foceiGradAnalytic.R:564), which is a workaround for the round trip, not a design.

### Feasibility (checked)

* The per-subject `VaeOuterE` structs already exist in C++ (`std::vector<VaeOuterE> Es`)
  and hold f/a/A/R/aR/AR/Rsig/RsigDir/trans.
* The DV transform R applies while stacking is ALREADY C++: `.foceiAnalyticTbsY()` is a
  thin `.Call(_nlmixr2est_powerD, ...)`, and `_powerD` is used directly in censResid.cpp.
  So `yB`, `limB` and the lambda dvSens/jacobian terms can all be formed C++-side.
* Everything else the kernel needs is small and per-call (Oi, dOiEst, tr28, dirTh,
  sigCol, neta/nth/nsg/nom, censOpt) or O(totObs) scalars (DV, CENS, LIMIT) -- cheap to
  pass.  Only the ndir^2 cubes are worth keeping in C++, and those are exactly the ones
  that would stop crossing.

### Steps (each committable and separately verifiable)

1a. Refactor `vaeOuterSolve_`'s body into a helper returning `std::vector<VaeOuterE>` +
    the ok flags, leaving `vaeOuterSolve_` as the R-facing wrapper over it.  No behavior
    change; gate is bit-identical gradients.
1b. New entry `foceiAnalyticGradPooled_`: calls the helper, closes the outer ES batch,
    runs the per-individual FD for flagged subjects, stacks into the arma matrices in
    C++, runs the existing per-subject kernel, SUBSTITUTES the FD row into `gmat.col(i)`
    for each flagged subject, and returns `list(g, etaP)`.
1c. `.foceiAnalyticGradCore` calls it and skips the R stacking block when it succeeds;
    the rxSolve route stays as the fallback.  The zero-fill workaround then goes away.

Remember `src/init.c` is MANUAL: a new `[[Rcpp::export]]` needs BOTH
`Rcpp::compileAttributes(".")` AND a hand-edited prototype + `callMethods[]` arity.

Gate: `g` bit-identical to the current path on theo_sd, the multiple-endpoint model and
a censored model, plus the existing 1-vs-4-thread identity.

### Phase 8E open decision: which solve tolerance the FD runs at

`OdeSolveTolGuard _tolGuard(tol)` in `vaeOuterSolve_` was introduced on THIS branch
(f65b73a5a), and its scope extends over the per-individual FD block that was added
later.  So the FD's `innerOpt1()` re-optimizations run at the ANALYTIC solve tolerance
(~1e-10), not the fit's.

That is a scope accident, not a decision, and it is not what was measured: every
accuracy number recorded above came from calling `foceiOuterFdInd_` directly from R,
where no tol guard is active.  The shipped path would differ from the gated path.

Decide it explicitly during the port, and measure both:

* the gradient should be of the objective the OPTIMIZER minimizes, which argues for the
  fit's tolerance;
* but the analytic gradient it is being substituted into is computed at 1e-10, which
  argues for consistency with its neighbours in the same sum.

Whichever is chosen, the FD must be gated at the tolerance it will actually run at.

### Phase 8E: the analytic gradient now solves at the FIT's tolerance

Done, both routes.  There is no longer a separate "analytic gradient" tolerance:
`.foceiGradSolveTol(ui)` returns the fit's own `rxControl` atol/rtol and is used at the
gradient call sites, while `.foceiAnalyticSolveTol()` (and `covSolveTol`) stays where it
belongs, on the COVARIANCE path.  The pooled route resets to the fit's values in C++
(`OdeFitTolGuard`), so a tolerance shift earlier in a fit cannot leak into the gradient.

**A vacuous verification, corrected.**  The first check of the multiple-endpoint model
reported "PASS, unchanged" -- because that model does not use the pooled route at all.
Probing `odeSwapInfo_()$pooledSolveN` across a gradient call settles it:

    theo_sd (1 endpoint)   pooledSolveN 0 -> 1   POOLED
    warfarin (2 endpoints) pooledSolveN 1 -> 1   rxSolve fallback

The two-endpoint augmented model is larger than the pool, so it takes the fallback, which
still carried 1e-10.  Changing only the pooled route left that case untouched and the
identical numbers looked like evidence when they were the absence of it.  Do not read an
unchanged result as a passing test without first checking the code ran.

Measured once BOTH routes were on the fit's tolerance:

    warfarin PKPD, analytic vs central difference of the OBJECTIVE (the test's own
    reference, itself taken at the fit's tolerance):

      max rel diff   0.008996 (tight 1e-10)  ->  0.009063 (fit tolerance)   [limit 0.01]
      emax           3e-5 either way

So the tight tolerance was not buying accuracy here: the change is neutral, and the
`emax 63 vs 569` blow-up recorded for f65b73a5a is NOT reproduced by this test.  The
worst term (`tcl`, 0.9%) sits near the 1% limit both before and after, so that tightness
is pre-existing and comes from somewhere else -- most likely the reference itself, an
h=1e-3 difference of a refit objective.

theo_sd is unchanged to three decimals in the FD/analytic ratios.

### Phase 8E steps 1b/1c: written, exact, and GATED OFF

`foceiAnalyticGradPooled_` does the whole FOCEI (f,R) gradient in one C++ region: solves
the augmented model in the shared pool, finite-differences the failed subjects, reads DV
/ CENS / LIMIT straight from the individual (`getIndDv` / `getIndCens` / `getIndLimit`,
the same way the inner problem reads them), applies the transform with `_powerD` /
`_powerDLambda` / `_powerDL`, stacks, runs `foceiGradSubjectFR_` (de-static'd and
declared in the new `src/foceiGrad.h`) and substitutes the FD rows.  R dispatches to it
before `.foceiAnalyticSolveAll` and returns early; anything unsupported falls through.

Verified: `max abs diff 0` against the staged path on theo_sd -- bit-identical.

**BUT: it leaves R's protect stack unbalanced by -10, so it is OFF unless
`NLMIXR2EST_GRAD_POOLED` is set.**  More UNPROTECTs than PROTECTs risks premature GC of
live objects; results being exact today does not make that safe.

What is known, so the next attempt does not redo it:

* Argument marshalling alone is clean -- returning immediately on entry gives no
  imbalance, so the 15 arguments (including the by-value arma types) are not the cause.
* The imbalance is already present by the end of the solve, i.e. within
  `OdeSwapEsBatch` + `OdeFitTolGuard` + the `cols` gates + `vaeOuterSolveFill`.
* `vaeOuterSolve_` performs that SAME sequence and is clean (checked directly by forcing
  the old route).  So it is not any one of those calls in isolation; it is something
  about this entry point.  That contradiction is the thing to explain first.
* Beware: `MODE=staged` (`.odeSwapNoPool$on <- TRUE`) does NOT exercise `vaeOuterSolve_`
  at all, so a clean staged run proves nothing about it.  Force the old pooled route
  instead.

Also worth recording, both self-inflicted:

* `nlmixr2est:::.odeSwapNoPool$on <- TRUE` is not a valid assignment (`:::<-` does not
  exist); it fails with "object 'nlmixr2est' not found" and produces spurious stack
  imbalance warnings of its own.  Bind the env first, then assign into it.
* Running `pkgbuild::compile_dll()` in this worktree while a test's `load_all()` is
  compiling there destroys both -- the gate run died with no result and no error.
  Serialize builds against test runs.

### Phase 8E is wired at the WRONG LAYER -- assessment

Checked at the user's prompting.  The 8E work is real but it optimized the middle of a
path that should not exist.  Per OUTER GRADIENT EVALUATION the live sequence is
C++ -> R -> C++ -> R -> C++:

`analyticOuterGrad` (inner.cpp:4551), all before it ever reaches the gradient:
  * `foceiEtas()` builds an R list/data.frame of EVERY subject's etas + OFV -> `etaObf`
  * `getOmega()` -> R matrix -> `omega`;  `.gradTheta` -> R vector
  * `Function _agf = ...".foceiCalcGradAnalytic"` -> **calls into R**

R `.foceiCalcGradAnalytic` -> `.foceiAnalyticGradFocei` -> `.foceiAnalyticGradCore`:
  * re-derives the direction set / Oi / dOiEst / tr28 EVERY call
  * **`.Call`s back down** (foceiAnalyticGradPooled_, or the staged kernel)
  * stashes etaP on the env as `.foceiGradEtaP` (neta x npars x nsub doubles)

back in C++:
  * `as<NumericVector>(_res)` for g
  * `loadAnalyticEtaP()` (inner.cpp:4515) reads `.foceiGradEtaP` BACK OUT of the env and
    copies it into `op_focei.getaP`

So etaP is computed in C++, wrapped to R, stashed on an environment, read back out, and
copied into a C++ array -- a full round trip to no purpose.  The etas make the mirror
trip: computed in C++, wrapped into a data.frame, parsed back into a matrix in R, sent
down again.

8E removed the INNER round trip (E structs -> R -> stacked -> back down) and left the
outer one, then added one more R object to it (the etaP cube in `List::create`).  That is
why the state is wrong, and it is the same class of error as the earlier foceiGradAllFR_
misfire: correct code at the wrong layer.

### What "right" looks like

* `analyticOuterGrad` calls the C++ gradient DIRECTLY -- no `Rcpp::Function`, no env
  writes, no `.Call` back down.  The only `.Call` on a fit's hot path should be the inner
  problem's.
* Per-fit setup (dirs, cols, Oi, dOiEst, tr28, dirTh, sigCol, lamDir) cached in C++ once,
  not rebuilt in R on every gradient call.
* `g` written into the caller's `double *g`; etaP written straight into `op_focei.getaP`
  with the `dUnscaleParDx` scaling, never becoming an R object.  `foceiEtas`/`omega`/
  `.gradTheta` env writes drop out entirely for the analytic path.
* Keep a thin R-facing `.Call` wrapper for TESTS only -- that is the one place building
  R objects is worth it.
* Where SEXPs are genuinely needed, use the existing `rxProtect` RAII guard
  (src/rxProtect.h) rather than hand-balanced PROTECT/UNPROTECT.

This very likely also disposes of the -10 protect imbalance: on the corrected path the
live gradient creates NO R objects, so the imbalance has nowhere to come from.  Chasing
it inside a design that is being removed is not worth the cycles -- fix the layer first,
then re-check.

### Phase 8E FIXED AT THE RIGHT LAYER: analyticOuterGrad calls C++ directly

`analyticOuterGrad` now tries `analyticOuterGradDirect()` first, and that path touches R
nowhere.  What it replaced, per GRADIENT EVALUATION: building `etaObf` (a data.frame of
every subject's etas via `foceiEtas`), `omega` and `.gradTheta` as R objects; an
`Rcpp::Function` call into R; R re-deriving the whole setup; a `.Call` back down; and
finally C++ reading etaP back OUT of the fit env into `op_focei.getaP`.

How each per-call input is obtained now, without R:

* theta      -- `op_focei.fullTheta`
* EBEs       -- `fInd->saveEta` (NOT `fInd->eta`: saveEta is the eta at which lik[0] was
                computed, and it is what `foceiEtas` fed the R route, so the two agree by
                construction)
* Omega      -- `getOmegaInv()`
* dOmega^-1 / tr.28 -- new `getDOmegaInvL()` / `getTr28V()` macros on the SAME `_rxInv`
                handle the inner problem already uses.  Nothing is recomputed and nothing
                is shuffled; this was the user's point.
* g          -- written into the caller's `double *g` (x dUnscaleParDx)
* etaP       -- written STRAIGHT into `op_focei.getaP` on the scaled parameterization,
                never becoming an R object

Per-fit constants (lhs column map, dirTh/sigCol/lamDir, dimensions, censOption) are
computed once by `.foceiGradPooledSetup()` before `foceiFitCpp_` and copied into a POD
(`FoceiGradPooledSetup`) by `loadGradPooledSetup()` when the outer optimizer starts.  The
POD holds plain vectors, so no R object stays alive across the fit, and
`vaeOuterSolveFill` was moved onto it (plus `std::vector<double>`/`arma::mat` instead of
NumericVector/NumericMatrix) so the solve touches no SEXP either.

Measured, theo_sd, full fits:

    objf fast=FALSE 130.2898192   fast=TRUE 130.2921878   rel diff 1.8e-5
    analytic gradient iterations 16, finite-difference iterations 0
    extra: "grad: analytic"

Comfortably inside test-focei-fast-grad.R's 0.02 objf / 1e-2 fixef tolerances.

**The -10 protect-stack imbalance is gone** -- zero warnings on the direct path.  That
confirms the earlier read: it was a symptom of building R objects in a nested
C++ -> R -> C++ path, not a bug to hunt in isolation.  Fixing the layer removed it, which
is why chasing it directly was wasted effort.

The R-mediated route remains as the fallback for everything the direct path does not
cover (FOCE, AGQ, ll(), mixtures, models whose augmented form outgrows the pool), and
`foceiAnalyticGradPooled_` stays as the R-callable wrapper for tests.

### Extending the direct C++ gradient to FOCE / AGQ / ll() -- scope

User direction: FOCE, AGQ and ll() should all run through the direct path; only mixtures
stay out (they need probability-weighted per-component contributions plus the derivative
of the weights -- solvable, not solved).  Then the R route is deleted.

Surveyed.  The three are NOT equal in cost.

**FOCE** -- kernel `foceiGradSubjectFoceFR_` / `foceiGradAllFoceFR_` already exists in
C++.  The assembly is the same shape as FOCEI with different blocks: `aRe`/`aRc` instead
of `aR`/`AR` (no AR cube), `R0`/`R0sig` instead of `R`/`Rsig`/`RsigDir`.  Two extra
pieces, and the second is real work:
  * a SECOND augmented solve at eta=0 (`E0all`) when `foceType==0` and `ef$dependsF0` --
    cheap, it is `vaeOuterSolveFill` with a zero ebes matrix;
  * `.foceiAnalyticFoceEbeBatch` -- an EBE RE-SOLVE at the frozen R0.  This is a Newton
    solve currently done in R and it has no C++ equivalent.  It is a genuine port, not
    wiring.  (`fInd->saveEta` is NOT a substitute: those are the inner problem's etas,
    and the FOCE gradient needs the EBEs of the frozen-R0 objective.)

**AGQ** -- kernel `foceiGradSubjectAgqFR_` / `foceiGradAllAgqFR_` exists.  Needs the
quadrature grid built from Ht's Cholesky FACTOR (`etaSolve[i,] + sqrt(2) * GinvL %*% x`)
and one augmented solve per node, then the same stacking.  Moderate, self-contained.

**ll()** -- `.foceiAnalyticGradCoreLL`, ~93 lines of R, and there is NO C++ kernel for it
at all.  It also carries the `fd2` route (directional central FD of the analytic 2nd-order
A).  This is the largest of the three.

### Sequencing, and why the R route cannot go first

Deleting the R route before these land would silently strip FOCE/AGQ/ll fits of their
analytic gradient (they would fall to finite differences), which is a behaviour
regression, not a cleanup.  Order: FOCE, AGQ, ll(), THEN delete the R path and the
`.foceiCalcGradAnalytic` hook with it.

Per the user's other instruction: when the analytic outer gradient cannot be computed it
must not be applied at all -- so the direct path returning false should mean the fit uses
its ordinary gradient, with no second, slower analytic attempt through R.
`foceiControl(fast=TRUE)` replaces the finite-difference outer gradient with an analytic
one (Almquist 2015).  As of `222bbb38c` the **FOCEI** case runs entirely in C++:
`analyticOuterGrad()` -> `analyticOuterGradDirect()` -> `gradPooledCore()` solves the
augmented model in the shared ODE pool, stacks the sensitivities and calls
`foceiGradSubjectFR_`, without touching R.  Measured on theo_sd: 16 analytic gradient
iterations, 0 finite-difference, objf within 1.8e-5 of `fast=FALSE`, and the -10
protect-stack imbalance the previous C++ -> R -> C++ design produced is gone.

Three shapes still fall back to an R route (`.foceiCalcGradAnalytic`, reached through
`Rcpp::Function` at `src/inner.cpp:4663`):

* **FOCE** (`interaction=0`), both `foceType` 0 "nonmem" and 1 "foce+"
* **AGQ** (`nAGQ > 1`)
* **general likelihood** (`ll()` endpoints)

That route rebuilds `etaObf`/`omega`/`.gradTheta` as R objects on **every gradient
evaluation**, re-derives the per-fit setup in R, `.Call`s back down, and reads `etaP`
back out of the fit environment -- the exact round trip Phase 8E removed for FOCEI.  It
also lets R run between the augmented solve and the assembly, where it can disturb the
shared solve pool.

Intended outcome: all four shapes compute the gradient in C++ with no R on the hot path;
the R gradient implementation is then deleted outright and the existing gradient tests
are re-pointed at the C++ code the fit actually runs.  Mixtures stay out of scope (they
need probability-weighted per-component contributions plus the derivative of the
weights -- solvable, not solved) and keep the ordinary finite-difference gradient.  When
the analytic gradient cannot be computed it is not applied at all: no second, slower
analytic attempt through R.

## Two decisions taken up front

1. **The R gradient implementation goes away completely**, including the post-fit
   `.foceiGradAnalyticCalc(fit)` entry.  That entry is currently what ~15 assertions in
   `test-focei-fast-grad.R`, `test-agq-fast-grad.R`, `test-subFinal.R` and
   `test-odeswap.R` call -- so today **none** of them exercise the direct C++ path.  They
   get re-pointed at a thin test-only wrapper over the same C++ code the fit runs, so the
   tests verify shipping code instead of a parallel R implementation that can drift.
2. **The FOCE Newton EBE re-solve may improve on the R version** rather than reproducing
   it bit-for-bit (per-subject early exit from the batched solve, better tolerances).
   FOCE fit results may shift slightly; the gradient-vs-central-difference assertions are
   the correctness bar, not equality with today's numbers.

## What already exists (do not rebuild)

| Piece | Where |
| --- | --- |
| FOCEI direct path | `src/inner.cpp` `analyticOuterGradDirect` (12322), `gradPooledCore` (12162) |
| Augmented pooled solve | `src/inner.cpp` `vaeOuterSolveFill` (11928) -- SEXP-free; takes the eta matrix as an argument, so a zero matrix gives the FOCE eta=0 solve |
| Per-fit setup POD | `src/inner.cpp` `FoceiGradPooledSetup` (4565), `loadGradPooledSetup` (4612) |
| Per-subject FD fallback | `src/inner.cpp` `foceiOuterFdInd_`, `calcOuterThetaHf`, `fdCentralAt`, `shi21LikTheta` |
| Pool slots for all three extra models | `src/odeSwap.h`: `odeSlotOuter`, `odeSlotOuterNode` (AGQ order-1), `odeSlotHess2` (ll d2) |
| ES-shape batch boundary / per-ind scope | `src/odeSwap.h`: `OdeSwapEsBatch`, `OdeSwapScope` |
| FOCE gradient kernel | `src/foceiGrad.cpp:471` `foceiGradSubjectFoceFR_` (**`static`**) |
| AGQ gradient kernel | `src/foceiGrad.cpp:1043` `foceiGradSubjectAgqFR_` (**`static`**, non-throwing, reports via `ok_out`) |
| Censored partials | `src/censEst.h` `censNormalPartials()` -- header-inline, already called from C++ |
| **AGQ quadrature grid** | already in C++: `op_focei.aqx` / `op_focei.aqw` (`_aqn x neta`), `_aqn`, loaded by `setupAq1_` (inner.cpp:9417) |
| **ll() per-subject Hessian reader** | already in C++: `calcEtaHessian`'s hess2 branch (inner.cpp:2278-2315) reads `rx__d2pred_` into `H` |
| **Method-agnostic per-subject FD** | `shi21LikTheta` re-optimizes eta at the perturbed theta via `innerOpt1`, i.e. it differences *the objective's own inner problem* -- failed FOCE/AGQ/ll subjects already fall back correctly with no change |

Both gradient kernels take **pure arma / plain scalar arguments only** -- no SEXP in any
parameter list -- so calling them from C++ constructs zero R objects.  They are `static`
only because they never had a C++ caller.

Two R-side complications evaporate in the pooled path and must NOT be ported:

* the AGQ **node-replicated event data** (`.foceiAgqRepData`, `.foceiAgqRepCache`) exists
  only because `rxSolve` needs pseudo-subjects; the pool solves one individual at a time,
  so nodes become a loop.
* the `innerHess2` **`ETA[k]`/`THETA[k]` bracket-naming** gotcha exists only because R
  sets parameters by name; C++ writes `par_ptr` positionally through `op_focei.etaTrans`
  / `op_focei.thetaTrans`.

## Phase 0 -- make the kernels callable

* `src/foceiGrad.cpp`: drop `static` from `foceiGradSubjectFoceFR_` and
  `foceiGradSubjectAgqFR_`.
* `src/foceiGrad.h`: declare both, beside the existing `foceiGradSubjectFR_`.

No `init.c` / `compileAttributes` churn -- these are not `[[Rcpp::export]]`.

When writing driver loops: `foceiGradSubjectFR_` and `foceiGradSubjectFoceFR_` use
**throwing** `arma::inv()` and the caller owns the `try/catch(...)` (`Makevars` sets
`ARMA_DONT_USE_OPENMP` but not `ARMA_DONT_USE_EXCEPTIONS`, so a throw escaping an OpenMP
structured block is `std::terminate`).  The AGQ kernel uses the non-throwing forms.

## Phase 1 -- generalize the augmented solve

`vaeOuterSolveFill` is hardcoded to `odeSlotOuter` / the global `rxVaeOuter` and always
reads the 2nd-order blocks.  All three new routes need it against a different slot, model
and column map.  Proposed signature (`src/inner.cpp`):

```cpp
// what to read back from the solve
enum OuterFillWhat {
  outerFillFull = 0,   // f,a,A,R,aR,AR,Rsig,RsigDir,Rsig2,trans  (odeSlotOuter)
  outerFillOrder1,     // f,a,R,aR,Rsig                            (odeSlotOuterNode, AGQ)
  outerFillEtaHess     // eta-eta 2nd derivative block only         (odeSlotHess2, ll)
};

static void outerSolveFill(int slot, rxSolveF *fns,
                           const std::vector<double> &thVals, const arma::mat &ebes,
                           const FoceiGradPooledSetup::Cols &C, int what,
                           int cores, rx_solving_options *op, int nsub, int neta,
                           std::vector<VaeOuterE> &Es);
```

`vaeOuterSolveFill(...)` becomes `outerSolveFill(odeSlotOuter, &rxVaeOuter, ...,
outerFillFull, ...)`.  Split the lhs column-map fields of `FoceiGradPooledSetup` into a
nested `Cols` struct so the POD can hold three (`cols`, `colsNode`, `colsHess2`).

The caller keeps ownership of `OdeSwapEsBatch` and `OdeFitTolGuard`, as it does today.
`odeSwapEsModelForSlot` already maps `odeSlotOuterNode` to the same `odeEsOuter` role as
`odeSlotOuter`, so an AGQ node loop stays inside one batch; `odeSlotHess2` has its own
role and needs its own.

**Commit and run the FOCEI tests here** -- this refactor must be provably a no-op for the
working path before anything is built on it.

## Phase 2 -- widen the per-fit setup, and add the test-only C++ entry

### 2a. Setup

`.foceiGradPooledSetup()` (`R/foceiGradAnalytic.R:530`) returns `NULL` for
`interaction != 1L`, `nAGQ > 1L` and `ef$isLL` (550-552) and requires `cols$hasR` (555).
C++ has a sixth, structural exclusion: `loadGradPooledSetup`'s `if (!G.hasR) return;`
(inner.cpp:4628), which must relax in step.  Note `.foceiAnalyticGradSetup` reports
**`interaction = 0L` for ll() as well**, so `interaction` alone cannot separate FOCE from
ll -- carry `isLL` explicitly and test it first.

Replace those gates with per-shape assembly:

* always carry `interaction`, `foceType`, `nAGQ`, `isLL`, `dependsF0`, `canVanish`,
  `censOpt` (`dependsF0`/`canVanish` come from `.foceiAnalyticErrFull(ui)`, per-fit
  constants)
* FOCE: no new columns -- the eta=0 solve reuses the same augmented model
* AGQ: add `colsNode` from `.vaeOuterCols(amNode)`, `amNode` preferring
  `ui$foceiModel$outerNode` + `outerNodeMeta`, else
  `.foceiAnalyticAugModelDirs(ui, dirs, order = 1L)`; preserve the current
  "fall back to the order-2 model" behaviour by emitting `colsNode = cols`
* ll(): drop the `hasR` requirement; add `thPos` (each structural theta's `ntheta`
  position, needed for the non-mu perturbation) and the hess2 column info
* keep every existing scope gate (mixtures, `linCmt`, IOV, bounded transforms, fixed
  mu-thetas, `cholSEOpt`, finite `agqLow`/`agqHi`, FOCE + `censOption="laplace"`,
  AGQ + censoring/lambda, ll + censoring) -- each becomes a per-shape `ok` flag rather
  than a blanket `NULL`.

Mirror the new fields in the C++ POD and in `loadGradPooledSetup`.  This stays a
**once-per-fit** read.

Everything else the routes need is already reachable in C++ and must come from there, not
R: theta (`op_focei.fullTheta`), EBEs (`fInd->saveEta`), Omega / `d.omegaInv` / `tr.28`
(`_rxInv` via `getOmegaInv()` / `getDOmegaInvL()` / `getTr28V()`), DV / CENS / LIMIT
(`getIndDv` / `getIndCens` / `getIndLimit`), the DV transform (`_powerD`, taken from the
augmented model's `E.trans`, **not** `getIndLambda`), the thread count (`getOpCores`), and
the AGQ grid (`op_focei.aqx/aqw`, `_aqn`).

### 2b. The test entry point, and the mechanism counter

Both are needed *before* Phase 3, not after, so each shape can migrate its own tests as it
lands.

* A test-only `[[Rcpp::export]]` that computes the gradient for a fit through
  `analyticOuterGradDirect`/`gradPooledCore` and returns `list(g, etaP)`.
  `foceiAnalyticGradPooled_` (inner.cpp:12440) is the existing shape to follow.  This is
  the one place where building R objects is worth it.  Needs `Rcpp::compileAttributes(".")`
  **and** a hand edit to `src/init.c` (prototype + `callMethods[]` arity).
* `op_focei.nAnalyticGradDirect`, incremented on a successful `analyticOuterGradDirect`.
  `nAnalyticGrad` / `nFDGradFast` (inner.cpp:349-350) already drive the `"; grad: analytic"`
  suffix in the fit's `extra`, but they count the R route and the direct route alike, so
  during Phases 3-5 a silent fallback to R would look exactly like success.

## Phase 3 -- FOCE

A FOCE branch in `gradPooledCore` (`G.interaction == 0 && !G.isLL`):

1. **eta=0 solve** when `foceType == 0 && dependsF0`: `outerSolveFill` with a zero `ebes`
   -> `E0all`.  Same slot and ES batch as the eta-hat solve.
2. **Newton EBE re-solve at the frozen R0** -- the port of `.foceiAnalyticFoceEbeBatch`
   (`R/foceiCovAnalytic.R:2304-2337`).  Per subject: `S = Oi*eta + sum(q0 * a)`,
   `Hf = Oi + sum(q1 * a a' + q0 * A)`, with `q0 = -(y-f)/R0`, `q1 = 1/R0` replaced on
   censored rows by `censNormalPartials(..., 2)` columns 1 and 3; iterate
   `eta -= solve(Hf, S)`.  Free to improve on the R structure -- notably by dropping
   converged subjects out of the batched solve instead of re-solving everyone each
   iteration.  Non-convergence returns false (ordinary gradient).  Put it in a helper
   `foceEbeNewton(...)`; the per-subject algebra parallelizes over subjects while the
   solve stays batched, so the ES batch remains outside.
   `fInd->saveEta` is **not** a substitute -- those are the inner problem's etas, not the
   frozen-R0 objective's.
3. **Stack + kernel**: the FOCEI loop with `aRe`/`aRc`/`R0`/`R0sig` in place of
   `aR`/`AR`/`R`/`Rsig`/`RsigDir` and no `AR` cube; `fp = (foceType == 1 || no E0)`; call
   `foceiGradSubjectFoceFR_`.

Then re-point the FOCE assertions in `test-focei-fast-grad.R` at the Phase-2b entry and
add the `nAnalyticGradDirect > 0` assertion for a FOCE fit.

## Phase 4 -- AGQ

Branch on `G.nAGQ > 1`:

1. eta-hat solve as FOCEI (full order-2 blocks).
2. Per subject `Ht = Oi + sum(a a'/R + 0.5 aR aR'/R^2)`, `chol`, **require
   `rcond(Ht) >= 1e-10`** (arma `is_sympd()` and R `chol()` disagree near the PD boundary;
   a mismatch means the objective took the `nmNearPD`/`cholSE` branch, a different
   non-smooth function), `GinvL = inv(trimatu(chol))`.
3. Node loop over the `_aqn` rows of `op_focei.aqx`: eta = `ehat + sqrt(2) * GinvL * x_k`
   -- the `sqrt(2)` must match the kernel's `etaCur` or the node sensitivities sit at the
   wrong eta -- with one
   `outerSolveFill(odeSlotOuterNode, ..., outerFillOrder1, ...)` per node over all
   subjects.  Drop the R chunking (`.maxPs = 2048`): it exists only to bound pseudo-subject
   count in one `rxSolve`.
4. Assemble node-major buffers (node `k` at rows `k*nobs`), call
   `foceiGradSubjectAgqFR_`; any subject with `ok_out == false` fails the whole gradient.

Out of scope, already gated in R and to stay gated: censoring, estimated DV-transform
lambda, finite `agqLow`/`agqHi`, `cholSEOpt`.

Then re-point `test-agq-fast-grad.R` and add its mechanism assertion.

## Phase 5 -- ll()

No new kernel is strictly required, but for symmetry and unit-testability put the
per-subject algebra in `src/foceiGrad.cpp` as `foceiGradSubjectLL_` with the same
out-param shape as the other three.  Port `.foceiAnalyticGradCoreLL`
(`R/foceiGradAnalytic.R:1101-1193`):

1. eta-hat solve on `odeSlotOuter` for `a` and `A`; per subject `H = Oi + (-sum A)`, PD
   check on the min eigenvalue, `Hi = inv(H)`.
2. Build `2*(nth+nom)` perturbed configurations, `hFD = 1e-4`.  The direction is **not** a
   coordinate axis: for theta `p` with `s = dirTh[p]`, `sh = Hi * d2` where
   `d2[l] = sum(A[,l,s])`, plus `sh[s] += 1` when `s <= neta` (mu-referenced); for omega
   `k`, `etaP = Hi * (-dOi * eta)`.  Mu-referenced thetas hold the theta value and move
   only the eta (perturbing an ODE-entering theta goes non-finite); non-mu thetas move
   `th[thPos[p]] +/- hFD`.
3. Solve each configuration on `odeSlotHess2` reading only the eta-eta block
   (`outerFillEtaHess`), falling back to `odeSlotOuter` when hess2 is unavailable.  Its own
   ES role, so its own `OdeSwapEsBatch`.
4. Assemble: `g[p] = -2*dl[p] + sum(Hi % dH)`;
   `g[nth+k] = dq[k] - 2*nsub*tr28[k] + sum(Hi % (dOi + dH))`.  `etaP` is NULL for ll(), so
   set `op_focei.etaPValid = 0` and skip the Eq-48 warm start, matching today.

Out of scope: censored observations (they enter as `-logPhi`, which the log-density
augmented model does not carry).

Then re-point `test-focei-ll-fast-grad.R` / `test-focei-ll-fast-grad-fit.R`.

## Phase 6 -- delete the R implementation

Only after 3-5 land and their tests pass on the C++ entry.  Deleting earlier would
silently drop FOCE/AGQ/ll fits to finite differences.

* `src/inner.cpp`: remove the `Rcpp::Function _agf` fallback block in `analyticOuterGrad`
  (~4649-4700), leaving `foceiOfv0` + `analyticOuterGradDirect`.  Re-check whether
  `restoreFitSolve_()` / `odeSwapRepin()` are still needed anywhere -- the direct path does
  not use them.
* `R/foceiGradAnalytic.R`: delete `.foceiCalcGradAnalytic`, `.foceiAnalyticGradFocei`,
  `.foceiGradAnalyticCalc`, `.foceiAnalyticGradCore`, `.foceiAnalyticGradCoreLL`, the
  config solvers, `.foceiAnalyticSolveAll`'s `rxSolve` fallback if nothing else uses it,
  the `NLMIXR2EST_GRAD_POOLED`-gated experimental block (818-846) together with its now-false
  comment about the protect-stack imbalance, and `.foceiAgqRepData`/`.foceiAgqRepCache`.
  Keep `.foceiGradPooledSetup`, `.foceiGradSolveTol`, `.foceiAnalyticGradSetup`.
* `R/foceiCovAnalytic.R`: `.foceiAnalyticFoceEbeBatch` / `.foceiAnalyticFoceEbe` are also
  used by the covariance path -- check before removing; the cov path is out of scope here.
* `R/vaeGrad.R:95` teardown of `.foceiGradHess2` / `.foceiGradAugNode` / `.analyticStarted`
  follows whatever survives.
* `test-odeswap.R:294-298` compares the pooled route against
  `.foceiAnalyticGradViaRxSolve` -- that comparison loses its second operand.  Decide then
  whether it becomes a pool-vs-`odeSwapNoPool` comparison or is dropped.

## Files touched

* `src/inner.cpp` -- the bulk: `outerSolveFill`, the FOCE/AGQ/ll branches, POD widening,
  the test entry, the direct-route counter
* `src/foceiGrad.cpp`, `src/foceiGrad.h` -- de-`static` two kernels, add `foceiGradSubjectLL_`
* `src/init.c` + `Rcpp::compileAttributes(".")` -- **required**, for the Phase-2b export
* `R/foceiGradAnalytic.R` -- setup widening (Phase 2), deletions (Phase 6)
* `tests/testthat/test-{focei-fast-grad,agq-fast-grad,focei-ll-fast-grad,focei-ll-fast-grad-fit,subFinal,odeswap}.R`
  -- re-pointed at the C++ entry, plus mechanism assertions
* `NEWS.md` -- one bullet under `## Bug fixes` / `### Estimation`
* `CLAUDE.md` -- the odeSwap/gradient-layer note owed by Phase 9 of the parent plan

## Verification

Per phase, in the worktree, with **no concurrent `compile_dll`** -- a build in the same
`src/` while a test's `load_all` compiles there silently destroys the run:

```r
testthat::test_file("tests/testthat/test-focei-fast-grad.R")      # Phases 0-3
testthat::test_file("tests/testthat/test-agq-fast-grad.R")        # Phase 4
testthat::test_file("tests/testthat/test-focei-ll-fast-grad.R")   # Phase 5
testthat::test_file("tests/testthat/test-focei-ll-fast-grad-fit.R")
testthat::test_file("tests/testthat/test-subFinal.R")             # also calls the gradient entry
testthat::test_file("tests/testthat/test-odeswap.R")
testthat::test_file("tests/testthat/test-cov-analytic.R")         # shares the aug-model machinery
```

These assert the analytic gradient against central differences and that `fast=TRUE` fits
match `fast=FALSE` fits -- that is the correctness bar, and it survives the FOCE Newton
being allowed to converge differently.  Each phase additionally asserts
`nAnalyticGradDirect > 0` on a fit of its shape, so a silent fallback fails loudly.

Before the PR: `devtools::test()` on the essential subset, plus an antigravity review of
each commit (`agy -p ... --add-dir /home/matt-fidler/src/nlmixr2est-dry-ode-swap`) with
findings triaged before pushing.

## Risks, worst first

1. **Re-pointing the tests is the highest-risk step, not the C++.** Fifteen assertions
   currently validate an R implementation; if the new C++ entry is not wired to exactly
   the same code path the fit uses, the suite goes green while proving nothing.  Assert on
   the direct-route counter, not just the numbers.
2. **The FOCE Newton EBE re-solve** is the only genuinely iterative algorithm being
   ported, and it is now allowed to differ from R.  It must fail closed on
   non-convergence.  Give it its own commit and check its `eta0Mat` against the R version
   for one fit before wiring the gradient on top.
3. **ES-shape batching.** The shape is a process global; AGQ adds a second slot in the
   same role, ll() a slot in a different role.  A batch opened in the wrong place
   mis-specifies jumps and produces plausible wrong numbers, not a crash.
4. **Pool sizing.** `odeSwapPlan()` sizes the pool for the largest-neq registered model;
   the AGQ order-1 node model and hess2 must be declared before `rxSolve_` builds the pool
   or `odeSwapCanPool` denies them.  Check `odeSwapCanPool(slot) == odeDenyNone`
   explicitly.
5. **`static` kernels and OpenMP exceptions** -- easy to get wrong once, fatal
   (`std::terminate`) when wrong.
6. **Silent scope loss.** Widening `.foceiGradPooledSetup` risks a shape reporting `ok`
   that the C++ branch does not handle.  Every gate removed in Phase 2 must be
   re-expressed per shape, not dropped.

## Sequencing

Phases below are 8F.0 .. 8F.6 of the overall plan:

8F.0 -> 8F.1 (prove no-op on FOCEI) -> 8F.2 -> 8F.3 (FOCE) -> 8F.4 (AGQ) -> 8F.5 (ll)
-> 8F.6 (delete R).  One commit per phase, `origin/main` merged in on each commit,
antigravity review at each commit with findings fixed before pushing, and
`plans/dry-ode-swap-phaseC.md` updated as each phase lands.

Then back to the overall plan: Phase 7, Phase 9, Phase 10, verification sweep, PR.

## Outstanding, independent of this work

`tests/testthat/test-focei-fast-grad.R` has no result against current HEAD -- the last
attempt was destroyed by a concurrent build.  Run it clean before 8F.0 so the baseline is
known.  Nothing has been pushed since `ef8a86a7a` (13 local commits); push once 8F.0/8F.1
are green.

### 8F.3 follow-up: why FOCE declines the direct gradient (measured, not guessed)

A FOCE `fast=TRUE` fit on theo_sd reports `direct 8 / fd 9` -- the direct route declines
about half its evaluations.  The gradient is verified correct whenever it IS produced
(4.1e-06 abs / 9.7e-07 rel against the R route at the initial estimates, against a FOCEI
control at 1.9e-04 / 6.5e-05), so this was never a correctness question.

Instrumented rather than guessed.  Decline reasons, then Newton sub-reasons:

    FOCE   declines -> newton: 9  e0: 0  other: 0
           newton   -> maxit: 9  solve: 0  singular: 0   worst |S| left 0.005498
    FOCEI  declines -> newton: 0  e0: 0  other: 0

So: every decline is the frozen-R0 Newton, and every Newton failure is iteration
exhaustion -- not a failed solve, not an unfactorable Hf.

TWO HYPOTHESES KILLED BY MEASUREMENT:

* "convTol 1e-9 is below the solve's noise floor at sigdig=3 (atol 1e-6)."  Refuted: a
  score stalled at the noise floor would leave |S| ~ 1e-6, not 5.5e-3.  Acting on this
  would have loosened a tolerance that is not the problem.
* "It converges, just slowly; 30 iterations is too few."  Refuted: maxit=200 leaves the
  failure count unchanged at 9 and moves worst |S| only 0.005498 -> 0.004977.  170 extra
  iterations bought 9%.

CONCLUSION: the step is computable (singular: 0) but ineffective -- a stalled undamped
Newton, not a budget problem.  The fix is step control (backtracking on |S|, or a trust
region on ||step||), which costs extra batched solves per Newton iteration and so wants
measuring before it is adopted.

NOT A REGRESSION FROM THE PORT: |S| ~ 0.005 exceeds even the loose first-iteration
skipTol of 1e-3, so R's `.foceiAnalyticFoceEbeBatch` returns NULL at these same thetas.
FOCE fits have always fallen back to finite differences here; the counters are new, the
behaviour is not.

Every one of these declines degrades cleanly to the finite-difference gradient, so this
is coverage and speed, not correctness.

### foceEbeTol: the Newton tolerance from sigdig, and what sweeping it proves

The FOCE frozen-R0 Newton had a hard-coded conv 1e-9 / skip 1e-3 sitting next to a
`sigdig` that drives every other tolerance in the fit.  Now derived:
`conv = 10^-(sigdig+6)`, `skip = 10^-sigdig`, exposed as `foceiControl(foceEbeTol=)`.
At the default sigdig = 3 that reproduces the historic pair exactly, so the change is
behaviour-neutral until a user moves sigdig.

Sweeping it settles the decline question.  theo_sd, FOCE nonmem, 15 outer iterations;
reference objf: FOCE fast=FALSE 118.1402, FOCEI 117.2569.

    foceEbeTol   objf        direct  fd   newtonMaxit
    1e-9         118.3220    8       9    9
    1e-4         117.8059    12      7    7
    1e-2         117.2552    8       0    0

Two things fall out, the second more important than the first:

* Loosening DOES reduce declines -- a prediction that it would not (because the worst
  |S| stalls at 5.5e-3) was wrong.  The stall is not uniform: most subjects converge
  somewhere between 1e-9 and 1e-4 and only the worst sits at 5e-3.

* Loosening is nonetheless NOT a fix, it is a corruption.  Watch the objf walk from the
  FOCE reference to the FOCEI one.  At 1e-2 the tolerance exceeds |S| at the INCOMING
  eta, so the Newton accepts the inner problem's mode immediately and never solves the
  frozen-R0 mode at all.  Zero declines because it stopped doing the work.  That is
  precisely the failure the skip test is only meant to permit for a genuinely stationary
  subject.

So tighter is more correct and the declines are the honest price -- each falls back to
finite differences, which is also correct.  This is why the fix is step control
(backtracking on |S| / a trust region) and not a looser target: a looser target buys a
faster fit that answers the wrong question.

### 8F.4 (AGQ): first attempt corrupted the heap -- ASAN named it, FIXED

Written and compiling, but a fit crashes the allocator on the node solve:

    malloc(): invalid size (unsorted)
    munmap_chunk(): invalid pointer
    free(): invalid next size (fast)

FOCEI and FOCE in the same run completed and matched the R route first (1.9e-04 and
4.1e-06), so the corruption is specific to the new AGQ path.  Stashed as
"8F.4 WIP: AGQ branch" rather than committed.

Found on the way, and correct independently of the crash: **odeSlotOuterNode was
DECLARED but never REGISTERED.**  odeSwapDeclare is metadata-only on purpose (loading a
sensitivity model's DLL before rxSolve_ builds the pool rebinds rxode2's ES globals), so
the node model's neq fed the pool sizing while rxOuterNode kept null entry points.  Any
node solve would have silently reported failure.

LEADING HYPOTHESIS for the corruption, to test first: the node solve runs INSIDE
`OdeSwapEsBatch(odeSlotOuter)`.  odeSwapEsModelForSlot maps odeSlotOuterNode and
odeSlotOuter to the same role (odeEsOuter), on the assumption that the augmented models
share one event-sensitivity shape.  But the node model is order 1 and the gradient model
order 2, so their ES shapes (_rxEsNState/_rxEsNParam/_rxEsNParam2) plausibly differ.
handle_evid sizes its scratch from the effective neq but calls the INSTALLED model's
dydt, and odeSwap.h states a solve may only be compacted when the two agree -- a
mismatch there writes past the buffer, which is exactly the signature seen.

Checks to run before writing more code:
  1. compare odeSwapNeq/odeSwapNlhs and the ES shape of odeSlotOuter vs odeSlotOuterNode
     for an nAGQ=2 fit -- if the ES shapes differ, the node solve needs its OWN
     OdeSwapEsBatch, which means closing the outer batch around it;
  2. call odeSwapCheckLhsWidth(odeSlotOuterNode, &rxOuterNode, rx, op) -- gradPooledCore
     does this for odeSlotOuter and the AGQ path never did;
  3. confirm _aqn and op_focei.aqx/aqw are sized (_aqn x neta) at gradient time before
     aliasing them with the advanced arma::mat constructor.

Do NOT chase this by reading the diff again -- run it under valgrind or ASAN, which will
name the write directly and costs less than another round of hypotheses.

**RESOLVED.**  ASAN named the write on the first run:

    ERROR: attempting free on address which was not malloc()-ed   thread T38
      handle_evid  ->  odeSwapSolveInd  ->  outerSolveFill (the node solve)

The hypothesis above was right, and the mechanism to prevent it already existed.
`OdeSwapEsBatch` keys on the SLOT, not the role -- odeSwap.cpp:174-178 says so in a
comment written for exactly this case: "thetaSens/outer/outerNode/outerCov all share the
odeEsOuter role but are DIFFERENT compiled models with different ES shapes, so a role
match would skip installing the one we are about to solve."  The node loop had no batch
of its own, so the order-1 node model was integrated under the order-2 gradient model's
ES shape and handle_evid freed jump scratch sized for the other model.

Fix: `OdeSwapEsBatch _nodeBatch(odeSlotOuterNode)` around the node loop (outside
outerSolveFill's OpenMP region, since the shape is a process global; the destructor
restores the outer slot for the assembly), plus the odeSwapCheckLhsWidth call for the
node slot that gradPooledCore already makes for odeSlotOuter.

Verified twice over -- ASAN clean, and the gradient matches the R route at the same point:

    FOCEI (control)   max abs 1.891e-04   max rel 6.472e-05
    FOCE  (8F.3)      max abs 4.110e-06   max rel 9.676e-07
    AGQ nAGQ=2        max abs 2.757e-06   max rel 8.190e-07
    AGQ nAGQ=3        max abs 8.355e-07   max rel 2.826e-07

Both were needed: ASAN alone proves only that it stops corrupting memory, and an
uninstrumented run can return plausible numbers while still being wrong.

LESSON, since this is the second time in 8F the fix was "use the existing mechanism
correctly" rather than write new machinery: the ES-shape constraint is documented in
odeSwap.h and was quoted earlier in this very phase, and the bug still got written by
assuming the role mapping implied shape sharing.  Read the slot/role distinction before
adding any further peer solve.
