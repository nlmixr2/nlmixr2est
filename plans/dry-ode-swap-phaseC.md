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
