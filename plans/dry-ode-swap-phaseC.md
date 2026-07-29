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
