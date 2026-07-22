# Exact outer gradient for the VAE non-mu thetas

Replace the `nonMuTheta="regress"` bobyqa M-step with the FOCEi analytic outer
gradient (Almquist 2015), solved as one complete augmented `rxSolve` per M-step --
the same way `foceiControl(fast=TRUE)` does it.  No parallelization, no
per-individual loop.

## Why

`.foceiCovAnalyticDirs` (`R/foceiCovAnalytic.R:173-193`) already gives a theta
with no eta its OWN true-sensitivity direction, so `d(OFV)/d(theta)` for exactly
the parameters the regressor hunts is exact and falls out of the same sensitivity
solve.  Measured on theo_sd (`v <- exp(tv)`, tv non-mu; scripts in
`scratchpad/envbias/`):

- **The eta-sensitivity objection does not differentiate the two.**  Perturbing
  eta by iid noise degrades the analytic gradient and the exact gradient of the
  frozen-eta objective bobyqa descends at the same rate (at 1 eta SD: tv -1480 vs
  -1425, tka 383 vs 383).  It is a property of working at frozen non-EBE etas,
  which the incumbent has in full.
- **At realistic encoder accuracy both are fine.**  A healthy VAE fit has encoder
  eta error 0.11 eta SD; there `cos(gA, truth) = 0.947`, `cos(gG, truth) = 0.960`,
  `cos(gA, gG) = 0.961`.  Real encoder error is structured, not iid, and far more
  benign than the synthetic sweep at the same RMS.
- **The gradient has the better fixed point.**  `gA(eta*)` is ~0 at the FOCEi MLE;
  the frozen-eta objective is not stationary there, and its minimizer is displaced
  (tv 3.4299 -> 3.4144 joint / 3.4236 tv-alone; the real VAE fit landed on
  3.4175).  Classic joint-vs-marginal inconsistency.  ~0.4% on tv here -- small,
  but systematic and directional.
- **Magnitudes are inflated ~100x** at encoder etas while direction survives, so
  the step rule must be scale-invariant (Phase 4).

## Scope

In scope for the VAE caller, beyond what focei allows today:

- **Bounded transforms -- IN (spike done).** `preProcessBoundedTransform.R` sets
  `boundedTransforms` on the ALREADY-REWRITTEN ui
  (`.newUi$boundedTransforms <- transforms`, `R/preProcessBoundedTransform.R:280`);
  by then the model is on the unconstrained `rxBoundedTr.*` scale.  The field is a
  record for the post-estimation back-transform hook, not a live constraint.
  focei gates on it because it REPORTS a natural-scale gradient that would need a
  Jacobian correction; the VAE consumes the gradient internally on the same
  unconstrained scale it optimizes, so no correction arises.  Droppable here.

- **IOV -- OUT (spike done, was risk 1).** Not a gate lift.  Two blockers:
  1. `.foceiAnalyticDirections` DOES yield a direction set
     (`ETA_1_|ETA_2_|ETA_3_|THETA_3_` on a 2-eta + 1-IOV model), but `ETA_3_` is
     the SINGLE ui-level `iov.cl`.
  2. The RUNTIME eta vector is expanded per occasion: a refit reports `etaMat`
     columns `eta.ka, eta.cl, rx.iov.cl.1, rx.iov.cl.2` -- 4 etas against the
     direction set's 3.  The augmented model would emit `df/d(iov.cl)` where the
     assembly needs per-occasion `df/d(rx.iov.cl.k)`, which differ (each applies
     only on its own occasion's records).

  Making IOV work needs the direction set expanded per occasion AND the augmented
  model taught the occasion structure -- real work in the builder, not a policy
  change.  Deferred; the IOV gate STAYS for both callers.

Still out of scope (unchanged): non-normal `ll()` endpoints, `linCmt()`, `fo`,
multiple estimated boxCox/yeoJohnson lambdas, a theta mu-referenced by more than
one eta.  Note the `foceiControl(fast=)` roxygen ("single additive/proportional
Gaussian endpoint") is STALE -- `.foceiAnalyticErrFull` gates on distribution only,
and multi-endpoint / combined / transformed error all route to the general (f,R)
assembler.  Fix that doc while here.

## Phases

### Phase 0 -- split the scope gate

There are TWO gate sites, and both must be threaded or the model never builds:
`.foceiAnalyticGradSetup` (`R/foceiGradAnalytic.R:856-895`) and `.foceiOuterDirs`
(`R/foceiGradAnalytic.R:987-998`, which `rxUiGet.foceiOuter` calls).  Add a caller
policy resolved from the ui's controls -- `fast=TRUE` -> "focei",
`nonMuTheta="grad"` -> "vae" -- and skip only the `boundedTransforms` gate for
"vae".  Every other gate stays shared, IOV included.

### Phase 1 -- build the augmented model at VAE setup

Mirror `R/focei.R:1294-1310`: call `rxUiGet.foceiOuter(list(ui))` once, stash
`augMod` at top level and the direction metadata separately so the qs2 `rxLoad`
reloads it.  Key the VAE model cache on the new control flag, exactly as focei
keys `fast` (`R/focei.R:1437-1439`) -- otherwise a non-gradient build gets reused
for a gradient fit.

### Phase 2 -- inner/outer solve-state save/restore

This is the "two models" handling, and it is save/restore, not parallel
discipline.  The augmented `rxSolve` calls `rxSolveFree()` and replaces the global
solve that `vaeInnerLikCore` reads through `getRxSolve_()`.

`storeCovSolveArgs_` is called from `foceiSetup_` (`src/inner.cpp:5345`) -- a path
`vaeInnerSetup_` does NOT take.  So:

1. Call `storeCovSolveArgs_(obj, rxControl, params, data)` in `vaeInnerSetup_`.
2. After the augmented solve, `restoreFitSolve_()` before any further
   `vaeInnerLikCore` call.
3. Treat a failed restore the way `analyticOuterGrad` does
   (`src/inner.cpp:4137-4146`): fall back for that M-step, warn once.

Because the solve is serial and happens between gradient steps, no per-thread lhs
scratchpad is needed.  If a later phase ever moves this inside the OpenMP region,
it would need the `op_focei.thetaSensNlhs`-style thread-local lhs buffer the imp
M-step required -- explicitly out of scope here.

### Phase 2b -- SUPERSEDES 2: pooled solve, no free/restore at all

**Done as Phase 2, and it works, but it is the wrong mechanism.**  Going through
`rxode2::rxSolve` for the augmented model calls `rxSolveFree()` and rebuilds the
global solve, so every M-step pays a full `rxSolve_` setup (~48ms x 250 M-steps
== the entire 12s regression measured against bobyqa on theo_sd).

Swapping the solve state instead is NOT possible from nlmixr2est: subject records
are views into the single process-global `_globals` (`ind->all_times =
&_globals.gall_times[...]` rxData.cpp:4135, `ind->cov_ptr = &_globals.gcov[...]`
:4264, `ind->solve = &_globals.gsolve[...]` :4816), and even the `indOwnAlloc`
path keeps `ind->linH = &_globals.gLin[...]`.  Stashing `rx->subjects`/`rx->op`
restores dangling views once the augmented solve reallocates `_globals`.

The right mechanism is the one imp and advi already use here: **size the shared
solve pool for the LARGER model (the augmented one) and swap FUNCTION POINTERS to
solve the smaller inner model inside it.**  The solve buffer is a scratchpad; what
matters is accumulated into the output vectors.  No second solve context, no free,
no restore.  `vaeInnerSetup_` ALREADY does exactly this for `est="advi"`
(`src/inner.cpp:10152-10170`), so this follows an in-file precedent, not a new
pattern:

1. **Pool sizing.** `.vaeInnerSetup` passes the augmented model as the pool model
   (the `_impPoolModel` slot, `src/inner.cpp:4897`, consumed at :5359) so
   `foceiSetup_`'s `rxSolve_` sizes every per-subject buffer for it.  Record the
   inner model's state count in `op_focei.innerNeq` so inner solves run with
   `ind->neqOverride` (`impSetInnerNeqOverride`, :9458) -- the existing call in
   `vaeInnerSetup_` already fires when `innerNeq > 0`.
2. **Function pointers.** Register the augmented model into a new `rxSolveF
   rxVaeOuter` peer of `rxInner`/`rxPred`/`rxThetaSens` via `rxUpdateFuns`, and a
   `vaeOuterOde(id)` macro alongside `thetaSensOde` (:66).  Record its `neq` and
   its true lhs width.
3. **The solve.** New export `vaeOuterSolve_()`: per subject, `IndNeqOverrideGuard
   neqGuard(ind, vaeOuterNeq)`, `iniSubjectE(..., rxVaeOuter.update_inis)`, solve,
   then `rxVaeOuter.calc_lhs(...)` into a **thread-local** buffer sized to the
   augmented lhs width.  This buffer is mandatory, not an optimization: the imp
   M-step bug was exactly `calc_lhs` overflowing the inner-sized per-thread lhs
   slice (`op_focei.thetaSensNlhs`).
4. **The seam.** `.foceiAnalyticSolveAll` already reduces the whole solve to one
   column matrix up front ("Extract every sensitivity column from the WHOLE solve
   as a matrix ONCE"); have it take that matrix from `vaeOuterSolve_()` instead of
   `rxode2::rxSolve` for this caller.  Everything downstream in
   `.foceiAnalyticGradCore` is untouched.
5. Drop `storeCovSolveArgs_`/`restoreFitSolve_` from the VAE path (Phase 2's
   `_vaeNeedSolveArgs` flag and the `restoreFitSolve_` call in the M-step go away).

Bonus: this also removes the per-M-step R round-trip, and makes the solve
parallelizable later under the same `_innerParallel` discipline imp uses.

### Phase 3 -- the M-step hook

In `vaeTrainCpp_`, replace the `regIdx` bobyqa block (`src/inner.cpp:13222-13250`)
with a single R callback per M-step:

- assemble theta with the existing `vaeBuildTh`
- pass `etaMat = last.mu - baseline`, current `omega`, data, ids
- R side calls `.foceiAnalyticGradSetup` + `.foceiAnalyticGradCore` -- **already
  prototyped**: `scratchpad/envbias/lib.R::gAnalytic` does exactly this with a
  caller-supplied eta matrix and works
- `restoreFitSolve_()`
- keep only the `regIdx` components (discard the omega/sigma block)

Once per M-step, not per gradient step: `nGradStep` steps still pay one augmented
solve.  Keep the `it > klWarmup` gate.

### Phase 4 -- step rule

Do NOT take a raw gradient step (~100x inflation).  Use a dedicated Adam block
reusing `vaeAdam` -- per-coordinate scale-invariant, so inflation is absorbed --
then project onto the `ini()` bounds with `vaeClampVec` and blend with the M-step
gain `gamma`, exactly as the bobyqa result is blended today.

### Phase 5 -- control + fallback

New `nonMuTheta` value (`"grad"`); `"regress"` keeps bobyqa.  Auto-downgrade to
bobyqa with a `$runInfo` note when out of analytic scope or on a solve/restore
failure.  Keep the note under 75 chars and unprefixed (CLAUDE.md).

### Phase 6 -- independent bug

`nonMuTheta="regress"` diverges without `ini()` bounds: `tv <- 3.45` gave
`tv = 7.9e306`, `add.sd = 5.72`; `tv <- c(2, 3.45, 5)` converged cleanly.
`gVaeThetaObjR` takes its bounds straight from `ini()` with no fallback for
infinite ones.  Fix regardless of this work.

### Phase 7 -- tests

Per CLAUDE.md batching: one quick core test stays essential and must assert the
**mechanism** was used (gradient path taken, not just an equivalent number) --
silent fallback to bobyqa must be loud.  Multi-iteration fit comparisons go into a
`.slowBatches` batch.  Do not use `skip_on_ci()` in batched files.

## Risks

Risks 1 and 2 were spikes and are now RESOLVED in the Scope section above (IOV out,
bounded transforms in).  Remaining:

1. **Re-baselining.** The gradient converges to the marginal fixed point, not the
   frozen-eta one, so existing `nonMuTheta="regress"` expectations will shift by
   roughly the displacement measured above.  Expect to re-baseline, and do it
   deliberately rather than loosening tolerances.
2. **Cost** -- one augmented solve per M-step vs dozens of full inner sweeps for
   bobyqa should be a clear win; measure rather than assume, since the augmented
   model is larger per solve.
