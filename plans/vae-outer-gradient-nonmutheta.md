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

- **IOV.** `.foceiAnalyticGradSetup` returns NULL on `.uiIovEnv$iovVars`; that is a
  FOCEi-inner concern.  The VAE handles IOV outside the inner problem, so the gate
  does not apply.  NEEDS A SPIKE (see risks) -- the gate's removal is only correct
  if `rxUiGet.foceiOuter` actually emits directions for the IOV etas.
- **Bounded transforms.** `preProcessBoundedTransform.R` rewrites the model before
  estimation, so the model reaching the gradient is already on the unbounded scale.
  The `ui$boundedTransforms` gate is bookkeeping, not a real restriction, for this
  caller.

Still out of scope (unchanged): non-normal `ll()` endpoints, `linCmt()`, `fo`,
multiple estimated boxCox/yeoJohnson lambdas, a theta mu-referenced by more than
one eta.  Note the `foceiControl(fast=)` roxygen ("single additive/proportional
Gaussian endpoint") is STALE -- `.foceiAnalyticErrFull` gates on distribution only,
and multi-endpoint / combined / transformed error all route to the general (f,R)
assembler.  Fix that doc while here.

## Phases

### Phase 0 -- split the scope gate

`.foceiAnalyticGradSetup` hard-gates on `rxGetControl(ui,"fast")`, IOV and
`boundedTransforms`.  Add a caller policy argument (`caller = c("focei","vae")`)
so the VAE can opt in without loosening focei.  Keep every other gate shared.

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

1. **IOV in the augmented model** -- the biggest unknown.  Spike BEFORE committing
   to Phase 0: build the augmented model on an IOV model and confirm the direction
   set covers the IOV etas and the assembled gradient matches an FD reference.  If
   it does not, IOV stays gated and this becomes a follow-up.
2. **Bounded transforms** -- verify what `ui$boundedTransforms` actually holds at
   the point the VAE reaches the setup.  If the hook leaves the field populated as
   a record while the model is already unbounded, the gate is stale and droppable;
   if it signals a live untransformed parameter, it must stay.
3. **Re-baselining.** The gradient converges to the marginal fixed point, not the
   frozen-eta one, so existing `nonMuTheta="regress"` expectations will shift by
   roughly the displacement measured above.  Expect to re-baseline, and do it
   deliberately rather than loosening tolerances.
4. **Cost** -- one augmented solve per M-step vs dozens of full inner sweeps for
   bobyqa should be a clear win; measure rather than assume, since the augmented
   model is larger per solve.
