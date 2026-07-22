# VAE + log-likelihood endpoints: `nonMuTheta="grad"` and the two-stage regress

## STATUS

**Gap 2 (two-stage regress) is DONE.**  Eligibility is decided per parameter --
an `err` parameter, or one no solve-defining expression reaches -- so an `ll()`
endpoint's log-density-only thetas reach stage 2, including in a multi-endpoint
model that also has `err` rows.  `gVaeFreezeObjR` gained the theta write path.
Verified by unit tests plus two fits (`test-vae-ll-grad-fit.R`).

**Gap 1 (`nonMuTheta="grad"` for `ll()`) is BLOCKED, and the blocker is not what
this plan assumed.**  Everything above the solve works: `.foceiLLGradInScope(ui,
"vae")` is TRUE, `.foceiOuterDirsLL` builds the direction set,
`.foceiAnalyticGradSetup` returns `ef$isLL`, and `.foceiAnalyticGradCore` routes
to the log-density core.  Widening `.vaeGradInScope` to accept `ll()` makes the
fit SEGFAULT.

Measured, on a one-compartment `ll()` model:

- The fault is heap corruption: R aborts inside its own byte-compiler
  (`compiler:::tryCmpfun`) on entry to `.foceiAnalyticGradCore`, i.e. BEFORE any
  gradient solve has run.
- It is not the pooled column layout -- forcing `.vaeOuterCols` to `NULL` (the
  plain `rxSolve` path) still faults.
- It is not the `dH/dtheta` batch -- stubbing `.foceiAnalyticSolveConfigsLL` /
  `.foceiAnalyticHess2ConfigsLL` to `NULL` still faults.
- The SAME model and data fit cleanly under `nonMuTheta="regress"`.

`"regress"` is the only one of those that does not set `poolModel`.  So the
trigger is the `"grad"` pool arrangement itself: `.vaeInnerSetup` sizes the
shared solve pool with the augmented outer model (42 states / 34 lhs here) and
pins the inner MAP to `neqOverride = 4`.  That is correct for a Gaussian model
and corrupts the heap for an `ll()` one.  The fix belongs in `vaeInnerSetup_`,
not in the R scope gate.

`.vaeGradInScope` therefore still declines `ll()`, with the reasoning inline; the
supporting work (caller-aware `.foceiLLGradInScope`, `startedEnv` reuse in
`.vaeGradEval`, `.vaeGradReset` hygiene, and a `predHess2Offset` reset in
`vaeInnerSetup_` so a stale offset from an earlier focei `ll()` `fast=TRUE` fit
cannot leak in) is in place and tested.  Re-enable the `ll()` branch and restore
the two `"grad"` fit tests together once the pool arrangement is fixed.


Two gaps left open when PR #807 (analytic `fast=TRUE` for `ll()`/generalized
endpoints) merged.  Both are about the VAE M-step, which #807 did not touch.

## Established facts

Probed on the merged tree with a Poisson model
(`cp ~ pois(lambda)`, one non-mu structural theta `lb`, `kd` fixed):

| probe | value |
| --- | --- |
| `.foceiLLGradInScope(ui)` | `TRUE` |
| `.foceiOuterDirsLL(ui)` | non-`NULL` (LL direction set builds) |
| `.foceiOuterDirs(ui, "vae")` | `NULL` (Gaussian `ErrFull` declines) |
| `.vaeGradInScope(ui)` | `FALSE` |
| `ui$iniDf$err` | all `NA` -- so `a` (error-param vector) is EMPTY |

## Gap 1 -- `nonMuTheta="grad"` declines every `ll()` model

### Diagnosis

The gradient machinery is already wired end to end for this case:

- `.foceiAnalyticGradSetup(ui, th, Om, caller = "vae")` takes its
  `.foceiLLGradInScope()` branch regardless of caller and returns
  `ef = list(isLL = TRUE)`.
- `.foceiAnalyticGradCore()` dispatches on `ef$isLL` to
  `.foceiAnalyticGradCoreLL()`, so `.vaeGradEval()`'s existing call already
  routes correctly.
- `rxUiGet.foceiOuter` falls back to `.foceiOuterDirsLL()` for an `ll()` model.

The single blocker is the scope probe.  `.vaeGradInScope()` is
`!is.null(.foceiOuterDirs(ui, "vae"))`, and `.foceiOuterDirs()` requires
`.foceiAnalyticErrFull()`, which is Gaussian-only.  `vae.R` (the
`identical(env$vaeControl$nonMuTheta, "grad")` block) therefore downgrades to
`"regress"` and warns "analytic gradient out of scope".

### Change

1. Widen `.vaeGradInScope()` to `Gaussian dirs OR .foceiLLGradInScope(ui)`.
   Keep it a cheap static probe (no symengine/gcc pass).
2. `.vaeGradInit()` sets `.vaeGradEnv$outerCols <- .vaeOuterCols(.am)` to enable
   the pooled `vaeOuterSolve_` read path.  `.vaeOuterCols()` resolves
   `rx_predf_`/`f1`/`f2` against the Gaussian `.foceiAnalyticCols()` layout and
   sets `hasR` from `am$hasRvar`.  **Required guard:** confirm what it returns
   for an LL augmented model.  A `NULL` is safe (the solve stays on the slower
   `rxSolve` path); a non-`NULL` but wrong layout would make the pooled reader
   write the wrong columns.  If the layout is not provably identical, gate
   `outerCols` to the Gaussian case and let LL use the `rxSolve` path.
3. `.foceiAnalyticGradCoreLL()` takes `startedEnv` to reuse `innerHess2` (the
   cheaper eta-only second-order model) for the `dH/dtheta` perturbation batch.
   `.vaeGradEval()` passes none, so LL vae would solve the full outer model
   `2*(nth+nom)` times per M-step.  Give `.vaeGradEnv` the same cached slot
   (`.foceiGradHess2`) and pass it -- this is the difference between usable and
   not on a real model.
4. `.foceiAnalyticGradCoreLL()` returns `etaP = NULL`; `.vaeGradEval()` reads
   only `.r$g`, so no change -- but assert it, since a future Eq-48 warm start
   in the vae path would silently get nothing.

### Verification -- exercised, not asserted

The existing `test-vae-nonmutheta-grad.R` checks scope gates and counters.  For
`ll()` the requirement is a real fit:

- New fit test: a Poisson (or `ll()`) model with at least one non-mu structural
  theta, fit twice -- `nonMuTheta="grad"` and `nonMuTheta="regress"` -- with the
  same seed.  Assert the two land on the same parameters within tolerance, which
  is the same shape as #807's `fast=TRUE` vs `fast=FALSE` assertion.
- Assert the gradient path was actually TAKEN, not silently fallen back:
  `nRegGrad > 0` and `nRegFallback == 0` in the fit's counters (the mechanism
  evidence the repo already requires of fast-path tests), plus absence of the
  "analytic gradient out of scope" `$runInfo` warning.
- Finite-difference cross-check of one `.vaeGradEval()` return against a central
  difference of the M-step objective, at a fixed theta/eta/omega.  This is the
  test that would catch a wrong pooled-column layout from item 2.

Placement: the FD cross-check and scope probes are cheap -- add to the essential
`test-vae-nonmutheta-grad.R`.  The two-fit agreement test is a multi-fit file --
new `test-vae-ll-grad-fit.R` added to a `.slowBatches` batch (batch 6 holds the
VAE internals; check its measured time before adding).

## Gap 2 -- the two-stage regress never engages for `ll()`

### Diagnosis

`residOptimize="twoStage"` (the DEFAULT) is block coordinate descent:

- stage 1: non-mu structural thetas, residuals held -- `gVaeThetaObjR`
- stage 2: residual parameters alone, ODE frozen -- `gVaeFreezeObjR`

Eligibility for stage 2 is `regErrMapV[j] >= 0`, i.e. the regressed theta has a
slot in `a`.  `a` is built in `.vaeDataPrep` from `iniDf[!is.na(iniDf$err), ]`.
An `ll()`/generalized model has NO `err` rows (confirmed above), so:

- `residByOpt` is `false` and `sub` is empty -- **stage 2 never runs**;
- every regressed theta falls into stage 1, so `"twoStage"` silently degenerates
  to a single joint solve, which is exactly the `"optimize"` mode the docs call
  EXPERIMENTAL and warn diverges on near-collinear parameters;
- `residOptimize` is a no-op for these models -- `"moment"`, `"twoStage"` and
  `"optimize"` all do the same thing;
- the `vaeUpdateErr` moment warm start is a no-op.

The stage-2 objective itself is NOT the blocker.  `gVaeFreezeObjR` calls
`vaeInnerLikCore` -- the same likelihood the fit reports, distribution-agnostic.
It is already correct for `ll()`.  Two things are wrong: the eligibility rule,
and the candidate write path.

### Change

**(a) Generalize the eligibility rule from "has an `err` slot" to "the ODE solve
is invariant to it".**  That is the actual precondition for stage 2: it pins the
ODE states from stage 1 and re-evaluates with only the candidate moving, so any
parameter the states do not depend on is eligible.  For a Gaussian model the
`err` thetas are exactly that set, so the rule reduces to today's behavior and
Gaussian fits are unchanged (assert this).  For an `ll()` model it picks up the
thetas that appear only in the log-density expression -- e.g. an overdispersion
or scale parameter -- which is the intended stage-2 population.

Detection: scan the `d/dt()` right-hand sides and state initial conditions for a
transitive dependency on the theta, via rxode2 model vars.  Decide the mechanism
first (a dependency helper may already exist); on any doubt classify as
structural, which is the current, safe behavior.

**(b) Give `gVaeFreezeObjR` the theta write path it lacks.**  It substitutes
candidates only into `a` (`gVaeElsMap[k] >= 0` guards the write) and then calls
`vaeBuildTh(..., gVaeRegErrIdx0, ac)`, which writes `a` over the ERROR theta
positions.  A plain non-`err` theta has neither, so the objective would be FLAT
in the candidate and bobyqa would not move it.  Mirror `gVaeThetaObjR`'s dual
write: always write the theta slot, additionally write the `a` slot when mapped.
This needs a `gVaeEls*` companion to `gVaeElsMap` carrying the theta indices.

**(c) Decide and document the `"moment"` case.**  There is no closed-form moment
estimator for an `ll()` parameter.  Either keep the parameter at its `ini()`
value (today's behavior for `pow`/Box-Cox, already documented) or route
`"moment"` to stage 2 for these models.  Prefer the first -- consistent with the
existing documented gap -- and say so in `?vaeControl`.

**(d) If stage 2 stays empty for a model** (no ODE-invariant regressed theta,
e.g. a plain Poisson whose only free theta is structural), that is legitimate,
not a failure -- but `residOptimize` is then genuinely inert.  Do not warn per
iteration; the roxygen note from (c) covers it.

### Verification -- exercised, not asserted

- New fit test: an `ll()` model with BOTH a structural non-mu theta and an
  ODE-invariant one (a scale/overdispersion parameter inside the `ll()`
  expression).  Fit with `residOptimize="twoStage"` and assert the
  ODE-invariant theta MOVED off its `ini()` value and the objective is no worse
  than the joint `"optimize"` solve on the same seed.  A test that only checks
  the fit completes would pass today, when stage 2 never runs.
- Mechanism evidence: assert stage 2 was entered -- surface a per-fit counter
  (number of stage-2 solves) the way `nRegGrad`/`nRegFallback` already work for
  the gradient path.  Without it a regression in the eligibility scan is silent.
- Regression guard: an existing Gaussian `residOptimize` fit
  (`test-vae-residopt.R`) must be bit-identical after (a).  This is the test
  that protects the reclassification.

## Sequencing

1. Gap 1 items 1-2 (scope widen + `outerCols` guard) + the FD cross-check test.
   Self-contained, R-only, no C++.
2. Gap 1 items 3-4 (`innerHess2` reuse) + the two-fit agreement test.
3. Gap 2 (b) -- the `gVaeFreezeObjR` write path -- with a Gaussian
   bit-identity check.  C++; remember `rm -f src/*.o src/*.so`.
4. Gap 2 (a) -- the eligibility scan -- plus the stage-2 counter and both fit
   tests.
5. (c) docs, `NEWS.md` under `## New features`, and `.slowBatches` placement.

Steps 1-2 and 3-4 are independent and can land as separate PRs.

## Notes

- The `est="focei"` path is unaffected throughout; every change is behind the
  vae M-step or a scope probe that focei does not consult.
- `.vaeGradInScope()` returning `TRUE` for `ll()` means `rxUiGet.foceiOuter`
  builds the LL augmented model for vae fits that previously built nothing --
  confirm the vae shared solve pool is still sized for it
  (`.vaeInnerSetup` `poolModel`), since that pool sizing is what the
  `neqOverride` inner Newton depends on.
