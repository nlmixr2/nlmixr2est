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
4. **The seam -- return the per-subject `E` list directly.**  Do NOT rebuild
   `.foceiAnalyticSolveAll`'s intermediate column matrix and re-slice it in R;
   `vaeOuterSolve_()` emits the `E` structures its loop produces and
   `.foceiAnalyticGradCore` consumes unchanged, dropping a whole representation.
   Per subject (`no` = n obs, `nd` = ndir, `nsig` = length(am$sigTh)):

   | field | shape | when |
   |---|---|---|
   | `f` | `no` | always |
   | `a` | `no x nd` | always (1st-order pred sens) |
   | `A` | `no x nd x nd` | always (2nd-order, symmetric) |
   | `R`, `aR` | `no`, `no x nd` | `am$hasRvar` |
   | `AR` | `no x nd x nd` | `am$hasRvar` (symmetric) |
   | `Rsig` | `no x nsig` | `hasRvar && nsig > 0` |
   | `RsigDir` | `no x nd x nsig` | as above |
   | `Rsig2` | `no x nsig x nsig` | as above (symmetric) |
   | `trans` | list of 4 vectors | `am$hasTrans` |

   `E$y` is attached by `.foceiAnalyticGradCore`, not here.
5. Drop `storeCovSolveArgs_`/`restoreFitSolve_` from the VAE path (Phase 2's
   `_vaeNeedSolveArgs` flag and the `restoreFitSolve_` call in the M-step go away).

**The lhs buffer is ours, not rxode2's.**  The augmented model's lhs width does
NOT match the inner model's, so `calc_lhs` must never write into rxode2's
per-thread lhs slice (that is precisely the imp M-step bug: theta-sens `calc_lhs`
overflowing an inner-sized slice on `dV`).  Allocate the lhs vector outside
rxode2, sized to the augmented model's own width, and pass it to `calc_lhs` --
per thread once parallel.  It also makes the column mapping inspectable, so a
wrong offset shows up as a debuggable index rather than heap corruption.

### Phase 2c -- OPEN CORRECTNESS ITEM: the inner objective needs the outer adjustment

Mixing the inner and outer problems means they must be on the SAME objective, and
right now they are not.  `likInner0` returns `fInd->llik`, the joint (data +
prior) negative log density:

```
fInd->llik = -trace(fInd->llik - 0.5*(etam.t() * omegaInv * etam));   // inner.cpp:1873
```

The FOCEi OUTER contribution `fInd->lik[0]` adds three terms the inner one does
not (`src/inner.cpp:2194-2197`):

```
lik = log(slik) + logH0diag + op_focei.logDetOmegaInv5;   // Laplace determinant + 0.5 log|Omega^-1|
lik += fInd->tbsLik;                                      // DV-transform Jacobian
fInd->lik[0] = lik;
```

So the VAE ELBO (built from `vaeInnerLikCore`'s `obj`, i.e. `llik`) and the
analytic outer gradient (which differentiates the Laplace/marginal objective)
are DIFFERENT functionals.  The theta step and the encoder are currently training
against different objectives.

This is consistent with the measurement that motivated the work: `gA(eta*) ~ 0` at
the FOCEi MLE (marginal), while the frozen-eta joint gradient was large there --
the gap between them IS these terms.  So the observed accuracy win is real, but
it comes from stepping theta on a different objective than the ELBO reports.

Decide before Phase 2b lands:
- **(a) Adjust the inner objective** -- add `logH0diag`, `logDetOmegaInv5` and
  `tbsLik` to what the VAE consumes, so ELBO and theta step share one objective.
  Correct, and the printed ELBO then means the marginal likelihood.
- **(b) Keep them split deliberately** -- theta at its marginal MLE, encoder on
  the ELBO -- and document that the reported ELBO is not the objective the non-mu
  theta is optimizing.

(a) is the coherent choice, and it also makes the reported ELBO comparable to a
focei OFV.  DECIDED: (a).

#### Audit of the other FOCEi-inner consumers (imp family, npag/npb, advi)

The three terms are NOT uniformly required -- which of them belongs depends on
what the algorithm is doing with eta:

| term | mode-based (focei, vae grad) | sampling/integrating (imp, advi) | fixed support point (npag/npb) |
|---|---|---|---|
| `logH0diag` (Laplace det) | REQUIRED | must NOT be there (double counts the integration) | must NOT be there (not a mode) |
| `logDetOmegaInv5` | required | required, unless the method adds its own prior normalization | N/A (no prior by design) |
| `tbsLik` (DV Jacobian) | required | **required** | **required** |

So `logH0diag` is correctly absent from imp/advi/npag -- not a gap.  `npag`
deliberately drops the prior too (`src/npag.cpp:538-540`: "npag ignores the prior
-- it sums llikObs -- but likInner0 always forms it"), which is right for a
nonparametric support point.

**The real candidate gap is `tbsLik`.**  It is accumulated separately
(`fInd->tbsLik += tbsJac`, `src/inner.cpp:1677`) and added ONLY inside
`LikInner2` (`:2197`).  Every consumer that builds its objective from
`likInner0`/`fInd->llikObs` instead therefore omits the DV-transform Jacobian:

- `npEvalCondLik` sums `fInd->llikObs[kk]` -- no `tbsLik`.
- the VAE ELBO path (`vaeInnerLikCore` with `adjOuter=false`) -- no `tbsLik`.
- the imp/advi joint-likelihood paths -- to be confirmed per call site.

This only bites a model with a both-sides transform (`tbsLik` is 0 otherwise), so
it is bounded, but for boxCox/yeoJohnson/lnorm it makes those objectives wrong on
the DV scale and non-comparable to a focei OFV -- and actively wrong when lambda
is estimated, since the Jacobian then moves with the parameter.

FOLLOW-UP (not this branch): add `tbsLik` to the imp/advi/npag objective
assembly, with a transformed-model objective comparison against focei as the
test.  Needs its own before/after baselines since it shifts reported OFVs.

**Parallelize the outer solve.**  Per-subject writes are disjoint, so the loop
parallelizes under the same discipline `vaeInnerLikCore` and the imp M-step use:
`cores = min2(cores, getOpCores(op))`, `solveMethodThreadSafe(op)` gate,
`sortIds(rx, 2)` + `_innerParallel` bracket, `setRxThreadId(omp_get_thread_num())`
inside, one lhs buffer per thread.  Deterministic for a fixed core count (disjoint
writes, no cross-subject fold).

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

### Phase 2b status -- BLOCKED on a parameter-NAME mismatch

`vaeOuterSolve_` is written, compiles, registered in `src/init.c`, and its column
resolution (`.vaeOuterCols`) verifies against the augmented model's real lhs
(`predf=0, f1=1,2,3, f2=4..9, rvarf=10`, 29 lhs).  The `.foceiAnalyticSolveAll`
branch is written.  What blocks it is pool sizing:

```
Error: The following parameter(s) are required for solving:
       ETA_2_, ETA_1_, THETA_4_, THETA_3_, THETA_2_, THETA_1_
```

`foceiSetup_`'s `rxSolve_` matches parameters BY NAME, and supplies the inner
model's `THETA[1]`/`ETA[1]`, while the augmented model declares `THETA_1_`/`ETA_1_`
(the generated-model convention in CLAUDE.md).  The positional LAYOUT already
agrees -- verified: same count, same order, only the spelling differs -- so this
needs the augmented model emitted with the inner model's parameter spelling, the
way `.impmapThetaSensModel` does for the imp pool.  It is not a layout change.

Until that lands both the pool wiring (`.env$poolModel`) and the
`.foceiAnalyticSolveAll` branch (`.vaeGradEnv$outerCols`) are deliberately left
OFF, and the M-step keeps the rxSolve + `restoreFitSolve_` path.

**Do not enable one without the other.**  Running `vaeOuterSolve_` against an
inner-sized pool writes 26 states / 29 lhs into buffers sized for 6 / 6 and dies
with `double free or corruption (!prev)` -- observed, not hypothetical.  That is
the same class of failure the plan warns about for the lhs buffer, and it is why
the buffer is ours; the STATE array is still rxode2's and is sized by the pool.

Next step: emit the augmented model with `THETA[n]`/`ETA[n]` parameter names (or
teach `foceiSetup_` to supply both spellings), then flip both switches together
and re-run the theo_sd comparison -- expect the ~15s gap over `regress` to close.

#### Attempted fix 1 -- rename only the `param()` line: FAILED

`.foceiAnalyticAugModelDirs` (`R/foceiCovAnalytic.R:1367`) deliberately rewrites
the declaration `params(THETA[1], .., ETA[1], ..)` to `THETA_1_`/`ETA_1_`.  Making
that rewrite optional (`innerParamNames=`) and having `rxUiGet.foceiOuter` request
the inner spelling for the vae caller does NOT work -- the same error comes back:

```
Error: The following parameter(s) are required for solving:
       THETA_4_, THETA_3_, THETA_2_, ETA_2_, THETA_1_, ETA_1_
```

Reason: the augmented model's BODY genuinely uses `THETA_1_`/`ETA_1_` as
variables -- they are the sensitivity DIRECTION names (`dirs` is
`ETA_1_|ETA_2_|THETA_3_`), so every generated expression and every generated
column name (`rx_f1_ETA_1_`, `rx_rsig1_4_THETA_3_`, ...) is built from them.
Declaring `param(THETA[1])` does not rebind a body variable called `THETA_1_`.
Renaming the body too would rename the output columns and break
`.vaeOuterCols` / `.foceiAnalyticCols`, so this is the wrong lever.

#### Remaining option: alias the parameters at pool setup

Do not touch the augmented model.  Instead, when `_impPoolModel` is set and its
parameter names differ from the inner model's, have `foceiSetup_` add the pool
model's names as ADDITIONAL columns of the `params` it hands `rxSolve_`, aliasing
the same values positionally (the orders are already verified identical).  That
is a contained `foceiSetup_` change and leaves both models untouched.

Until then: `.env$poolModel` and `.vaeGradEnv$outerCols` stay OFF and the M-step
keeps the rxSolve + `restoreFitSolve_` path.  Enabling either alone is fatal
(26 states / 29 lhs into 6 / 6 buffers -> "double free or corruption").

#### Attempted fix 2 -- alias the params in `foceiSetup_`: WORKS

Neither model is touched.  `foceiSetup_` builds `params` as `THETA[1]`/`ETA[1]`
just before the pool `rxSolve_`; when `_impPoolModel` is set, every pool-model
parameter not already present and matching `THETA_n_`/`ETA_n_` gets an extra
column aliasing the corresponding bracket column (same SEXP, so no copy).  The
orders were already verified identical, and the alias is name-driven, so a pool
model whose names already match is a no-op (imp/advi unaffected).

Result on theo_sd: the pooled path runs every M-step (`nRegGrad=250`,
`nRegFallback=0`) and gives `tv = 3.429396` -- BIT-IDENTICAL to the rxSolve path
it replaces, which is the correctness check that matters: same gradient, different
solve plumbing.  `restoreFitSolve_` and the `_vaeNeedSolveArgs` stash are gone
from the VAE path.

Cost: grad 43.4s -> 36.4s (regress 26.6s), so roughly half the rebuild overhead is
recovered.  The rest is the augmented solve itself (26 states vs 6), which is real
work, not setup.  With one regressed theta on 12 subjects bobyqa is at its most
favorable; the crossover should move with more regressed thetas, since bobyqa's
cost scales with their number and one augmented solve does not.  Not yet measured.

##### Blast radius of the `foceiSetup_` alias -- verified no-op for imp/advi

The alias block runs for any `_impPoolModel`, so imp/advi reach it too.  Checked
directly rather than inferred: `.impmapThetaSensModel` emits params

```
THETA[1], THETA[2], THETA[3], THETA[4], ETA[1], ETA[2]
```

i.e. the bracket spelling, which is exactly what `params` already carries.  Every
one hits the "already present" `continue`, nothing matches
`^(THETA|ETA)_[0-9]+_$`, `addNm` stays empty, and no column is added -- a
structural no-op for imp/advi.  Only the VAE's augmented model (whose parameter
names ARE its sensitivity directions) has the underscore spelling that triggers
it.  focei never sets a pool model at all, so it cannot reach the block.

##### REGRESSION FOUND AND FIXED: the pooled branch leaked into focei

`.foceiAnalyticSolveAll` is SHARED with focei's own fast gradient, and
`.vaeGradEnv` is a package-level environment that lives for the whole session.
Gating the pooled branch on `outerCols` alone meant any focei `fast=TRUE` fit
AFTER a vae grad fit in the same session took the pooled path -- against a pool
sized for ITS OWN inner model.  Caught by running the suites in one session:
`test-focei-fast-grad.R` went 47/47 -> pass=24 fail=23.

It only shows up in that ORDER, which is why the earlier standalone runs were
green and why this needed a same-session test rather than more reasoning.

Fix: the branch additionally requires `isTRUE(.vaeGradEnv$active)`, which
`.vaeGradEval` sets around its own gradient call only (`on.exit` clears it), plus
`.vaeGradReset()` on exit from `.vaeTrain` to drop the fit-specific state.
Regression test added (`a vae grad fit does not leak into a later focei fast
fit`): focei fast -> vae grad -> focei fast, asserting the two focei fits agree
and that `active`/`outerCols` are cleared.

Lesson for any future shared-entry-point branch here: gate on an ACTIVE-call
flag, never on the presence of cached state.

#### tbsLik: what was actually wrong, and what was NOT

The earlier audit claimed npag/npb omit the DV-transform Jacobian.  That was
WRONG, and checking before changing it avoided shipping a double-count:
`likInner0` already folds `tbsJac` into `llikObs` per observation, gated on
`op_focei.isNpag || op_focei.isNpb` (`src/inner.cpp:1806` and `:1834`), with a
comment saying exactly why.  `npEvalCondLik` sums `llikObs`, so the Jacobian is
already in it -- adding `fInd->tbsLik` there would count it twice.  No change
made to npag/npb.

The VAE ELBO path IS missing it (the same gate excludes it), so
`vaeInnerLikCore(adjOuter=false)` now subtracts `fInd->tbsLik` from `obj` (which
is -log p).  Zero for a model without a both-sides transform.

##### An lnorm VAE fit is still far from focei -- but NOT because of this

Checked, because a 500x OFV gap looks like a Jacobian bug: focei OFV 685.6 vs vae
359315.7 on the same lnorm theo_sd model.  It is neither the Jacobian nor the
objective change -- the fit DIVERGED, `tv` reaching 9.1e68.  For this data
`sum(-log(DV))` over the 123 observations is only -180.8, three orders of
magnitude too small to explain the gap.

This is the SAME unbounded-`ini()` divergence recorded as Phase 6 (`tv <- 3.45`
with no `ini()` bounds runs away; `tv <- c(2, 3.45, 5)` converges).  It is a
pre-existing bug, independent of everything on this branch, and it is why an
lnorm/focei OFV comparison cannot currently be used to validate the Jacobian.
Fix Phase 6 first, then re-run this comparison as the real check.

##### Phase 6 fixed -- and what it did NOT fix

The fallback bounds stop the runaway on lnorm too: that fit went from OFV 359315
(`tv` 9.1e68) to OFV 13412 (`tv` 2.446).  So the divergence is gone.

But the lnorm VAE fit is still POOR -- focei gets OFV 685.6 at `tv` 3.857, the VAE
13412 at `tv` 2.446, a 19.6x OFV gap.  That is NOT the transform Jacobian:
`tbsLik` is -180.8 for this data, two orders of magnitude too small to account for
a 12727 gap.  It is fit quality on a log-scale error model, a separate issue from
anything on this branch.

Consequence: **the tbsLik change is still not validated against a focei OFV.**  It
is correct by construction (the term is present for focei/npag/npb and was absent
only from the VAE ELBO path) and the suites are green, but the end-to-end
"comparable to a focei OFV" claim needs a transformed model the VAE fits WELL.
Either find/tune one, or validate the term directly: assert
`ELBO_with - ELBO_without == sum(log|dy'/dy|)` at a fixed parameter vector, which
isolates the Jacobian from convergence entirely.  The latter is the better test
and does not depend on lnorm fit quality.

##### tbsLik in the VAE ELBO: RESOLVED -- the term was right, my reference was wrong

Applying the transform-both-sides Jacobian to the VAE ELBO was first reverted on
the strength of a bad comparison, then restored after instrumenting.  Recorded
because the trap is easy to fall into twice.

`vaeInnerLik` at eta=0 on theo_sd, per subject, `fInd->tbsLik` vs
`-sum(log(DV))`:

| subject has DV==0 | agreement |
|---|---|
| no (subjects 1, 7, 10) | EXACT, ratio 1.0000 |
| yes (the other 9) | off by exactly +18.0218 each |

`18.0218 = -log(sqrt(DBL_EPSILON))`: the transform FLOORS a non-positive DV at
`sqrt(eps) = 1.49e-8`, so a `DV==0` row contributes `-log(1.49e-8)`.  theo_sd has
9 such rows, and

    -180.848  +  9 * 18.0218  =  -18.6516      (measured total: -18.652)

So `fInd->tbsLik` is exactly correct.  The apparent "10x mismatch" came entirely
from a reference that DROPPED the zero-DV rows (`dv[dv > 0]`) while the C code
keeps them at the floor.  The change is restored, and
`test-vae-nonmutheta-grad.R` now pins the identity with the floor applied.

(Whether a `DV==0` under `lnorm()` should be an epsilon-floored density at all is
a separate modeling question -- it is really a below-LLOQ/censoring case.  What
matters here is only that the VAE now treats the term exactly as focei's
`LikInner2` does.)

Still true and unaffected: npag/npb already fold `tbsJac` into `llikObs`
(`inner.cpp:1806/1834`), so nothing must be added there -- that one WOULD
double-count.

##### The non-TBS objective is unchanged, exactly

`add` gives 212.0768500568 with and without the term -- byte-identical -- so the
change is a strict no-op absent a both-sides transform.

LESSON: when validating a C-side accumulator against an R reference, reproduce the
C code's DOMAIN HANDLING (flooring, censoring, dropped rows), not just its formula.
Three separate wrong conclusions on this branch came from a plausible-looking
reference that differed from the implementation in exactly one such detail.

### est="vae" + IOV: localized to .vaeToFit, NOT to nonMuTheta

An earlier note on this branch said "est=vae + IOV is broken on main for every
nonMuTheta mode".  That was too broad and the mechanism was wrong.  Measured:

| configuration | result |
|---|---|
| `test-vae-iov.R` (all-mu, `covariateSelection=FALSE`, `returnVae=TRUE`) | **9/9 PASS** |
| all-mu, `covSel=TRUE`, `returnVae=FALSE` | ERROR |
| all-mu, `covSel=FALSE`, `returnVae=FALSE` | ERROR |
| non-mu `tv`, `covSel=TRUE/FALSE`, `returnVae=FALSE` | ERROR |
| all-mu, `covSel=FALSE`, **`returnVae=TRUE`** | **OK** |

So it is NOT the `nonMuTheta` mode (`eta`/`fix` fail too), NOT the mu-structure
(all-mu fails), and NOT covariate selection (fails with it off).  The single
discriminator is `returnVae`: with `TRUE` the raw VAE object is returned and the
run succeeds; with `FALSE` the fit-object assembly `.vaeToFit` runs and throws
`invalid second argument of length 0`.

**TRAINING WITH IOV WORKS.**  The bug is in building the nlmixr2 fit object from
an IOV VAE result.  That is also why `test-vae-iov.R` is green -- it asserts on
the raw object (`fit$zPop`, `fit$omega`, `fit$prep`) and never exercises
`.vaeToFit`, so the whole fit-assembly path is untested for IOV.

The error carries no call (`conditionCall` is NULL), so bisect `.vaeToFit`
directly: the per-occasion expanded etas (`rx.iov.cl.1/.2`) are the obvious
suspect -- the eta-collapse/output code paths are written against ui-level etas,
and the runtime vector is longer.  Note `.vaeNonMuThetas`/`isFree` handling is a
red herring here; `isFree` is already respected by the covariate B&B
(`inner.cpp:13480`).

Follow-up: add a `returnVae=FALSE` case to `test-vae-iov.R` so the assembly path
is covered.

##### Why the IOV error has no usable traceback

`nlmixr2Est0` catches the estimation error and re-raises it with
`stop(paste(ret$error, collapse="\n"), call.=FALSE)` (`R/nlmixr2Est.R:296`), so by
the time it reaches the caller the original call and stack are gone --
`conditionCall` is NULL and `sys.calls()` under `withCallingHandlers` shows only
the re-raise frame.  Do not spend time on traceback tricks.

To find the line, either instrument `.vaeToFit` directly (it is reachable
standalone: fit with `returnVae=TRUE` to get the raw object, then call the
assembly path on it), or temporarily bypass the catch in `nlmixr2Est0`.

##### IOV: the actual failure site

Stashing `.vaeToFit`'s arguments with `trace()` (no source edit) and replaying it
at TOP LEVEL -- where `nlmixr2Est0`'s catch cannot swallow the condition -- gives:

```
MSG: invalid second argument of length 0
1: some etas defaulted to non-mu referenced, possible parsing error: rx.iov.cl.2
2: In foceiFitCpp_(.ret)
```

So `.vaeToFit` runs a FOCEi pass (`foceiFitCpp_`) to assemble the fit object, and
in THAT pass the per-occasion expanded IOV etas are not paired with a theta --
`rx.iov.cl.2` "defaulted to non-mu referenced" -- and something downstream then
receives a zero-length argument.

Note it names occasion **2**, not `.1`: the first occasion's eta apparently pairs
and later ones do not, which points at the mu-reference map being built from the
ui-level `iov.cl` (one entry) while the runtime has one eta per occasion.  This
matches the mechanism described from memory: the IOV etas need to read as a
`theta.iov + eta.iov` mu term (small held omega, `theta.iov` absorbing
`mean(eta.iov)`), and that pairing is what is missing on the assembly path.

Training is unaffected because `vaeTrainCpp_` treats the IOV etas as `isFree`
(theta forced to 0) and never needs the mu map.

Repro recipe (fast, no source edit):
    trace(nlmixr2est:::.vaeToFit,
          tracer=quote({assign("DBG_env", env, envir=globalenv())
                        assign("DBG_fit", fit, envir=globalenv())}), print=FALSE)
    <run the fit>          # error is caught and discarded as usual
    nlmixr2est:::.vaeToFit(DBG_env, DBG_fit)   # replays it uncaught

##### IOV assembly: two more hypotheses DISPROVEN

1. **The "traceback" I first reported was the deferred WARNING list**, not a call
   stack (`In foceiFitCpp_(.ret) :` is a warning prefix).  The
   "some etas defaulted to non-mu referenced ... rx.iov.cl.2" line is therefore
   EXPECTED noise -- `preProcessBoundedTransform.R` already has
   `.isSyntheticIovMuWarning` / `.filterSyntheticIovMuWarnings` /
   `.getSyntheticIovEtaNames` precisely to muffle it.  It is not the fault.

2. **It is not an eta count/dimension mismatch.**  After `.vaeUpdateModel`, with
   the stashed objects:

       prep$etaNames / ui2 eta rows : eta.ka, eta.cl, eta.v, rx.iov.cl.1, rx.iov.cl.2
       ui2 nEta 5   vs etaMat cols 5   vs ncol(mu) 5   vs omega 5x5

   Everything lines up -- the updated model DOES declare both per-occasion etas.

So `.vaeToFit` hands `foceiFitCpp_` a self-consistent 5-eta problem and the
length-0 error happens INSIDE that call.  Next step is to bisect within
`foceiFitCpp_`/`foceiSetup_` for the IOV case (compare against a working
`est="focei"` IOV fit on the same model, which succeeds -- so the difference is in
what `.vaeToFit` puts on `.ret` versus what the focei path builds).

Remaining budget note: I disproved covariateSelection, the nonMuTheta mode, the
mu-structure, and the dimension mismatch.  Do NOT re-test those.

##### Fifth elimination: `.uiIovEnv` is NOT cleared

`.vaeToFit` snapshots/restores `.muRefTrans$cur` because the focei covariance
recompute re-runs preprocessing and clears staged hook state, so `.uiIovEnv`
(staged by the IOV preprocess hook, consumed by the `.uiFinalizeIov` post-hook)
looked like the same bug.  It is not -- traced:

    AT .vaeToFit ENTRY     iovVars=[iov.cl] ui=TRUE
    AFTER .vaeUpdateModel  iovVars=[iov.cl] ui=TRUE

State is live throughout.

### STOP GUESSING -- change approach

Five hypotheses eliminated by black-box experiment: covariate selection, the
nonMuTheta mode, the mu-structure, eta dimensions, `.uiIovEnv` lifetime.  Each
cost a full fit cycle.  `.vaeToFit` demonstrably hands `foceiFitCpp_` a
self-consistent 5-eta problem with live IOV state, and it still throws
`invalid second argument of length 0` inside that call.

Further black-box guessing is not converging.  The next person should bisect
INSIDE the call instead:

1. `est="focei"` fits this exact IOV model fine.  Dump both envs
   (`.ret` from `.vaeToFit` vs the focei fit env) and diff field by field --
   `names()`, then class/length/dim of every shared name.  The difference is a
   field `.vaeToFit` sets differently or omits, not anything ruled out above.
2. Failing that, the message is R-level ("invalid second argument of length 0" is
   not an Rcpp stop), so it comes from an R callback invoked by the C++ during
   setup/finalize.  Grep the R functions `foceiFitCpp_` calls back into and look
   for a two-argument primitive (`seq_len`/`rep`/`match`/`%in%`/matrix indexing)
   fed a zero-length IOV-derived vector.

##### The env diff -- candidate list is now SHORT

Stashing `foceiFitCpp_`'s argument from both paths (its argument IS the `.ret`
env) and diffing:

    ONLY-IN-FOCEI:  (nothing)
    ONLY-IN-VAE  :  method, vaeControl, parHistData, vae, adjObf, extra
    DIFF etaMat     focei=NULL:0        vae=matrix:12x5
    DIFF ui         focei=rxUi:38       vae=rxUi:39

The VAE-only fields are cosmetic (labels/history that `.vaeToFit` presets to stop
the C++ finalize writing focei labels).  Two substantive differences remain:

1. **`etaMat`** -- focei passes NULL, the VAE passes a 12x5 matrix whose column
   names are the RUNTIME etas `eta.ka, eta.cl, eta.v, rx.iov.cl.1, rx.iov.cl.2`.
   The ui declares `iov.cl` as ONE eta, so any focei-side code that matches
   etaMat columns against ui-level eta names gets no hit for the expanded pair --
   the obvious source of a zero-length second argument.  PRIME SUSPECT.
2. **`ui`** -- one extra binding in the VAE's ui (39 vs 38); worth diffing
   `ls()` on both if (1) does not pan out.

Test in flight when this was written: null out `etaMat` in a `trace()` tracer on
`foceiFitCpp_` (the arg is an environment, so the tracer can mutate it) and see
whether the fit completes.  Baseline reproduced the error in the same script, so
the comparison is valid.  If nulling etaMat fixes it, the real fix is to map the
expanded IOV columns to what the focei setup expects (or drop them from the
starting-eta matrix), NOT to stop passing starting etas.

##### etaMat DISPROVEN too (sixth elimination) -- read this, not the note above

The `etaMat` "prime suspect" call in the preceding section is WRONG.  Nulling it
out at the `foceiFitCpp_` call changes nothing:

    BASELINE          : ERROR: invalid second argument of length 0
    WITH etaMat=NULL  : ERROR: invalid second argument of length 0

(Same script, so the baseline is a valid control.)

That leaves exactly ONE substantive candidate from the env diff: **`ui`, which has
39 bindings on the VAE path vs 38 on the focei path.**  Next step is mechanical --
`setdiff(ls(vaeUi), ls(foceiUi))` on the two stashed envs to name the extra
binding, then work out why the focei setup chokes on it for an IOV model.

Running tally of eliminated hypotheses (do NOT re-test): covariate selection,
nonMuTheta mode, mu-structure, eta dimensions, `.uiIovEnv` lifetime, `etaMat`.

### THE ACTUAL FIX DIRECTION: assemble like SAEM, not like focei

Seventh elimination first: the remaining env-diff candidate (`ui` 39 vs 38) is
NOT structural.  Both uis are identical where it matters --

    focei iniDf rows 10   vae iniDf rows 10
    focei etas  eta.ka,eta.cl,eta.v,rx.iov.cl.1,rx.iov.cl.2
    vae   etas  eta.ka,eta.cl,eta.v,rx.iov.cl.1,rx.iov.cl.2
    focei thetas tka,tcl,tv,add.sd,iov.cl      (same for vae)
    no extra meta bindings either way

-- and `est="focei"` fits this very ui fine.  So nothing about the ui or the etas
is wrong.

**The problem is the assembly ROUTE.**  `.vaeToFit` calls `foceiFitCpp_`, which is
focei's OWN fit path.  SAEM -- a foreign method with exactly the VAE's problem
shape (it computes its own results, then has to produce an nlmixr2 fit, and it
supports IOV) -- never calls `foceiFitCpp_`.  It ends `.saemFamilyFit` with

    .nlmixr2FitUpdateParams(.ret)
    .saemMixFix(.ret, .ui)          # must run before the create call
    .saemAddParHist(.ret)
    .saemCalcLikelihood(.ret)
    .ret$theta <- ...; .ret$model <- ...; .ret$est <- "saem"
    .saemControlToFoceiControl(.ret)
    .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData,
                                      control=.ret$control, table=.ret$table,
                                      env=.ret, est="saem")

`nlmixr2CreateOutputFromUi` is the supported output-construction entry point for a
non-focei method, and it is what makes SAEM+IOV work.

ACTION: rework `.vaeToFit` to follow that sequence instead of driving
`foceiFitCpp_`.  This is an architectural change to the VAE output path, not a
one-line patch, so it wants its own branch and a before/after on the non-IOV VAE
fits (which currently work through the focei route and must not regress).

Full eliminated list (do NOT re-test): covariate selection, nonMuTheta mode,
mu-structure, eta dimensions, `.uiIovEnv` lifetime, `etaMat`, ui contents.

##### Caveat on the SAEM route: code evidence strong, empirical check INCOMPLETE

The route argument rests on code reading -- `.saemFamilyFit` uses
`nlmixr2CreateOutputFromUi` and never `foceiFitCpp_`, and `nlmixr2Est.saem`
carries `attr(..., "iov") <- TRUE`.  That is solid.

What was NOT confirmed before this session ended: an actual `est="saem"` fit of
the SAME IOV model/data succeeding.  The run was still going when work stopped.
Do that first -- it is one fit, and it decides whether the route hypothesis holds:

  * SAEM fits it        -> route confirmed, rework `.vaeToFit` onto the SAEM
                           sequence.
  * SAEM fails too      -> the bug is upstream of the assembly route entirely and
                           this eighth hypothesis dies with the other seven.

Treat the section above as a lead, not a verdict, until that runs.
