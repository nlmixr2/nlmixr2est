# Plan: an ELS residual step for `est = "vae"`

Replace the VAE's moment-based residual M-step with an extended-least-squares
(ELS) step that builds the objective from `f` and `r` and optimizes it, following
what SAEM and npag already do.

## Why

Three separate problems, in descending order of how much they matter to users:

1. **Most error models are not estimated at all.**  `.vaeErrTypeCode`
   (`R/vaeFit.R:127`) maps `add -> 0`, `prop -> 1`, **everything else -> 2**, and
   `vaeUpdateErr` leaves code 2 "kept at current value".  A user writing `pow()`,
   `boxCox()`, `yeoJohnson()` or `lnorm()` with `est = "vae"` silently gets their
   `ini()` value back as the estimate.  This is the real defect.
2. **The combined estimator is not maximum likelihood.**  For `add + prop`,
   `vaeUpdateErr` fits `res^2 ~ a^2 + b^2 f^2` by ordinary least squares on
   squared residuals.  The likelihood is `sum[(y-f)^2/sigma^2 + log sigma^2]`.
   The moment fit is a different estimator -- biased and inefficient, and
   sensitive to heavy-tailed residuals.
3. **No endpoint dimension.**  Residual statistics are pooled over all
   observations, so a multi-endpoint model cannot have per-endpoint error
   parameters estimated correctly.

SAEM already solves exactly this: per-endpoint sufficient statistic, a switch on
`res_mod(b)`, closed form where one parameter determines it, and `newuoa`
(default, `saemControl(type=)`) against
`cur = (ytr-ft)/g; sum += cur*cur + 2*log(g)` -- the ELS objective -- otherwise.
npag reached the same design independently (`residOptimize = "alternate"`).

## The parameter set, and what this actually is

**Optimized:** every parameter that is NOT mu-referenced --
(a) structural thetas with no associated eta, and (b) the residual error
parameters.  All of them, uniformly, through the optimizer.

**Left alone:** mu-referenced thetas and the omegas.  Those are what the encoder
and the closed-form M-step exist to estimate; nothing here touches them.

That framing makes this **not a separate ELS step**.  The VAE already optimizes
set (a) -- `nonMuTheta = "regress"` runs a bounded optimizer over exactly those
thetas against the full FOCEi outer objective, re-solving per candidate.  The
work is to **extend that same optimizer to set (b)** and delete `vaeUpdateErr`,
rather than to build a second estimator beside it.

This dissolves the npag concern quoted below rather than working around it.  Its
warning is against mixing *objectives* -- an ELS objective at fixed etas for
residuals against a marginal objective for structural shifts.  Here both sets are
scored by the SAME full outer objective, which already accounts for residual
parameters correctly (it is the FOCEi objective).  One objective, one optimizer,
a larger parameter vector.

It also means `nonMuTheta = "grad"` covers residual parameters for free: the
analytic outer gradient is already defined over "structural theta + residual
sigma + Omega", and `.foceiAnalyticSubjectGradFR` treats a residual sigma as a
pseudo-direction with `df/dsigma = 0` and `dR/dsigma = E$Rsig`, noting that the
`rho(f,R,y)` partials are model-independent closed forms so "ANY variance
structure works".  So the gradient path gets `pow`/Box-Cox/Yeo-Johnson residual
estimation with no new sensitivity columns -- the same models frozen today.

## What is NOT in scope, and why

**This will not close the residual gap against the reference on the neonatal case
study** (`a` 34.9 vs 27.9).  That model is pure additive with `prop.err` fixed at
0, and for pure additive the ELS optimum is exactly `a^2 = SSE/n` -- the moment
estimator already in place.  The two implementations use the same estimator
there; the difference lies in the residuals themselves.  Do not sell this work as
fixing that.

## Design constraints, taken from npag

`R/npCommon.R:261-263` records the constraint the hard way:

> optimize it as a "regressor" against the MARGINAL likelihood over the whole
> support (a structural shift is not identified by the residual ELS step at fixed
> etas) ... Kept SEPARATE from the residual (ELS) set -- mixing the two
> objectives mis-identifies both.

Read carefully, this warns about mixing OBJECTIVES, not about optimizing the two
parameter sets together.  npag's residual step uses an ELS objective at fixed
etas while a structural shift needs the marginal likelihood; scoring both with
the residual objective mis-identifies the structural parameter.  The design above
avoids that by scoring BOTH sets with the full outer objective, so the hazard
does not arise.  What must not happen is reintroducing a separate residual-only
objective alongside it.

npag also needed **warm starts**, because an ELS step evaluated at fixed etas
wanders when the posterior is still poor.  The VAE has that structure by
construction (the M-step runs at the encoder's current posterior means, which are
bad early), so warm starts are required here too, not optional.

## Phases

Each phase is independently committable and leaves the package working.

### Phase 0 -- plumbing, no behavior change

* `vaeInnerLikCore` currently keeps only column 0 of `grabRFmatFromInner`
  (`src/inner.cpp:10384`).  That call already computes **both** `retF` and
  `retR`, so carry `r` alongside `f` at no extra solve cost.
* Add a per-observation endpoint index and a per-error-parameter endpoint index
  to the prep, from `ui$predDf` (`cond`/`var`/`dvid`), which the VAE prep
  currently discards.
* Extend `.vaeErrTypeCode` to a real `res_mod`-style classification per endpoint
  rather than a flat per-parameter 0/1/2.

*Verification*: existing results bit-identical; a unit test asserting the new
`r` matches `(add.err)^2` for an additive model and `(f*prop.err)^2` for a
proportional one.

*Constraint found while doing this*: the stored `r` column is **exactly zero**
wherever `f` is zero (12 of 132 observations on `theo_sd`, the predose records).
The f-throttle that keeps the likelihood finite (`handleF`) is applied when the
likelihood is EVALUATED, not baked into `rx_r_`.  So the ELS objective must apply
the throttle itself rather than consume the column raw, or it divides by zero and
takes `log(0)` at every predose record of a proportional model.

**The rule is `if (r == 0) r = 1`** -- the same convention `handleF` uses for a
zero prediction (`if (adjustF && fa == 0.0) fa = 1.0`).  A unit variance leaves
`(y-f)^2/r + log r` equal to the plain squared residual at that observation,
which is the intended degenerate behavior: the record contributes, but with no
scale information.  Pinned by `test-vae-rvar.R`.

Also measured: the residual parameter round-trips through the inner problem's
scaling with a ~2e-5 ABSOLUTE offset on the SD (0.7 -> 0.7000153, 2.5 ->
2.500021), which roughly doubles in the variance.  Irrelevant numerically, but it
means Phase 1's "reproduces the closed form" check needs a tolerance near 1e-3,
not exact equality.

### Phase 1 -- add the residual parameters to the existing regressor

* Extend the `nonMuTheta = "regress"` parameter vector (`gVaeRegIdx`) to include
  the `err`-tagged thetas with their `ini()` bounds.
* Optimizer choice needs a decision rather than a default: `saemControl(type=)`
  now defaults to **newuoa**, which is unbounded (SAEM keeps its residual
  parameters positive by optimizing their square roots), while npag uses bounded
  **`minqa::bobyqa`** precisely because residual parameters carry `ini()` bounds,
  which the VAE's regressor already honors.  Prefer bobyqa for the bounds and
  offer newuoa as the alternative, rather than inheriting SAEM's default
  unexamined.
* Note npag's rationale for the `log(r)` term: it "keeps the residual from
  drifting to zero on a flexible support (which the marginal likelihood would
  reward), giving the saem/focei residual".  The VAE's encoder is a flexible
  posterior in the same sense, so that term is load-bearing here too -- it is not
  merely a normalizing constant to be dropped for speed.
* **Every parameter in the set is informed by the optimization -- there is no
  closed-form assignment branch.**  A deliberate departure from SAEM, which
  assigns `ares(b) = sqrt(sig2)` for a single-parameter endpoint and optimizes
  only the multi-parameter ones.  One code path, every error model.
* It follows that every residual parameter takes the stochastic-approximation
  step, the way SAEM treats only its optimized cases.  The "form it from the
  sufficient statistic and assign it" rule used for `omega`
  (`omegaUpdate = "suffStat"`) does NOT carry over here -- do not reintroduce it
  by analogy.
* **Keep `vaeUpdateErr`, do not delete it.**  npag keeps its moment estimator for
  two jobs, and both apply here:

  1. **Warm start.**  npag warm-starts each variance scale from the per-endpoint
     moment (additive SD from `sqrt(mean(err^2))`, proportional from
     `sqrt(mean((err/f)^2))`, on the transform-both-sides scale) before running
     the optimizer.  That is what `vaeUpdateErr` already computes, and the VAE
     needs it more than npag does, since the M-step runs at posterior means that
     are poor early in training.
  2. **A closed-form option for the single-additive-error case.**  With exactly
     one additive error the ELS optimum IS `sqrt(SSE/n)`, so the optimizer buys
     nothing: the option returns the identical answer without an optimizer call
     per M-step.  Offer it as a control, defaulting to the optimizer for
     uniformity, and pin the equivalence with the Phase 1 oracle test -- if the
     two ever disagree on a pure-additive model, one of them is wrong.

  It is a shortcut and a warm start, never a second estimator running beside the
  optimizer on the same parameters.

*Verification*: the closed form is an **oracle, not a code path**.  For a
pure-additive endpoint the optimum is `sqrt(SSE/n)`, so the optimizer's answer
must agree to tolerance (~1e-3, per the round-trip note above); a proportional
endpoint gives a second oracle.  If they disagree the objective or the parameter
plumbing is wrong, and that is far cheaper to find here than in Phase 2.
Separately, assert that `nonMuTheta = "grad"` moves a residual parameter, which
exercises the pseudo-direction path.

#### Phase 1 outcome (measured)

Implemented and **opt-in only** (`residOptimize = "optimize"`; the default stays
`"moment"`).  The oracle did its job and found a real limit.

With ONE free residual parameter the optimizer is a clear improvement -- on
`theo_sd` with a combined model and `prop.err` fixed, the objective drops from
143.6 (moment) to 122.6 (optimize), and on a pure-additive model it lands on the
closed form as required (0.79845 against 0.80222, objective 131.757 against
131.812).

With `add` and `prop` BOTH free it diverges: 0.062 / 0.402 at objective 320.7,
against the moment estimator's 0.520 / 0.113 at 122.5.  The two parameters are
near-collinear, and because the M-step gain is 1 until `gammaIter` (60 of 80
iterations in that run) the optimum found at FROZEN etas is adopted undamped and
walks to a corner.  This is not a plumbing bug -- the single-parameter case
proves the objective, the substitution into `a` and the write-back are all
correct.

**Resolved by a two-stage split** (`residOptimize = "twoStage"`), which is what
npag's `residOptimize = "alternate"` does:

* **Stage 1** optimizes the non-mu-referenced structural thetas with the residual
  parameters HELD, so it is driven by `(dv - f)` and the residual does not enter.
* **Stage 2** holds those and optimizes the residual parameters alone against the
  extended least-squares objective `sum[(y-f)^2/r + log r]` over the CACHED
  `(y, f)` pairs.  Since stage 1 fixes `f`, stage 2 needs **no ODE re-solve** --
  the same structure SAEM uses with its cached `_saemYptr`/`_saemFptr`.

Measured on `theo_sd`: the combined model goes from diverging (320.7) to the best
of the three (121.03, against the moment estimator's 122.47), and the additive
model still lands on the closed form (131.79 against 131.81).

The diagnosis was that routing stage 2 through the full OUTER objective lets the
Laplace terms move with the residual at frozen etas; the pure ELS objective at
fixed `f` does not, which is why the near-collinear add/prop pair stops walking
to a corner.  Damping was not needed after all.

### Phase 2 -- combined models (first real behavior change)

* `add + prop` needs no new code once Phase 1 lands -- it is another set of
  bounded parameters in the same optimizer, which is the point of not
  special-casing the single-parameter cases.  What changes is the RESULT: the
  OLS-on-squares moment estimator is replaced by the objective's own optimum.
* Add warm starts (previous iteration's estimate as the start) and a
  burn-in gate, following SAEM's `nb_fixResid` and `warmStartResid`.
* Honor `addProp` (`combined1` vs `combined2`) inside the objective rather than
  in a separate estimator, which removes the special case added in
  `9fc517339`.

*Verification*: `test-vae-errmodel.R` combined case; compare against a FOCEi fit
of the same model, which should now agree more closely than the moment estimator
did.  Expect changed results -- this is the phase that needs release notes.

### Phase 3 -- the frozen error models

* `pow`, `add + pow`, `lnorm` and friends: currently code 2, never estimated.
  With the ELS objective in place these are new cases in the switch, not new
  machinery.

*Verification*: a regression test that the parameter actually moves off its
`ini()` value (the silent-freeze bug would pass any test that only checks
finiteness).

#### Phase 3 outcome (measured, `theo_sd`)

Done for `pow` and `lnorm`, and it did fall out of Phase 1 -- the work was
extending the type classification (`.vaeErrTypeCode`: `pow` -> 3, `pow2` -> 4,
`lnorm` -> 5) and the variance computation inside the stage-2 ELS objective.  No
new machinery.

| model | moment | twoStage |
|---|---|---|
| `pow(prop.err, pw)` | 0.30000 / 0.80000 (= `ini()`), objf 154.390 | 0.41511 / 0.44626, objf 134.848 |
| `lnorm(add.err)` | 0.50000 (= `ini()`), objf 26163.221 | 3.78183, objf 848.856 |

The moment column is the bug: exactly the starting values, untouched.  For
`lnorm` the frozen value is badly wrong for a log-scale residual, hence the 30x
objective difference.

`lnorm` transforms both sides to the log scale inside the objective.  Its
Jacobian does not depend on the residual scale, so it is constant under this
optimization and drops -- which is why `lnorm` lands in Phase 3 rather than
Phase 4.  A Box-Cox/Yeo-Johnson `lambda` that is ESTIMATED does move the
Jacobian, so it must be carried; that is what makes Phase 4 a separate step.

### Phase 4 -- transform-both-sides  [DONE]

Box-Cox and Yeo-Johnson `lambda` are estimated alongside the scale parameters, in
stage 2 (a lambda is an `err`-tagged theta, so it is a residual parameter and
belongs in the residual block).  Nothing was reimplemented: the transform is
rxode2's `_powerD` and the log-Jacobian is `_powerL`, the same pair FOCEi uses
through its `tbs()`/`tbsL()` macros.  Lambda is bounded to `(-2, 2)`; it is
unbounded in `ini()` and only meaningful on a narrow interval, and SAEM likewise
maps it through a bounded transform.

**The objective transforms `dv` ONLY:**

    (tbs(dv) - f)^2 / r + log(r) + jacobian

`f` leaves the solve ALREADY on the transformed scale, because a TBS model
predicts the transformed quantity -- which is why FOCEi applies `tbs()` to `dv0`
and never to the prediction.  Transforming both double-transforms and misfits
badly:

| model | double-transform (wrong) | dv only (correct) | moment |
|---|---|---|---|
| `boxCox` | 158.453 | **43.105** | 181.627 |
| `yeoJohnson` | 383.647 | **120.054** | 131.812 |

Yeo-Johnson under the double transform regressed so far that it looked like a
feature that could not work, and was briefly gated out on a guess about `DV == 0`
handling.  The guess was wrong; the transform was.  Worth remembering that a
plausible-looking regression was a bug in the objective, not a property of the
model -- and that the tests written against the wrong objective passed 25/25 and
would have locked it in.

### Phase 5 -- multiple endpoints  [NOT A TASK]

There is nothing to wire.  The objective already runs across ALL residual
contributors simultaneously: `likInner0` sums a subject's observations, each
carrying its own endpoint's `r` from `rx_r_`, and every endpoint's residual
parameters are in the one optimizer vector.  A model with two endpoints and
different error structures estimates both without per-endpoint statistics.

I had started threading a per-error-parameter endpoint index through the prep,
which was solving a problem that does not exist; it has been reverted.

### The stage-2 objective should use the existing ODE-freeze, not hand-rolled variance

`gVaeElsObjR` currently recomputes `r` itself from the candidate parameters, with
a hardcoded formula per error model (add, prop, combined1/2, pow, lnorm, and the
TBS transform).  **That is a reimplementation of machinery that already exists**
(`src/inner.cpp:10640-10713`):

    npResidFreezeBuild()  normal solves fill a per-subject state cache, then set
                          op_focei.freezeOde = true
    npEvalCondLik()       reuses the cached states and recomputes ONLY r
    npResidFreezeClear()

with the comment already stating the principle: *"params are optimized, f does
not change, so the states can be pinned and only r recomputed -- exactly like
saem's ODE-freeze."*  `npResidELS` is the ELS objective built on it.

Routing stage 2 through `likInner0` with the ODE frozen would be strictly better
than the current hand-rolled form:

* `r` comes from the model's own `rx_r_`, so EVERY error model works with no
  per-form code -- including the ones not yet handled (`type 2`).
* Multiple endpoints and censoring are handled by construction.
* The transform-both-sides bug could not occur: the model applies the transform,
  so there is no opportunity to transform `f` as well as `dv`.  That bug cost two
  rounds here and passed a 25-assertion test suite.

The np helpers are np-coupled (`npEvalCondLik`, `npMixCondLik`, `impNmix`), but
the freeze MECHANISM is generic -- it operates on `op_focei` and the solve cache.
The VAE analogue is to freeze the ODE and evaluate `vaeInnerLikCore` at the fixed
posterior-mean etas over the residual parameters only.  Note this is NOT the
joint solve that diverged: that failed because structural thetas moved at frozen
etas, whereas here `f` is pinned and only `r` varies, which is exactly what npag
does.

### Phase 6 -- control surface and documentation

* `vaeControl(residOptimize=)` mirroring npag's `"alternate"`/`"none"`/`"final"`,
  and a `type=` for the optimizer choice mirroring `saemControl()`.
* Update the `est = "vae"` article: the residual-error section currently
  documents the smoothing-order difference and says matching SAEM would need
  per-endpoint statistics plus an optimizer branch.  That text is the
  specification for this work and should be replaced by what was built.

## Risks

* **Silent behavior change.**  Phase 2 changes estimates for every combined-error
  VAE fit.  It needs to be called out in NEWS, not slipped in.
* **Early-iteration instability.**  The ELS step at bad posterior means is the
  known failure mode from npag.  If warm starts prove insufficient, gate the ELS
  step until after `itersBurnIn` and keep the moment estimator before that.
* **Cost.**  An optimizer call per endpoint per M-step, against cached `f`/`r`
  (no ODE re-solve, unlike the `nonMuTheta` regressor).  Should be minor, but
  measure it -- `nonMuTheta="regress"` is already the expensive part of a VAE fit
  and this must not compound it.
* **Scope creep into the regressor.**  Resist merging the ELS set with the
  `nonMuTheta` regressor set; npag's comment is there because that was tried.

## Suggested order

Phases 0-1 first and on their own: they are a no-op by construction, and they
either validate the plumbing or expose it cheaply.  Phase 3 arguably delivers
more user value than Phase 2 (a silently-frozen parameter is worse than an
inefficient estimator), so if effort is limited, 0 -> 1 -> 3 is a defensible
path that leaves the combined estimator alone.
