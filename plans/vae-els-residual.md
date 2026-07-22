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

So: the ELS set (`err`-tagged parameters, at fixed etas) and the regressor set
(structural non-mu thetas, against the outer objective) stay separate.  The VAE
already honors this -- `nonMuTheta = "regress"`/`"grad"` is a distinct mechanism.
Do not merge them.

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

### Phase 1 -- ELS objective + optimizer driver, proven as a no-op

* Implement the ELS objective `sum[(y-f)^2/r + log r]` over an endpoint's
  observations, and a `newuoa` driver (mirroring `_saemOpt`, whose default is
  now `newuoa`; keep bobyqa/nelder as fallbacks the way SAEM does).
* Apply it **only** to pure-additive and pure-proportional endpoints, where the
  ELS optimum is known in closed form.

*Verification*: this is the de-risking step.  The optimizer must reproduce
`sqrt(SSE/n)` (additive) and the proportional analogue to optimizer tolerance on
real fits.  If it does not, the objective or the `r` plumbing is wrong, and that
is far cheaper to find here than in Phase 2.

### Phase 2 -- combined models (first real behavior change)

* Route `add + prop` through ELS instead of the OLS-on-squares normal equations.
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

*Verification*: a fit that recovers a known `pow` exponent from simulated data;
a regression test that the parameter actually moves off its `ini()` value (the
current silent-freeze bug would pass any test that only checks finiteness).

### Phase 4 -- transform-both-sides

* Box-Cox / Yeo-Johnson `lambda` estimated jointly with the scale parameters,
  as SAEM's `*Lam` variants do.  Needs `_powerD`-equivalent handling inside the
  objective so `ft`/`ytr` are on the transform scale.

### Phase 5 -- multiple endpoints

* Per-endpoint statistics and per-endpoint optimization throughout, so a model
  with two endpoints and different error structures estimates both.

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
