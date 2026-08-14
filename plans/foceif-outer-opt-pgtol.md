# `foceif` stops early: the outer stopping rule, not the analytic gradient

Origin: discussion #924 ("Expectations for `foceif` estimation results compared
to `focei`").  Reported symptom: `foceif` takes a *deterministic* number of outer
steps that tracks the parameter count (8 for a 1-cmt model, 16 for a 2-cmt one)
while `focei`'s count varies with the data, and across six candidate models the
worst model under `foceif` was the best under `focei`.

## What it is not

Two plausible explanations were checked and ruled out by measurement:

- **Not the FOCE-FAST phenomenon** (premature termination because an analytic
  gradient walks into a local minimum, per the NONMEM comparison in
  PMID 33884580).  Plain `focei` with ordinary finite-difference gradients also
  stalls at 8 evaluations the moment it is given `outerOpt="lbfgsb3c"`.  No
  analytic gradient is involved in that run.
- **Not the #868 whole-gradient FD fallback.**  The `foceif` run's
  `parHistData` is 9/9 `"Analytic Gradient"`, with no `(finite difference)`,
  `(relaxed)` or `(Chartrand)` rows.  The fallback never fired.

## What it is

`.foceiFastCtl()` (`R/foceiFast.R`) does not only set `fast=TRUE`; it lets the
outer optimizer re-default, and `R/foceiControl.R` splits on that:

```r
outerOpt <- if (isTRUE(fast)) "lbfgsb3c" else "bobyqa"
```

So `focei` vs `foceif` is **bobyqa vs lbfgsb3c**.  The two stop on different
criteria, which is exactly the reported contrast:

- **bobyqa** stops when its trust-region radius reaches `rhoend = 10^-sigdig` in
  scaled parameter space -- a *parameter-space* criterion, so its cost depends on
  the surface.
- **lbfgsb3c** ran with `lbfgsPgtol = 0`, which disables the projected-gradient
  test and leaves `factr` as the ONLY stopping rule.  `factr` tests the objective
  reduction of a **single step**, so it fires as soon as one step is small,
  regardless of how far the gradient still is from zero.  The count tracks the
  parameter count because that is roughly how long L-BFGS-B takes to build its
  limited-memory curvature model.

`pgtol` is the stationarity test that catches this.  It is only supportable on an
analytic gradient: a finite-difference gradient never drives `max|proj g|` below
its own noise floor, so the test would either never fire or fire on noise.

## Phase 1 -- enable pgtol exactly where the gradient is known (DONE)

`.sigdigPgtol(sigdig) = 10^-sigdig` in `R/sharedControl.R`, alongside the
existing `.sigdigFactr`.  Applied where the gradient is analytic, left at `0`
everywhere else:

| surface | gradient | pgtol |
|---|---|---|
| `foceiControl()` inner (`innerOpt="BFGS"`) | sensitivity equations, always | `lbfgsPgtol` derived |
| `foceiControl()` outer, `fast=TRUE` | analytic (Almquist) | `lbfgsPgtolOuter` derived |
| `foceiControl()` outer, `fast=FALSE` | finite differences | `0` |
| `lbfgsb3cControl()` | `gr=.nlmixrOptimGradC`, always | derived |
| `optimControl()`, `L-BFGS-B` + `solveType="grad"` | `gr=.nlmixrOptimGradC` | derived |
| `optimControl()`, anything else | optim's own FD | `0` |

`lbfgsPgtol` had to split into inner/outer: one `op_focei` field served both, but
`FdInnerTolGuard` tightens the inner one and the outer one has to follow `fast`.
`fast` can also be downgraded AFTER the control is built (linCmt(), out-of-scope
`ll()`, mixtures), so `.foceiFinalizeFast()` re-suppresses `lbfgsPgtolOuter` when
that happens and the default was taken -- tracked by `lbfgsPgtolOuterDefault`,
mirroring the existing `outerOptDefault` idiom so a built control round-trips.

## Phase 2 -- benchmark before touching the default (IN PROGRESS)

The first measurement was `theo_sd`, a toy that does not stress the optimizer, so
it cannot carry a default change on its own.  Phase 2 measures OFV / outer
evaluations / wall time / gradient type over problems that do:

- theo 1-cmt (baseline, for continuity)
- theo 1-cmt with a **block omega** (correlated etas)
- **2-cmt oral** (the 16-step case from the report)
- 1-cmt **Michaelis-Menten** (ill-conditioned Vm/Km ridge)
- **pheno_sd** (sparse, few observations per subject)

crossed with `focei`+bobyqa, `focei`+nlminb, `foceif`+lbfgsb3c at `pgtol=0` (the
old behavior), `foceif`+lbfgsb3c at the new derived `pgtol`, and `foceif`+nlminb.

## Phase 3 -- the default decision (BLOCKED on Phase 2)

Deliberately deferred.  Enabling `pgtol` may repair lbfgsb3c on its own, and
lbfgsb3c is the optimizer that actually consumes the analytic gradient -- so
flipping the default to nlminb first would mask that.  Decide from Phase 2:

- lbfgsb3c reaches the same optimum with `pgtol` on -> keep it, change nothing.
- lbfgsb3c still stops short WITH a stationarity test -> flip
  `R/foceiControl.R`'s `fast=TRUE` branch to `nlminb`.

## Phase 4 -- regression test

Pin that `foceif` reaches the same optimum as `focei` across outer optimizers, so
a future default change cannot silently move it.  Per the repo convention the
test must also assert the mechanism ran (the `parHistData` gradient type is
analytic), not just that the numbers agree.
