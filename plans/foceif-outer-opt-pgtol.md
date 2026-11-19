# `foceif` stops early: the outer stopping rule, not the analytic gradient

Origin: discussion #924 ("Expectations for `foceif` estimation results compared
to `focei`").  Reported symptom: `foceif` takes a *deterministic* number of outer
steps that tracks the parameter count (8 for a 1-cmt model, 16 for a 2-cmt one)
while `focei`'s count varies with the data, and across six candidate models the
worst model under `foceif` was the best under `focei`.

## What it is not

Three explanations were checked and ruled out **by measurement**:

- **Not the FOCE-FAST phenomenon** (premature termination because an analytic
  gradient walks into a local minimum, per the NONMEM comparison in
  PMID 33884580).  Plain `focei` with finite-difference gradients also stalls at
  8 evaluations when given `outerOpt="lbfgsb3c"`.  No analytic gradient involved.
- **Not the #868 whole-gradient FD fallback.**  Every `foceif` run's
  `parHistData` is 100% `"Analytic Gradient"` -- no `(finite difference)`,
  `(relaxed)` or `(Chartrand)` rows.  The fallback never fired.
- **Not a missing projected-gradient test.**  This was the first fix attempted
  and it was WRONG; see below.

## The `pgtol` dead end (implemented, measured, reverted -- 2b0cae85e / 58ea5715d)

The first hypothesis was that `lbfgsPgtol = 0` suppressed the stationarity test
and left `factr` as "the only stopping rule", so restoring `pgtol` would let the
fit run to stationarity.

**That does not follow.**  `pgtol` is an ADDITIONAL L-BFGS-B stopping rule, not a
replacement for `factr` -- the optimizer stops when EITHER fires.  Enabling it
can only end a fit sooner, so it cannot cure premature termination.

Measured anyway, with the plumbing verified live (`lbfgsPgtolOuter=1e6`
terminates after one evaluation, so the value does reach the optimizer): at
`10^-sigdig` it is not binding.  Bit-identical to `pgtol=0` across five problems
-- same objective to seven digits, same evaluation counts (8, 15, 6, 21, 7).

Reverted.  Keep the lesson: a stopping rule that only ADDS a termination
condition cannot fix late-stopping OR early-stopping in the direction wanted.

## What it actually is

`.foceiFastCtl()` (`R/foceiFast.R`) does not only set `fast=TRUE`; it lets the
outer optimizer re-default, and `R/foceiControl.R` splits on that:

```r
outerOpt <- if (isTRUE(fast)) "lbfgsb3c" else "bobyqa"
```

So `focei` vs `foceif` is **bobyqa vs lbfgsb3c**.  Two effects, of very different
size:

1. **`lbfgsFactr` is too loose** (`10^-sigdig` relative, on the objective
   reduction of a SINGLE step).  Real, and fixable by tightening.
2. **Derivative-free beats gradient-based on hard surfaces**, by far the larger
   effect, and NOT fixable by any tolerance.

Effect 2 is not about the gradient's provenance.  On theo_sd with a block omega:

| config | gradient | outer optimizer | OFV |
|---|---|---|---|
| `focei` bobyqa | finite difference | derivative-free | **95.69** |
| `foceif` nlminb | analytic | gradient-based | 107.03 |
| `foceif` lbfgsb3c | analytic | gradient-based | 108.30 |
| `focei` nlminb | **finite difference** | gradient-based | **113.19** |

`focei`+`nlminb` uses FD and lands worst of the four.  The optimizer class is the
discriminator.

## Phase 2 -- tolerance sweep (DONE)

Sweeping `lbfgsFactr = 10^-(sigdig+k)/eps` and nlminb's `rel.tol`/`x.tol` =
`10^-(sigdig+k)` for k = 0,1,2,3,4,6,8, against the `focei`/bobyqa reference, on
five problems: theo 1cmt, theo 1cmt with a block omega, pheno_sd (sparse),
1-cmt Michaelis-Menten, 2-cmt oral.  `pgtol` pinned to 0 throughout so this
measures the factr/rel.tol axis alone.

Results (`dOfv` vs the bobyqa reference, `evals` = outer evaluations):

| problem | k=0 (was) | k=1 | k=2 (new) | k>=3 |
|---|---|---|---|---|
| theo 1cmt, lbfgsb3c | +0.403 / 8 | -0.004 / 14 | **-0.003 / 21** | -0.003 / 21 |
| pheno, lbfgsb3c | +0.078 / 7 | +0.011 / 8 | **+0.008 / 9** | +0.004 / 16 |
| oral 2cmt, lbfgsb3c | +8.63 / 6 | +2.45 / 8 | **-0.096 / 13** | -0.174 / 16 |
| theo corOm, lbfgsb3c | +12.61 / 15 | +12.61 / 15 | +12.29 / 29 | +12.09 / 52 |
| oral MM, lbfgsb3c | +2173 / 21 | +2173 / 21 | +2173 / 21 | +2173 / 21 |
| theo 1cmt, nlminb | +0.098 / 43 | -0.003 / 51 | -0.005 / 55 | -0.006 / 72 |
| pheno, nlminb | +0.020 / 13 | +0.003 / 17 | +0.003 / 31 | +0.003 / 31 |
| oral 2cmt, nlminb | +0.083 / 22 | +0.104 / 34 | +0.104 / 34 | +0.104 / 34 |
| theo corOm, nlminb | +11.34 / 41 | +11.34 / 41 | +11.34 / 41 | +11.34 / 41 |
| oral MM, nlminb | +1953 / 50 | +1953 / 50 | +1953 / 50 | +1953 / 50 |

bobyqa reference evaluation counts: 65, 47, 96, 278, 57.

Two clean groups:

- **Tolerance-sensitive** (theo 1cmt, pheno, oral 2cmt).  k=2 matches or beats
  bobyqa on all three, at a fraction of its evaluations.  k=1 is NOT enough --
  it still leaves +2.45 on the 2-cmt model, the class the report was about.
- **Tolerance-insensitive** (block omega, Michaelis-Menten).  Flat across six
  orders of magnitude, for BOTH optimizers -- MM does not move by a single digit.
  The optimizers converge; they converge somewhere else.  Not a tolerance bug.

nlminb gains nothing from tightening and is slightly worse on the 2-cmt fit
(+0.083 -> +0.104), so its `rel.tol`/`x.tol` stay as they are.

## Phase 3 -- the default decision (DONE)

`lbfgsFactr` defaults to `10^(-sigdig-2)/.Machine$double.eps` (was
`10^-sigdig/...`).  Nothing else changes.

Blast radius is narrow, and narrower than it looks:

- `op_focei.factr` is read once and used at an inner and an outer call site, but
  the inner one is DEAD: `src/inner.cpp` sets `bool n1qn1Inner = true;` and never
  assigns it, so the inner eta problem always uses n1qn1 and the `lbfgsb3C`
  branch below it cannot run.  `factr` reaches only the outer problem.
  (`FdInnerTolGuard` divides `factr` by `k` for the same reason inertly; it is
  `epsilon`, which n1qn1 does use, that does the work there.)
- With `fast=FALSE` the default outer optimizer is `bobyqa`, which ignores
  `factr` entirely.  So the change lands on `fast=TRUE` (what was measured) and
  on anyone who explicitly asks for `outerOpt="lbfgsb3c"`/`"L-BFGS-B"`.

Deliberately NOT flipping the `fast=TRUE` outer optimizer to nlminb: it wins 5/5
against lbfgsb3c, but both sit well above bobyqa on the two hard surfaces, so
that change would trade one poorly-evidenced default for another without
addressing effect 2.

### Known wrinkle: the constant sits below the solver noise floor

`sigdig` also sets the ODE solver's `rtol = 10^-sigdig`, so a `factr` target of
`10^(-sigdig-2)` asks for a relative objective reduction two orders below the
precision the solve itself claims.  In principle that is chasing noise.

The predicted symptom -- the optimizer spinning on noise -- does NOT show up in
the measurements.  Evaluation counts saturate (13 -> 16 -> 16 -> 16 on the 2-cmt
fit; 21 flat on theo) rather than growing, because L-BFGS-B stops once it cannot
make progress, and the objective improves rather than degrading.  Kept as a note
to revisit if a model ever does show the outer optimizer spinning at large
`sigdig`, not as an open defect.

### Floor at factr = 1

`factr` is a MULTIPLE of machine epsilon, so `factr < 1` requests a reduction
smaller than eps itself -- unreachable, and with `lbfgsPgtol = 0` by default it
would leave L-BFGS-B with no active stopping rule at all.  The two-order
tightening moves that boundary from `sigdig >= 16` to `sigdig >= 14`, so the
default is floored at 1.  (Found by the Antigravity review; the boundary existed
before this change, at a `sigdig` nobody uses.)

### Dead code noticed in passing (not touched)

- `src/inner.cpp:3203` -- `bool n1qn1Inner = true;` is never assigned, so the
  inner `lbfgsb3C` branch and the `innerOpt="BFGS"` option cannot execute.
- `struct FdInnerTolGuard` (`src/inner.cpp:12220`) is defined but never
  instantiated anywhere.

Both are out of scope here; recorded so the next person does not re-derive them.

## Phase 4 -- regression test

Pin the converged objective across outer optimizers so a future default change
cannot silently move it.  Per the repo convention the test must also assert the
mechanism ran (the `parHistData` gradient type is analytic), not just that the
numbers agree.

## Open question

The five benchmark models were written for this exercise.  The Michaelis-Menten
and block-omega fits may be badly posed or badly initialized rather than
revealing an optimizer defect, which would change the reading of effect 2.  Worth
a second opinion from someone who knows these surfaces before drawing conclusions
about `foceif` accuracy in general.
