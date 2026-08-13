# DRY ODE swapping -- summary

Working record: `plans/dry-ode-swap-phaseC.md` (long, chronological, includes the dead
ends).  This file is the short version.

## What it does

A FOCEi fit compiles several rxode2 models -- the inner model, the 2nd-order `innerHess2`,
the theta-sensitivity model, the augmented outer-gradient model and its AGQ node sibling,
and the prediction model.  Before this work each was solved through its own global
`rxSolveF` struct and its own `ind_solve()` macro, with pool sizing decided ad hoc.

Now they are registered in SLOTS (`src/odeSwap.*`), one pool is sized for the largest, and
peers are swapped per individual via `ind->neqOverride`.

## DRY counts

| | before | after |
| --- | --- | --- |
| per-model `ind_solve()` macros | 5 (`innerOde`, `predOde`, `thetaSensOde`, `hess2Ode`, `nlmOde`) | 0 -- all go through `odeSwapSolveInd(slot, rxId)` |
| ad-hoc "is this model bigger than the pool?" blocks | 2 (second silently won) | 1 `odeSwapPlan()` |
| bad-solve retry loops | 4 copies | 1 `odeSwapRetryCore` |
| analytic outer-gradient implementations | 2 (C++ and a parallel R one) | 1, in C++ |
| R lines in `R/foceiGradAnalytic.R` | -- | -791 |

## Tolerance tiers used to accept a change

- **Tier A (bit-identical).** Refactors that must not move a number: the macro conversion,
  the counter relocation, `outerSolveFill` generalization.  Accepted only at max relative
  delta `0.000e+00`.
- **Tier B (`all.equal` 1e-8, stop above 1e-6).** Changes that alter which buffer or stride
  a solve uses but should not change the answer: the nlm pred override (measured 0).
- **Tier C (vs central differences).** New mathematics: every gradient shape is compared
  against central differences of the objective, NOT merely against the previous
  implementation.  FOCEI 4.3e-15, FOCE 2.6e-15, AGQ 8.1e-11, ll() 5.9e-09, ODE-free 5.6e-07.

The tiers exist because a mechanism counter (`nAnalyticGradDirect`, `pooledSolveN`) proves
only that code RAN.  Multi-endpoint `ll()` reached "grad: analytic" with every component
wrong by up to ~370x (issue #838); `pooledSolveN` was satisfied for months by the very
fallback it existed to rule out.  Numbers, not counters, decide.

## Memory-safety findings

No valgrind table: the two real corruptions were found under ASAN, which named both
directly and cost less than the hypotheses it replaced.

| symptom | cause |
| --- | --- |
| `attempting free on address which was not malloc()-ed`, `handle_evid` -> `odeSwapSolveInd` | AGQ node model solved under the gradient model's event-sensitivity shape.  `OdeSwapEsBatch` keys on the SLOT, not the ROLE -- slots sharing `odeEsOuter` are different compiles with different shapes. |
| `free(): invalid next size` after many fits (historic; had forced multi-endpoint onto the `rxSolve` route) | the same bug.  Fixing it let the multi-endpoint exclusion be lifted, so those fits now get the analytic gradient. |

## Known gaps

- Multi-endpoint general-likelihood (`ll()`, `pois()`, `binom()`) gradients are wrong and
  gated off -- issue #838.
- `est="foce"`/`"focep"` get no analytic gradient at the default `sigdig = 3`; the
  frozen-R0 EBE Newton cannot reach its score target -- issue #836, needs step control.
- Phase 10 (lifting the `imp` serial forcing for the pool-sizing case) is blocked on
  finding a model that actually pins the pool.
- `test-subFinal.R` has not been run on this branch.
