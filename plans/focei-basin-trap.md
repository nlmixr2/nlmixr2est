# FOCEi basin trap + small-variance objective cliff (#1132)

Origin: ferx-core epic FeRx-NLME/ferx-core#1564, row A1 (ferx#1565).  ferx's
default FOCEI fit of a seeded Michaelis-Menten ODE benchmark reports
`Converged: YES` at OFV 790.92 when the optimum is 775.39.  Five subjects'
warm-started EBEs stay in a secondary (high VMAX / low V / low KA) mode for the
whole trajectory, and ferx's final cold re-solve check is one-sided: a cold
solve that is BETTER than best-seen passes as "consistent".

## Does it apply here?  Measured

Benchmark: ferx's `sim_bench.py` `mm_n200` (200 subjects, 3 doses q12h, 8 obs,
seed 31), same perturbed inits (VMAX .8, KM 8, V 15, KA .5, omega .2, prop .2),
log-scale mu-referenced thetas, `est="focei"`, `covMethod=""`.  At each final
point: a cold re-solve (`finalUi`, `maxOuterIterations=0`, etas from 0) and a
full restart from the fit's own estimates.  Before = main @ 737ee7efe.

| variant | before | after Phase 1 | cold re-solve diff (after) | restart diff (after) |
|---|---|---|---|---|
| default (`resetEtaP=.15`) | 808.55 | 772.99 | -0.04 | -0.23 |
| `resetEtaP=0` | 1048.18 | 773.01 | +0.03 | -0.23 |
| `mceta=0` | 1030.37 | 772.73 | +0.03 | -0.03 |
| tight ODE (1e-8) | 825.84 | 772.98 | +0.02 | -0.07 |
| `fast=TRUE` | 881.41 | 835.94 | -0.05 | 0.00 |

Before the fix, `resetEtaP=0` and `mceta=0` converged to omega^2 VMAX ~.48
(truth .10), with only generic `$runInfo` notes: the ferx A1 picture.  After the
fix, every variant except `fast=TRUE` lands on the same optimum, and no variant
shows a warm/cold gap.  The "trap" was the cliff below, not flip-flop EBE modes.

### The cliff: residual variance REPLACED by 1 below sqrt(eps)

`src/inner.cpp` `likInner0` (from 85a603e7e, 2022 "NP fixes"):

```c
if (r <= sqrt(std::numeric_limits<double>::epsilon())) { r = 1.0; }
```

With a proportional error (sd .1), any prediction below ~1.2e-3 flips that
observation's `log(R)` from -18 to 0: a +16 step, located by theta, eta and the
ODE tolerance.  Subject 113 (obs t=42.49, DV .0012):
- `rx_r_` 1.4768e-8 gives OBJI -11.88;
- 0.005 away in eta, `rx_r_` 1.4903e-8 gives OBJI -27.95;
- an independent R FOCEi computation (rtol 1e-12, true R) gives about -28.0 at
  both points.

## Plan

### Phase 1 -- remove the cliff (DONE)

- Floor `r` in `(0, sqrt(eps))` at `sqrt(eps)`.  `r <= 0` (a structural-zero
  prediction) keeps the legacy `r=1`; NaN is untouched, as before.
- Where floored, `rp = 0` (a floored R is flat in eta), so the inner gradient
  and the FOCEi Hessian match the objective.
- Left as is: the `sqrt(eps)` gates in `conditionalInnerPartials` and the
  analytic Hessian expansion.  They DECLINE (falling back to FD) rather than
  replace, which is right, because the analytic derivatives do not model the
  floor.
- Test `test-focei-variance-floor.R`: a differential pair straddling the
  threshold (old build: jump 31.97; fixed build < 0.05), plus a closed-form
  FOCEi objective on the unfloored side.

### Phase 2 -- end-of-fit basin check (safety net; re-justify first)

No measured failure remains that Phase 2 would fix.  Before building it, find a
benchmark that still traps after Phase 1 (e.g. ferx's pk2/pk1 twins, or
flip-flop designs).  If one exists:
1. At the final theta, cold-solve each subject from eta=0 plus the `etaRestart`
   Omega draws, and compare on the per-subject MARGINAL objective.
2. Use a two-sided gate with an ABSOLUTE tolerance.  On a gap, run one restart
   leg and keep it only if strictly lower.
3. Add a `$runInfo` warning (< 75 chars) and `$env` counters so tests can
   assert the mechanism ran.

### Phase 3 -- in-trajectory basin re-verification (ferx A2)

Only if Phase 2 is justified.

### Phase 4 -- `fast=TRUE` ended at a wrong point (DONE)

Before: 835.94 (KM 4.9, omega^2 VMAX .18); cold re-solve agreed, restart did
not move.  The analytic gradient was wrong, not stalled: at the TRUE optimum
(772.955) it read lvmax +877 / lkm -755 against near-zero central differences.
It was not inner convergence (`trustFterm=1e-12`, `innerOpt="n1qn1"` left it
unchanged), the endpoint form, or the ODE tolerance.  Per subject it was ONE
subject: subject 5's last obs (DV .001 clip, IPRED 8.5e-4, R 7.6e-9 < floor)
carried the whole error.  `outerSolveFill` read the raw R and its derivatives
while the objective floors R.

Fix: `outerSolveFill` applies the same floor (a shared `foceiRFloor`) and zeros
`aR`/`AR`/`Rsig`/`RsigDir`/`Rsig2` for that row, which is the derivative of the
floored objective.  The same fill feeds the FOCE inner Newton, the LL gradient,
the analytic outer Hessian (still declines at `R <= floor`) and the analytic
covariance through R.  After: the `fast=TRUE` fit reaches 772.87 in 8.3 s
(restart diff -0.04).  The test fails on the Phase-1-only build (analytic -200 /
1888 / -254 vs CD -9 / -0.5 / 86).

Open: `.foceiAnalyticSolveAllFD3` (covariance FD helper) builds its own R
without the floor; only matters for floored rows in `covMethod="analytic"`.

Known Phase 1 residual: `rp=0` where floored drops the 0.5*c^2 term from log|H|,
so log|H| still steps by about 0.02 at the threshold (log(214/210) in the test).
This is the price of a gradient consistent with a flat R.
