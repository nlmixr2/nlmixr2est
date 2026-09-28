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

### Phase 4 -- `fast=TRUE` ends at a wrong point (835.94 vs 772.99)

KM 4.9, omega^2 VMAX .18.  Its cold re-solve agrees and a restart does not
move, so it is an OUTER stop (analytic-gradient optimizer), not a basin trap.
Diagnose separately: check `parHistData`, the stopping rule, and the gradient
vs central differences at the stopped point.
