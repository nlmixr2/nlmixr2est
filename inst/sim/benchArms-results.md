# Four Bauer arms x four estimators (plan phase 4.1)

Produced by `benchArms.R`, one pass, cells run sequentially so the wall clock is
comparable. `mceta = 0` throughout. Truth is Bauer's `_sim.ctl`; starts are his
`_imp.ctl`, displaced on every arm. The arms differ ONLY in the declared
relative variance, so they are a dispersion ladder:

| arm | lclrv | rv | CV |
|---|---|---|---|
| g1 | -2.408 | 0.09 | 30% |
| g3 | -0.693 | 0.5 | 71% |
| g2 | ~0 | 1.0 | 100% |
| g4 | +0.693 | 2.0 | 141% |

`mareMean` is over CL and V1, `mareRv` over the two relative variances, each
compared as what it means (exp of a log-mean is a clearance). They are separate
because the rv denominator moves 22x across the ladder, so a combined figure is
not comparable across arms -- on g1/focei the combined number was 261% and hid
29.4% on the means against 493.1% on the variances.

| arm | cv | method | mareMean | mareRv | rho | objf | secs |
|---|---|---|---|---|---|---|---|
| g1 | 30% | **saem** | **1.7** | 50.6 | 0.463 | -22220.9 | 218 |
| g1 | 30% | focei | 29.4 | 493.1 | 0.989 | -9271.0 | 258 |
| g1 | 30% | vae | 127.8 | 99.8 | 0.998 | 1240.3 | 347 |
| g1 | 30% | imp | 488.8 | 40.0 | 1.000 | 1898.3 | 1151 |
| g3 | 71% | **saem** | **1.2** | **12.7** | 0.483 | -20936.2 | 325 |
| g3 | 71% | focei | 24.3 | 40.1 | 0.324 | -5485.8 | 156 |
| g3 | 71% | vae | 37.7 | 100.0 | 0.998 | 9357.7 | 241 |
| g3 | 71% | imp | 45.5 | 40.1 | 0.453 | 1115449.8 | 1216 |
| g2 | 100% | **focei** | **2.3** | 30.3 | 0.265 | -3559.4 | 179 |
| g2 | 100% | saem | 7.2 | **15.2** | 0.457 | -17132.2 | 297 |
| g2 | 100% | imp | 28.3 | 29.2 | 0.456 | -3814.2 | 1155 |
| g2 | 100% | vae | 52.7 | 90.2 | -0.999 | 28829.5 | 290 |
| g4 | 141% | vae | 100.0 | 662.6 | -0.999 | 29322.8 | 195 |
| g4 | 141% | saem | 265.5 | 4701.5 | 0.932 | 43257.8 | 322 |
| g4 | 141% | imp | 2358.5 | **64.9** | 0.600 | 1392890665.9 | 978 |
| g4 | 141% | focei | 2590.3 | 66.3 | 0.590 | 4219212.3 | 222 |

truth: CL 5.104, V1 4.711, rho 0.500.

## What it shows

**saem owns the first three arms and loses the fourth.** Means to 1.7 / 1.2 /
7.2 % on g1-g3-g2 with rho within 0.04 of truth every time, then 265% on g4
with rv out by 4701% and rho drifting to 0.932. Performance on the easy arms
does not predict g4.

**g4 (CV 141%) defeats every method.** Best mean is vae's 100% -- which is what
an estimate near zero scores, not a good fit -- and the best rv is imp's 64.9%.
Nothing here recovers this arm, which is why it is the arm this work kept
returning to.

**No method is uniformly best.** saem wins g1 and g3, focei wins g2's means
(2.3% against saem's 7.2%) while losing its variances, and the ordering
reshuffles again on g4. A single "which estimator" recommendation is not
supportable from this table.

**vae's copula is broken independently of dispersion.** `rho` is pinned at
|1| on all four arms -- 0.998, 0.998, -0.999, -0.999 -- and SIGN-FLIPS between
g3 and g2. Its rv sits near 100% on g1-g3-g2 because it drives the variance to
~1e-4 (measured: rvCL 9.96e-05 on g1, 7.16e-05 on g3) and moves that variance
into the residual, which inflates to 3-4x truth. This is the one failure here
that is a clean, arm-independent defect rather than a hard-arm symptom.

**imp costs 4-7x the wall clock of the others** (978-1216s against 156-347s)
without buying accuracy, and its objective is unusable on declared models --
1898, 1.1e6, -3814, 1.4e9 across the ladder while its parameters stay
comparable to focei's. Read imp's MARE, never its objf.

## What it does NOT show

`objf` is comparable down a column, not across a row: focei reports its own
objective and saem's comes from Gaussian quadrature. Cross-method comparison
belongs to the MARE columns.

Two mechanisms were guessed at from partial rows during the run and are NOT
established by this table:

* "focei fails at low dispersion" -- true for g1-g3-g2 (29.4 -> 24.3 -> 2.3)
  and refuted by g4 (2590.3). It is a U-shape, best in the middle.
* "the shared copula estimator collapses for every method but saem" -- drawn
  from the g1 row, where focei/imp/vae all show rho near 1. False in general:
  on g3 and g2, imp's rho is 0.453/0.456 and focei's 0.324/0.265. Only vae
  collapses it everywhere. Swapping in saem's estimator on focei/g1
  (`etaDistCorSuff = TRUE`) stops rho running away (0.989 -> 0.234) and makes
  the fit WORSE (mareMean+rv combined 261% -> 406%, objf -9271 -> -6614), so
  the copula is not the operative cause even where it is visibly wrong.

Single realizations, one seed, no replication -- these are one dataset per arm,
so a cell is an observation, not an estimate of a method's expected error.
