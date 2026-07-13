# RPEM -- Evaluation Criteria

Precise, measurable acceptance criteria for a high-quality final product. Each
criterion states: the metric, how it is measured, and the pass threshold. Every
criterion is classified:

- **[BLOCK]** -- must pass to merge the milestone it belongs to.
- **[TRACK]** -- measured and reported every run; a regression beyond tolerance
  raises an alarm but does not block early milestones (becomes [BLOCK] at M5).

`11-validation.md` defines the four statistical bars; this file adds the
quantitative gates and the non-statistical dimensions.

## 1. Likelihood-model correctness (foundational, `13-likelihood-model.md`)

- **C1.1 [BLOCK]** For a fixed model + dataset + supplied `theta_i`, the RPEM
  likelihood model's `log p(Y_i | theta_i)` equals an independent reference
  (hand-computed for a 1-cmt analytic model, and rxode2-simulated residual
  density for the general case) to within `1e-8` relative error per subject.
- **C1.2 [BLOCK]** Supplying `ETA` as input reproduces, for `eta = 0`, the
  population prediction of the same model to `1e-10`; for `eta != 0` it matches
  a manual `theta = mu + eta` substitution.
- **C1.3 [BLOCK]** Residual-only re-scoring (change `beta`, reuse cached
  `rx_pred_f_`) equals a full re-solve to `1e-10` (proves the ODE-skip is exact).

## 2. Statistical accuracy

- **C2.1 [BLOCK] EM trend monotonicity.** Because `lnL` is a per-iteration Monte
  Carlo estimate (fresh draws, D7), do NOT require sample-to-sample monotonicity.
  Require: the least-squares slope of `lnL` over the rolling `slopeWindow` is
  positive until convergence, and within `[-slopeTol, +slopeTol]` after. A
  strictly decreasing smoothed (windowed-mean) `lnL` trend before convergence is
  a failure.
- **C2.2 [BLOCK] Parameter agreement vs reference.** For each shared parameter,
  `|theta_RPEM - theta_ref| <= 2 * sqrt(SE_RPEM^2 + SE_ref^2)` against SAEM and
  FOCEI reference fits (warfarin, theophylline). Reference values from cached
  fixtures.
- **C2.3 [BLOCK] Objective sanity (not cross-method equality).** Do NOT gate on
  matching the FOCEI OFV (it is an approximation) or the SAEM importance-sampling
  likelihood directly. Instead: RPEM `-2 lnL` agrees within `2` points with an
  independent high-accuracy marginal likelihood (adaptive Gaussian quadrature or
  a large-`m` importance-sampling estimate) on a small model.
- **C2.4 [BLOCK] Simulation recovery / coverage.** Over `S >= 100` simulated
  datasets: mean relative bias per parameter within `+/-5%` for fixed effects and
  `+/-10%` for variance components; empirical 95% CI coverage
  (`estimate +/- 1.96*SE`) within `[0.90, 0.98]` (binomial tolerance for `S`).
- **C2.5 [BLOCK] Paper reproduction.** Reproduce paper Table 1 (analytic, M2) and
  Table 2 voriconazole means/SEs within the paper's reported ranges (RPEM
  column, incl. the `+/-` SEs).

## 3. Numerical robustness

- **C3.1 [BLOCK] Determinism.** Identical results (bitwise on estimates, or
  `<1e-10`) for a fixed `seed` and fixed thread count across repeated runs.
- **C3.2 [TRACK] Seed stability.** Over `R >= 20` seeds on one dataset, the
  across-seed SD of each estimate is `<=` that parameter's reported Monte Carlo
  SE (between-seed noise consistent with the algorithm's own uncertainty).
- **C3.3 [BLOCK] Positive-definiteness.** `Sigma^(k)` stays PD every iteration;
  any `nearPD` repair is logged (not silent) and occurs in `<1%` of iterations on
  the benchmark models.
- **C3.4 [BLOCK] No non-finite leakage.** No `NaN`/`Inf` in estimates, `lnL`, or
  predictions; log-sum-exp stabilization asserted in tests with extreme-scale
  inputs.
- **C3.5 [BLOCK] Edge cases.** Passes with: fixed parameters, zero-omega params,
  bounded (logit) transforms, and subjects with a single observation.

## 4. Performance and scalability

- **C4.1 [TRACK -> BLOCK at M5] Speed vs SAEM.** Matched-accuracy protocol: same
  machine, same thread count, median wall-clock over `R >= 5` runs. Target on
  voriconazole: RPEM/SAEM wall-clock ratio `<= 0.4` (paper reports ~0.25-0.33).
- **C4.2 [TRACK] OpenMP efficiency.** Parallel efficiency of the E-step batch
  `>= 70%` at up to physical-core count (measured as `T(1)/(p*T(p))`).
- **C4.3 [TRACK] Memory.** Sample store within the documented budget for the
  benchmark models; chunked mode engages above the configured threshold without
  changing results (`C3.1` still holds).
- **C4.4 [BLOCK] No collateral regression.** Adding RPEM does not change SAEM/
  FOCEI/nlm results and does not increase their compile or runtime beyond noise.

## 5. Feature parity (the stated goal)

- **C5.1 [BLOCK per feature]** A parity matrix, each row with a passing test that
  matches a SAEM or FOCEI reference within `C2.2` tolerance:

  | Feature | Milestone | Reference |
  |---|---|---|
  | Single + multiple endpoints | M1 | FOCEI |
  | Additive / proportional / combined / power error | M1 | FOCEI |
  | Autocorrelated residuals | M1 | FOCEI |
  | Covariates (mu-ref and non-mu-ref) | M1 | SAEM/FOCEI |
  | Bounded parameter transforms | M1 | SAEM |
  | K>1 mixtures | M2 | paper / SAEM |
  | IOV | M3 | SAEM |
  | Censoring M2/M3/M4 | M4 | FOCEI |

  A feature is "done" only when its row's test is green; parity is the AND over
  all rows in the shipped milestone.

## 6. Integration / API quality

- **C6.1 [BLOCK]** The fit is a valid `nlmixr2FitData`: passes the same
  `print`/`summary`/accessor (`nmObjGet`) surface as a SAEM fit; `addCwres()`,
  `addNpde()`, IPRED/PRED, and VPC all produce without error.
- **C6.2 [BLOCK]** `est="rpem"` dispatch and `rpemControl()` follow
  `sharedControl()`; option names reuse SAEM/FOCEI vocabulary where an analog
  exists (verified by grep before naming).
- **C6.3 [BLOCK]** Fit round-trips: refit, `update`, `nlmixrCbind`, and
  simulate-from-fit succeed.

## 7. Code quality / maintainability

- **C7.1 [BLOCK] Reuse boundaries intact** (`02-architecture.md`): no forked
  residual-error math, no forked accessors, no new ODE solver, no MPI.
- **C7.2 [BLOCK] Style** (CLAUDE.md): camelCase R; ASCII-only source/docs/NEWS;
  terse comments/roxygen/NEWS; no AI-verbose prose.
- **C7.3 [BLOCK] Clean checks.** `devtools::document()` clean; `R CMD check` 0
  errors / 0 warnings, notes justified; C++ compiles under CXX17 with project
  flags and no new warnings.
- **C7.4 [BLOCK] Thread safety.** No data races (disjoint writes or explicit
  reductions); respects `ARMA_DONT_USE_OPENMP`/`EIGEN_DONT_PARALLELIZE`; honors
  the `tests/testthat.R` thread policy (single worker CI/CRAN, capped within-solve
  threads on CRAN).

## 8. Test coverage

- **C8.1 [BLOCK]** Unit tests exist for: likelihood model (C1.x), E-step
  estimators (`n_ik`, `N_i`, `tau`), M-step estimators + MH acceptance ratio
  (Eq 32), stopping rule, and parameter classifier.
- **C8.2 [BLOCK]** Integration tests cover the four validation bars; slow ones
  (voriconazole) are skipped on CRAN but run in CI.
- **C8.3 [BLOCK]** Invariants asserted in tests: `C2.1`, `sum_k w^(k)=1`, PD,
  determinism (`C3.1`).
- **C8.4 [TRACK]** RPEM fits registered in the fixture pre-fit cache once stable.

## 9. Documentation

- **C9.1 [BLOCK]** Roxygen for the `est` method and `rpemControl()` (every
  `@param`/`@return`, terse `@details`); at least one runnable example.
- **C9.2 [BLOCK]** Terse past-tense `NEWS.md` bullet per merged milestone.
- **C9.3 [BLOCK]** `design/rpem/` specs kept in sync: any decision change lands in
  `01-decisions.md` first, then dependent specs.

## Milestone gate summary

- M1 merges only when all `[BLOCK]` items in sections 1-3 and 5(M1 rows), 6-9
  pass, plus an initial `C4.1` timing recorded.
- M2-M4 add their feature rows (C5.1) and re-run sections 2-3.
- M5 promotes the `[TRACK]` performance items (C4.1-C4.3) to `[BLOCK]`.
