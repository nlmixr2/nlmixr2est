# RPEM -- Likelihood Model (phi/eta supplied, no optimization)

Foundational component -- read alongside `04-estep.md`. Everything else depends
on this. RPEM needs `p(Y_i | theta_i)` for an externally supplied `theta_i`
(= mu-referenced population params + `eta`), evaluated **directly with no inner
eta optimization**.

## Why neither existing generator is used as-is

- `foceiModel0ll` (`R/focei.R:532`): has the `THETA + ETA` parameterization
  (`.uiGetThetaEta`) we need, but it exists to *optimize* etas in the inner
  Newton problem. Forcing `foceiControl(maxOuterOpt=0, maxInnerOpt=0)` to
  suppress optimization technically evaluates a likelihood, but it is the wrong
  tool -- heavy inner machinery, awkward control path, and (per project owner)
  "terrible". Do not go this route.
- `nlmModel0` (`R/nlm.R:390`): evaluates the population log-likelihood
  **directly at supplied `THETA[]` with no optimization** (exactly the property
  we want) via `.uiGetThetaDropFixed` + `rxGetDistributionFoceiLines` +
  `rxCombineErrorLines` under `nlmixr2global$rxPredLlik=TRUE`, emitting
  `rx_pred_` (obs log-likelihood) and `rx_pred_f_`/`rx_r_` (censoring). But it
  has **no eta dimension** -- all params are `THETA` -- so it cannot represent
  `theta_i = mu + eta`.

RPEM's model is the hybrid: **nlm's direct-evaluation llik model carrying
focei's `THETA + ETA` parameterization**, where per-subject etas enter as inputs
and are never optimized.

## Construction (new `rxUiGet.rpemModel` / `rpemRxModel`)

Template on `nlmModel0` and `foceiModel0ll`:

- **Prefix lines**: `THETA[]` for population fixed effects **and** `ETA[]` for
  random effects -- use the focei-style `.uiGetThetaEta` mapping, NOT nlm's
  THETA-only `.uiGetThetaDropFixed`. Drop/substitute fixed thetas as nlm does.
- **Likelihood lines**: `rxGetDistributionFoceiLines` + `rxCombineErrorLines`
  with `nlmixr2global$rxPredLlik=TRUE` (the same switch both generators flip),
  so the model emits the observation log-likelihood and the
  `rx_pred_f_`/`rx_r_` pieces. This is what makes general error structures,
  multiple endpoints, autocorrelation, and censoring reusable rather than
  re-derived.
- **Output**: per-observation log-likelihood; RPEM sums per subject to get
  `log p(Y_i | theta_i)`. Fix and document the sign convention (nlm negates
  `rx_pred_`; align with the objective reporting in `09`/`10`).
- **Keep it pred/llik-only** (no eta gradient/sensitivity expansion) for the
  base E-step evaluation -- RPEM supplies etas and does not need their
  gradients. Gradients are added only if the numeric M-step (`05`) or SEs (`08`)
  require them.

## How RPEM drives it

- **E-step**: for each MC sample, set `THETA` (population mu / fixed effects) and
  `ETA` (the drawn eta) as inputs, solve once inside the batched `par_solve`,
  and read back the summed log-likelihood -> `p(Y_i | theta_i)`. Etas are
  **supplied data** per (subject, sample), never optimized. In production this
  solve is driven from C++ via rxode2's C API, reading `rx_pred_` from the solve
  buffer (D17); the R-level `rxSolve` used in the spec-13/C1.x tests is a
  test-only convenience.
- **M-step numeric re-scoring**: re-evaluate the same model with updated
  residual/fixed `THETA` against the cached `theta` samples. When only residual
  params change, the structural solution is unchanged -- so cache the structural
  prediction (`rx_pred_f_`) per sample and recompute only `rx_r_`/log-likelihood,
  skipping the ODE. This reuses nlm's existing `rx_pred_f_`/`rx_r_` split.

## Symengine load / prune / cache

- Reuse the load-prune pattern (`loadPruneNlm` `R/nlm.R:442`, `loadPruneSaem`) to
  prune `if/else` branches and load into symengine, producing compiled rxode2
  model text.
- Cache the compiled model by digest, mirroring `foceiModelCache`.

## Open items

- OI-1: RESOLVED -- `ETA` can be supplied in the model's solve path when the
  model emits `rx_pred_`/`rx_r_` (a pred/r output) instead of a DV prediction,
  which is exactly this model's output form. Remaining: confirm `.uiGetThetaEta`
  is the right eta-mapping helper to reuse verbatim vs. a thin RPEM variant.
- OI-2: Confirm residual-only re-scoring can truly skip the ODE by caching
  `rx_pred_f_` per sample; if not, the numeric M-step costs one extra
  llik pass (still cheap vs the E-step). Benchmark.
- OI-3: RESOLVED -- the generated model emits `rx_pred_` as the **negative**
  observation log-likelihood (`llikNorm(...)` then `rx_pred_ <- -rx_pred_`,
  matching nlm), and `rx_pred_f_` as the structural prediction. So the E-step
  uses `log p(Y_i | theta_i) = -sum(rx_pred_)`, and residual-only re-scoring
  reuses cached `rx_pred_f_`. Confirmed by the `rpemModel0` smoke test on a
  1-cmt add-error model.
