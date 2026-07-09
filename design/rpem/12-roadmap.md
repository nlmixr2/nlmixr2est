# RPEM -- Roadmap and Milestones

Ordering per D11 (breadth-first likelihood core) reconciled with the wider scope
gathered in the interview. Each milestone is independently mergeable and has its
own acceptance gate. Milestones bias toward small, compartmentalized PRs.

## M1 -- Breadth-first likelihood core (first mergeable increment)

Scope (in build order):
- **First task**: the dedicated RPEM likelihood model (`13-likelihood-model.md`)
  -- phi/eta supplied, no optimization. Nothing else can be validated until this
  returns `p(Y_i | theta_i)` correctly, so build and unit-test it against a known
  hand-computed likelihood before the engine.
- E-step + M-step engine, K=1 (single Gaussian).
- All continuous error structures + multiple endpoints + autocorrelated
  residuals, ALL via the rxode2 `llik` path (no bespoke residual math).
- Hybrid M-step: conjugate `mu`/`Sigma`, numeric fixed effects (residual +
  non-mu-ref structural).
- OpenMP + rxode2 `par_solve`; reproducible RNG.
- Full `nlmixr2FitData`: EBEs, PRED/IPRED, CWRES/NPDE (lazy), objective.
- Fisher-score SEs.
- Stopping rule + `rpemControl()`.

Gate: Bar 2 (match SAEM/FOCEI on warfarin + theophylline), Bar 3 (sim
recovery), K=1 smoke of Bar 1, and an initial Bar 4 timing. EM monotonicity and
reproducibility invariants hold.

Excluded from M1: K>1 mixtures, IOV, censoring (M2/M3/M4 residual censoring).

## M2 -- Mixtures (K > 1)

- Both backends (RPEM-specific + SAEM-representation), benchmark-and-choose
  default (`07-mixtures.md`).
- Gate: Bar 1 full K=2 reproduction of the paper's Table 1; label-switching and
  empty-component guards.

## M3 -- Inter-occasion variability (IOV)

- Reuse `R/iov.R` and the mu-referencing IOV handling SAEM uses; extend the
  parameter classifier and the Gaussian sampling to occasion-level etas.
- Gate: match SAEM on an IOV model.

## M4 -- Censoring (M2/M3/M4)

- Route censored-observation likelihood through `censResid`/`censEst` in the
  E-step likelihood call (may partly come free -- see `04-estep.md` OI-2).
- Gate: match SAEM/FOCEI M2/M3/M4 handling on a censored dataset.

## M5 -- Speed hardening (optional, post-parity)

- Chunked/disk-backed sample store for large `n*K*m_Gauss`.
- MH warm-start, burn-in/thin tuning, load-balanced solve layout.
- Full Bar 4 benchmark write-up vs SAEM/QRPEM/FOCEI.

## Cross-cutting, every milestone

- Keep the reuse boundaries in `02-architecture.md` intact (no forked residual
  math, no forked accessors, no new ODE solver).
- Add/extend tests per `11-validation.md`; register fixtures once stable.
- Terse `NEWS.md` bullet per merged milestone (house style).

## Dependency order

```
M1 --> M2 --> (M3, M4 parallelizable) --> M5
```

Paper Table 1 full reproduction is blocked on M2 (K=2). Everything else
validates at M1.
