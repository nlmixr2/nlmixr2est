# RPEM -- Roadmap and Milestones

Ordering per D11 (breadth-first likelihood core) reconciled with the wider scope
gathered in the interview. Each milestone is independently mergeable and has its
own acceptance gate. Milestones bias toward small, compartmentalized PRs.

## M1 -- Breadth-first likelihood core (first mergeable increment)

Progress (as of 2026-07-09): the K=1 core is built, validated, and callable.
DONE: likelihood model (C1.x); E-step (threefry rxRmvn draw -> in-process solve
-> log p -> n_i/lnL via log-sum-exp, with sample storage); M-step conjugate
mu/Omega via joint Metropolis-Hastings + additive residual (Eq 17); R scaffold
(rpemControl, parameter classifier .rpemClassify, .rpemFit E-M loop); and the
`est="rpem"` dispatch (nlmixr2Est.rpem). `nlmixr2(model, data, est="rpem")` runs
end-to-end and, on well-identified data, matches a FOCEI fit on mu, Omega AND
add.sd. Six RPEM test files pass on top of origin/main.
REMAINING: full nlmixr2FitData (residuals/tables/SEs); numeric fixed-effect
update for non-additive residuals + non-mu-ref structural params; move the
iteration loop into C++; OpenMP via par_solve (task #7); then M2+ (mixtures, IOV,
censoring, multi-endpoint).

Scope (in build order):
- **First task DONE**: the dedicated RPEM likelihood model
  (`13-likelihood-model.md`) -- phi/eta supplied, no optimization. Built,
  compiles, and C1.x-verified at machine precision
  (`tests/testthat/test-rpem-llik-model.R`).
- E-step + M-step engine, K=1 (single Gaussian), **in C++/C in `src/rpem.cpp`**
  (D17), driving rxode2's C-level solve in-process with no per-iteration R
  round-trip.
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

## Git cadence (every small step)

Work happens in the `rpem` worktree (`~/src/nlmix2est-rpem`, branch `rpem`).
At the end of each small, self-contained step:

1. Commit early and often -- one focused commit per step (not one giant commit
   per milestone). Terse, past-tense messages; ASCII only.
2. `git fetch origin` and merge `origin/main` into `rpem`. Resolving small
   conflicts continuously keeps the branch close to upstream so the eventual
   merge back is cheap, rather than one large end-of-project reconciliation.
3. `git push` (`-u origin rpem` the first time) so work is backed up and visible.
4. If a change to rxode2 is needed, use a separate rxode2 worktree and apply the
   same commit/merge-main/push cadence there.

Never merge `rpem` back into `main` without an explicit request; the cadence
pulls main INTO rpem, not the reverse.

## Dependency order

```
M1 --> M2 --> (M3, M4 parallelizable) --> M5
```

Paper Table 1 full reproduction is blocked on M2 (K=2). Everything else
validates at M1.
