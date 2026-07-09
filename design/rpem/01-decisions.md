# RPEM -- Decision Log

Authoritative record of decisions from the design interview. Later specs must
conform to these; if a spec needs to deviate, change it here first.

| # | Decision | Choice | Notes |
|---|---|---|---|
| D1 | Primary goal | First-class production method; must also reproduce and benchmark speed/accuracy vs SAEM and FOCEI on known models | Correctness first, speed second |
| D2 | Likelihood machinery for `p(Y_i\|theta_i)` | Build a dedicated RPEM likelihood model -- a hybrid of `nlmModel0` (direct, no-optimization llik eval) and `foceiModel0ll` (THETA+ETA parameterization). Etas supplied per sample, never optimized. NOT `foceiControl(maxInnerOpt=0)` (wrong tool) | Reuses `rxGetDistributionFoceiLines`/`rxCombineErrorLines`/`rxPredLlik`, so still buys endpoints, error structures, censoring, autocorrelation. See `13-likelihood-model.md` |
| D3 | Parallelization | OpenMP + rxode2 `par_solve`; no MPI | Matches house style + CRAN portability |
| D4 | Parameter space | Mu-referenced + existing bounded transforms as the primary path, AND also support non-mu-referenced params | Non-mu-ref params become numerically-estimated fixed effects |
| D5 | Truncated Gaussian | Avoid via nlmixr2 transforms (sample etas on unconstrained scale); do not implement paper's truncation for the transformed path | Truncation handling only if raw-normal params are later required |
| D6 | No-BSV / fixed-effect estimation (residual error, non-mu-ref structural) | Hybrid: conjugate closed-form where available, numeric optimization (nlm) otherwise | General error structures need the numeric path |
| D7 | E-step Monte Carlo | Redraw fresh samples each iteration (paper-faithful "randomized"); expose `m_Gauss` and M-step MH trial count as controls | Store per-iteration likelihoods for M-step reuse |
| D8 | Mixtures (K>1) | Build an RPEM-specific mixture backend AND a SAEM/NONMEM-representation backend; benchmark for correctness; pick the better as default | See `07-mixtures.md` |
| D9 | Standard errors | Include Fisher-score SEs in the first increment | See `08-standard-errors.md` |
| D10 | Fit object | Full residual/prediction parity (nlmixr2FitData with CWRES/NPDE/IPRED, population + individual predictions) | See `09-fit-object.md` |
| D11 | First increment scope | Breadth-first likelihood core: all continuous error structures + multiple endpoints + full fit object + SEs, K=1. Defer K>1 mixtures, IOV, censoring to later specs | See `12-roadmap.md` |
| D12 | Correctness bars | All four: reproduce paper analytic case; match SAEM/FOCEI on an nlmixr2 model; recover simulation truth; benchmark speed vs SAEM/FOCEI | See `11-validation.md` |
| D13 | Benchmark models | Warfarin PK; Theophylline 1-cmt; Voriconazole ODE (paper); Paper 1-cmt analytic | Voriconazole is K=1; analytic is K=2 |
| D14 | Deliverable of design phase | Numbered compartmentalized specs under `design/rpem/`, iterated and committed | This directory |
| D15 | Git workflow | Work in the `rpem` worktree; commit often (one focused commit per small step); at each step end fetch + merge `origin/main` into `rpem` and push. Keeps the branch close to upstream for a cheap eventual merge | See `12-roadmap.md` "Git cadence" |

## Assumptions (house style; not explicitly asked -- flag if wrong)

- A1: `est = "rpem"`, control constructor `rpemControl()`, S3 handler
  `nlmixr2Est.rpem` in `R/rpem.R`, following `sharedControl()` for common opts.
- A2: RNG uses rxode2's thread-safe parallel generator (threefry) with a
  `seed`/`sim` control for reproducibility across threads.
- A3: C++ kernel in `src/rpem.cpp`; UI translation in `R/rpemRxUiGet*.R`;
  pre-processing reuses `R/hook.R`.
- A4: Autocorrelated residuals are treated as an rxode2 `llik` error-structure
  feature and land in the first increment alongside the other error structures
  (only censoring and IOV/mixtures are deferred).

## Open tensions to resolve during implementation

- T1: D2 (minimal core) vs D11 (breadth-first) -- reconciled by defining the
  first increment as breadth-first over error structures/endpoints but K=1 and
  no IOV/censoring. Confirm the exact cut line in `12-roadmap.md`.
- T2: Paper analytic model is K=2; a literal Table 1 reproduction therefore
  needs the mixture backend. A K=1 reduction is the M1 smoke test; full Table 1
  reproduction is an M2 acceptance item.
