# RPEM in nlmixr2est -- Overview

Implements the Randomized Parametric Expectation Maximization (RPEM) algorithm
of Chen et al. 2024 (CPT:PSP 13:759-780) as a first-class estimation method in
`nlmixr2est`, dispatched via `est = "rpem"`.

RPEM is a Monte Carlo parametric EM (MCPEM) method. It differs from other MCPEM
engines by using an *unbiased* Metropolis-Hastings estimator in **both** the
E-step and the M-step, sampling the discrete labels (subject `i`, mixture
component `k`) and the continuous parameter vector `theta_i` in one randomized
chain. This is claimed to be ~3-4x faster than SAEM and more accurate than
SAEM/QRPEM on the paper's voriconazole model.

## What RPEM computes (paper notation -> nlmixr2 mapping)

Two-stage nonlinear mixed-effects mixture model:

- Stage 1 (data): `Y_i | theta_i ~ N(h_i(theta_i), G_i(theta_i))`, where `h_i`
  is the structural (PK/PD) prediction and `G_i` is the residual covariance.
- Stage 2 (population): `theta_i ~ sum_k w^(k) N(mu^(k), Sigma^(k))`, a
  K-component Gaussian mixture over subject parameters.

Paper's split of `theta_i` maps to nlmixr2 parameter classes:

| Paper term | Meaning | nlmixr2 class |
|---|---|---|
| `alpha_i` | random effects with between-subject variability | mu-referenced pop param + `eta` |
| `mu^(k)`, `Sigma^(k)` | mixture means/covariances | fixed-effect `theta` (typical value) + `omega` |
| `zeta_i` | single-Gaussian (non-mixed) random component | mu-referenced param, single component |
| `beta` (e.g. `sigma^2`) | pure fixed effects, no BSV | residual-error params + non-mu-ref structural params |
| `w^(k)` | mixture weights | mixture weight params |

Key consequence: nlmixr2 residual-error parameters and any non-mu-referenced
structural/covariate coefficients are the paper's "fixed effects" -- they are
NOT updated by the Gaussian conjugate formulas; they are estimated numerically
in the M-step (see `05-mstep.md`, `03-parameter-space.md`).

## Core algorithm (one iteration)

1. **E-step** (ODE-bound): for each subject `i` and component `k`, draw
   `m_Gauss` samples `theta_i ~ N(mu^(k), Sigma^(k))`, solve the model, and
   average the data likelihood `p(Y_i | theta_i)` to get `n_ik` (Eq 23-24).
   Form `N_i = sum_k w^(k) n_ik` (Eq 22) and `tau_i(k) = w^(k) n_ik / N_i`.
   **Store the per-sample likelihoods for reuse in the M-step.**
2. **M-step** (cheap, no new ODE solves): run a joint Metropolis-Hastings chain
   over `(i, k, theta_i)` targeting `pi(s) = g_ik(theta_i)/N` (Eq 30), with
   acceptance ratio Eq 32. From the accepted samples, form the population
   updates (Eq 15-21) for `mu^(k)`, `Sigma^(k)`, `w^(k)`; update fixed effects
   (residual error, non-mu-ref params) numerically.
3. **Objective**: `lnL = sum_i ln(N_i)` (Eq 26).
4. **Stop**: least-squares slope of the last 30 `lnL` values first turns
   non-positive (see `10-stopping-control.md`).

## Reuse map (do not rebuild what exists)

- Individual data likelihood `p(Y_i | theta_i)` -> **a dedicated RPEM likelihood
  model**, built as a hybrid of `nlmModel0` (direct, no-optimization evaluation)
  and `foceiModel0ll` (THETA+ETA parameterization). Etas are supplied per sample,
  never optimized. Reuses rxode2's `rxGetDistributionFoceiLines` /
  `rxCombineErrorLines` / `rxPredLlik` machinery, so multiple endpoints, all
  error structures, autocorrelation, and censoring come "for free" (see
  `13-likelihood-model.md`, `04-estep.md`).
- Population ODE batch -> **rxode2 `par_solve`** (threaded).
- Numeric M-step optimization of fixed effects -> **nlm-family optimizer**
  (`src/nlm.cpp`, `R/nlm.R`).
- Mixture UI/control representation -> **SAEM plumbing** as one candidate
  backend (see `07-mixtures.md`).
- Pre-processing -> existing **hook system** (`R/hook.R`).
- Fit object, residuals, predictions -> **nlmixr2FitData**, `addCwres()`,
  `addNpde()`.

## Spec index

- `01-decisions.md` -- decision log from design interview (authoritative)
- `02-architecture.md` -- dispatch, files, data flow, reuse boundaries
- `03-parameter-space.md` -- mu-ref, transforms, fixed effects, truncation
- `04-estep.md` -- sampling, ODE batch, storage, `m_Gauss`
- `05-mstep.md` -- joint MH, conjugate vs numeric (hybrid) updates
- `06-parallelization.md` -- OpenMP + par_solve, RNG, determinism
- `07-mixtures.md` -- K>1: two backends, benchmark-and-choose
- `08-standard-errors.md` -- Fisher-score SEs (in v1)
- `09-fit-object.md` -- fit integration, residuals, predictions, VPC
- `10-stopping-control.md` -- convergence rule, `rpemControl()` options
- `11-validation.md` -- correctness bars and test plan
- `12-roadmap.md` -- milestones and cut lines
- `13-likelihood-model.md` -- **foundational**: the dedicated phi/eta-supplied,
  no-optimization likelihood model (read with `04`)
- `14-evaluation-criteria.md` -- precise, measurable acceptance gates
  (blocking vs tracked) across accuracy, numerics, speed, parity, code quality
