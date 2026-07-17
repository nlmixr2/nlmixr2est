# SAEM kernel unification onto the shared inner driver (sharedInner="shared")

Prototype-and-checkpoint (user-agreed).  Goal: route SAEM's residual /
prediction / likelihood computation through the shared FOCEI inner driver
`likInner0` instead of the duplicated per-endpoint machinery, gated by
`saemControl(sharedInner="shared")` (default "classic" stays until proven).

## What is duplicated today (removal targets, under the flag)

- C++ `src/saem.cpp`: `res_mod` per-endpoint residual codes, `arResk`
  (SSR accumulation), `get_resMat`/`get_resInfo`, and the residual M-step.
- R `R/saemRxUiGet.R`: the `res_mod` switch-statement plumbing.
- The marginal `-2LL` (`calc.2LL`, `R/saem_fit_aux.R`) consumes `resMat`.
- Post-fit residuals (CWRES/IPRED/NPDE) ALREADY go through the focei bridge
  (`.calcCwres`/`predOnly`), so only the estimation-time path is duplicated.

## The shared driver already exists and is proven

- `likInner0(eta, id)` (`src/inner.cpp:1267`): solves the `rxInner` model at
  `eta` for subject `id`; fills `fInd->llik` (objective), `fInd->llikObs[]`
  (per-obs conditional log-density), and f/r read via
  `grabRFmatFromInner(id)` (`inner.cpp:1065`, cols f, r on the tbs scale).
- `npResidELS(postEta)` (`inner.cpp`) is ESSENTIALLY the shared-inner SAEM
  residual driver already: it loops subjects, calls `npEvalCondLik` ->
  `likInner0`, reads f/r via `grabRFmatFromInner`, and sums the ELS normal
  negative log-likelihood `0.5*(f-dv)^2/r + 0.5*log(r) + const` -- endpoint- and
  error-model-agnostic (the residual model is baked into `rxInner`'s `rx_r_`).
- `vaeInnerSetup_(env)` (`inner.cpp:9815`) sets up the FOCEI inner (op_focei +
  rxInner + a per-subject solve) ONCE with no outer optimizer.  `fsaem`
  (`R/fsaemInner.R:.fsaemInstallStep`) already installs this alongside SAEM and
  calls into it per subject.

## The N*nmc-vs-N reconciliation

Classic SAEM's `user_fn(phiM,...)` solves `N*nmc` rows (all MCMC chains).
`likInner0` is per-subject (`N`).  The shared path evaluates the residual /
likelihood at the per-subject CONDITIONAL MEAN eta (`mpost_phi` -> etas), i.e. N
evaluations, not N*nmc -- the MCMC chains are only needed for the stochastic
suff-stat accumulation, which the shared path does NOT change (Omega, mixture,
etc. stay classic).  So "shared" swaps only the f/r/loglik computation, at the
phi means, via `likInner0`.

## Prototype scope (ONE error model first)

1. R: when `sharedInner=="shared"`, install the focei inner in the SAEM path
   (mirror `.fsaemInstallStep`: build `foceiOptEnv`, `.foceiPreProcessData`,
   `vaeInnerSetup_`), storing the env in `cfg$sharedInnerEnv`.  Re-set it at the
   final theta/omega via `.fsaemInnerUpdate` before the residual read.
2. C++: add `saemSharedResid_(etaMat, cores)` = a thin wrapper over the
   `npResidELS`/`grabRFmatFromInner` pattern that returns, per subject/obs, the
   f, r, and per-obs loglik at the given etas (the SAEM conditional-mean etas).
3. Wire the SAEM finalization so that under "shared" the reported per-obs f/r and
   the residual objective come from `saemSharedResid_`, and add an equivalence
   assertion vs the classic `resMat`/`arResk` values.
4. Restrict the prototype to a single additive-error continuous endpoint
   (`.fsaemSupported`-style gate); broaden after the checkpoint.

## Equivalence gate (before generalizing)

On `tools/saemRefactorEval.R` models (add first, then prop/combined), assert the
shared path's per-subject f, r, residual objective, and CWRES match the classic
path within a stated tolerance; confirm the shared path actually ran (not a
silent fallback).  CHECKPOINT with the user before extending to
multi-endpoint / censoring / TBS and before retiring `arResk`.

## Order of increments (commit each)

- [ ] I1: R scaffolding -- `sharedInner="shared"` installs the focei inner
  (inert; sets up `cfg$sharedInnerEnv`, no behavior change).
- [ ] I2: C++ `saemSharedResid_` returning per-subject f/r/loglik at given etas.
- [ ] I3: equivalence test (shared vs classic f/r/objf on the add model).
- [ ] I4: checkpoint; then broaden error models + retire `arResk` under the flag.
