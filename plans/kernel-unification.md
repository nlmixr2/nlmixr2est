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

- [x] I1: R scaffolding -- `sharedInner="shared"` installs the focei inner
  (inert; `cfg$sharedInnerEnv`, shared==classic bit-for-bit). MERGED.
- [x] I2: C++ `saemSharedResid_(etaMat)` returning per-obs f/r/loglik via
  likInner0 at given etas. MERGED.
- [x] I3: equivalence at eta=0 vs an independent structural solve (f to ODE
  tol, r exact). MERGED.
- [x] I4: full-fit equivalence -- shared f reproduces a converged fit's IPRED
  (~2e-4), r exact. Add + prop error validated. MERGED.
- [x] G1: `saemSharedResidUpdate_(theta, omega, etaMat)` -- cheap in-loop inner
  re-parameterize (vaeInnerUpdateParCore) + residual; the per-iteration callable.
  Reproduces IPRED via the fast path. MERGED.
- [ ] G2 (NEXT, deepest/riskiest -- C++ SAEM-loop surgery):
  1. In `src/saem.cpp` under `sharedInner`, at the M-step: switch the global
     solve to the FOCEi inner (the fsaem pattern -- `vaeInnerSetup_`/
     `vaeInnerUpdateParCore` makes the inner current), call the
     `saemSharedResidUpdate_` logic at the conditional-mean etas
     (`mpost_phi`->etas) to get per-obs f/r, then `setupRx(fsaemSaemOpt, ...)`
     to RESTORE the SAEM (N*nmc) solve before the rest of the M-step (mirror
     `saem.cpp:716` and `:3652`).
  2. FIRST do it as a DIAGNOSTIC (compute shared f in-loop, assert
     max|f_saem - f_shared| ~ ODE tol each iteration) to prove the in-loop
     solve-switch is correct WITHOUT changing behavior.
  3. THEN estimate the residual params from the shared residuals (y - f) via the
     per-endpoint moment (like `warmStartResid`), replacing `arResk`'s SSR for
     the residual M-step -- under the flag only.
  4. Equivalence-gate vs classic on the eval set; keep the classic path until the
     default flips.
- [ ] G3: route `calc.2LL` -- careful (SAEM's -2LL is MARGINAL/GQ, the shared
  llikObs is CONDITIONAL; needs the GQ integrator to consume the shared per-obs
  density, not a drop-in swap).
- [ ] G4: flip default once equivalence holds across the eval set.

RISK NOTE: G2 is the first change that mutates the live SAEM kernel state
(solve switching inside the loop).  The full callable chain (I1-I4, G1) is proven
to reproduce SAEM's f/r, so the numerical risk is retired; G2's risk is
state-management (restoring the SAEM solve on every path incl. early returns).
