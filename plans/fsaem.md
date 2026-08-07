# fsaem: fast SAEM (f-SAEM) for nlmixr2est

Implementation plan for a new SAEM extension based on the nlme-IMH / f-SAEM
algorithm of Karimi, Lavielle and Moulines, "f-SAEM: A fast stochastic
approximation of the EM algorithm for nonlinear mixed effects models",
Computational Statistics and Data Analysis 141 (2020) 123-138
(doi:10.1016/j.csda.2019.07.001). Paper PDF: `~/src/private/fsaem.pdf`.

Worktree: `~/src/nlmixr2est-fsaem` (branch `fsaem`). Base: `main` @ f1147783.

## 1. Goal

Add a fast variant of SAEM whose MCMC simulation step samples each subject's
random effects from an independent Metropolis-Hastings (IMH) proposal centered
at the per-subject MAP estimate, with a Laplace/linearization covariance. This
converges in far fewer SAEM iterations than the current random-walk-Metropolis
(RWM) SAEM. Exposed two ways:

- `saemControl(fast=TRUE, ...)` -- the real switch and all tuning knobs.
- `est="fsaem"` -- a thin sugar est method that forces `fast=TRUE`, exactly as
  `mufocei` is sugar over `focei` (`R/mufocei.R`).

Delivered in two phases (see Section 6). Phase 1 is a faithful reproduction of
both numerical examples in the paper; phase 2 hardens it into a production
accelerator across model types.

## 2. Verified decisions (from interview -- do not silently revisit)

1. Goal is BOTH a faithful paper reproduction AND, later, a production
   accelerator. Phased: reproduction is the phase-1 acceptance gate.
2. Proposal covariance: implement BOTH variants with auto-selection --
   - Jacobian linearization (paper Eq 17) for continuous data:
     `Gamma_i = (J' J / sigma^2 + Omega^-1)^-1`.
     General heteroscedastic form (Eq 19): `Gamma_i = (J' Sigma^-1 J + Omega^-1)^-1`.
   - Hessian / Laplace (paper Eq 13) for the general case (incl. non-continuous):
     `Gamma_i = (-H + Omega^-1)^-1`, `H` = Hessian of `log p(y_i | psi_i)` at MAP.
   The Hessian is taken from the FOCEI inner likelihood machinery. This same
   coupling is extended so that the foce / focep / focei inner likelihoods can
   be evaluated inside SAEM -- giving full likelihood coverage in SAEM as a
   side benefit. Reuse `likInner0` / `lpInner` (`src/inner.cpp:1165` / `:1700`).
3. All three kernel schedules are built in and selectable, documented in the
   roxygen2 for `saemControl()`:
   - `"firstN"`  -- IMH for the first N iterations, then revert to RWM (paper).
   - `"throughout"` -- IMH every iteration for the whole run.
   - `"additive"` -- IMH appended alongside the RWM kernels every iteration.
4. MAP estimate and its covariance are computed by REUSING the FOCEI inner
   optimization (`inner.cpp`), not a fresh optimizer in `saem.cpp`. Precedent:
   `vaeInnerSetup_` (`src/inner.cpp:8702`) already sets up the FOCEI inner
   problem once, with no outer optimizer, and evaluates it per subject (and per
   mixture component) via `likInner0`/`lpInner` in a parallel OpenMP loop.
5. Phase-1 acceptance gate = reproduce BOTH paper examples: continuous 1-cmt
   warfarin PK (Jacobian path) AND Weibull repeated-events TTE (Hessian path).
   So the FOCEI-inner Hessian coupling is in scope for phase 1.
6. Control API = flat `fast=` flag plus separate discrete args on
   `saemControl()` (matching the existing `nmc` / `nu` style), see Section 5.
7. Defaults when `fast=TRUE`: `firstN` schedule (N ~ 20) with `auto` covariance
   (Jacobian for continuous endpoints, Hessian otherwise), then revert to RWM.
8. `est="fsaem"` is PURE sugar: `fsaemControl()` inherits `saemControl()` with
   `fast` forced TRUE and un-overridable; no other default changes.
9. The fast kernel is integrated INTO the real SAEM C++ loop, not a parallel
   silo: it must run in both the plain-SAEM MCMC path (`do_mcmc`) AND the
   mixture path (`do_mcmc_msaem`, MSAEM).  For the non-fast (tail) iterations of
   the `firstN` schedule -- and whenever `fast=FALSE` -- the loop degrades to
   the existing plain-SAEM random-walk kernels bit-for-bit.  (User directive,
   2026-07-10.)

## 3. Algorithm background

Standard SAEM (paper Alg 3) alternates: (1) simulate `psi_i ~ p(psi_i | y_i;
theta)` by M MCMC steps, (2) stochastic-approximation update of the complete-data
sufficient statistics, (3) maximize. The current code implements the simulation
step with three kernels in `do_mcmc` (`src/saem.cpp:3004`):

- method 1: independent proposal from the prior marginal `N(mprior, Gamma)`
  (MH ratio uses only the observation term, prior cancels).
- method 2: block random-walk over all phi.
- method 3: component-wise random-walk.

f-SAEM adds the nlme-IMH kernel (paper Alg 2). For subject `i`, given current
chain state:

1. Compute the MAP: `psi_hat_i = argmax p(psi_i | y_i)` (paper Eq 12/22). This is
   exactly the FOCEI inner problem (minimize the individual objective).
2. Compute the proposal covariance `Gamma_i` (Eq 13 Hessian, or Eq 17 Jacobian
   for continuous data; Eq 19 for heteroscedastic error).
3. Propose a candidate `psi_c ~ N(psi_hat_i, Gamma_i)`.
4. Accept with the full IMH ratio (Eq 23):
   `alpha = [p(psi_c|y) q(psi_prev|psi_hat)] / [p(psi_prev|y) q(psi_c|psi_hat)]`
   where `q` is the Gaussian proposal density (does NOT cancel like the prior
   proposal in method 1 -- it must be evaluated explicitly).

For transformed parameters (Remark 5, Eq 20-21) the proposal is formed in the
transformed (normal) space `phi = u(psi)` and mapped back. SAEM already works in
a transformed phi space, which aligns with this.

Practical recipe from the paper: M ~ 10 MCMC steps per SAEM iteration; run the
IMH kernel only during the first ~20 SAEM iterations, then revert to RWM (Kuhn &
Lavielle 2004: early iterations need only an approximate posterior). Result:
convergence in <10 iterations (PK) / ~30 (TTE) vs ~50+ for standard SAEM.

## 4. Architecture

### 4.1 Two model representations to bridge

SAEM works in phi-space: `phiM` rows are per-subject/per-chain individual
parameters, `mprior_phiM` the prior means, `Gamma2_phi1` the prior covariance
(`set_mcmcphi`, `src/saem.cpp:1138`). FOCEI inner works in eta-space around a
theta. The bridge (mu-referencing) between `phi` and `theta + eta` must be made
explicit; the MAP and `Gamma_i` produced by the inner solve (eta-space) have to
be mapped into the phi-space the SAEM chain lives in. This mapping is the single
most important correctness risk -- specced as its own milestone (P1-2).

### 4.2 Reusing the FOCEI inner (the `vaeInnerSetup_` pattern)

fsaem sets up a FOCEI inner problem once, alongside the SAEM solve, with no outer
optimizer -- following `vaeInnerSetup_`/`vaeInnerEval` (`src/inner.cpp:8702+`).
Per SAEM iteration (while the IMH kernel is active) and per subject:

- `innerOpt`/`innerOpt1` to get the MAP eta (`src/inner.cpp:2110`).
- inner Hessian (already computed for FOCEI cov) -> `Gamma_i` via Eq 13, or the
  Jacobian/sensitivities -> `Gamma_i` via Eq 17/19.
- feed `psi_hat_i` and `Gamma_i` into the new SAEM kernel.

The inner theta is set from the current SAEM parameter estimate each iteration
(the analogue of FOCEI's `updateTheta`).

### 4.3 New SAEM kernel (method 4) -- integrated, mixture-inclusive, degrading

Add a new IMH kernel to the SAEM loop implementing steps 3-4 above: draw from
`N(psi_hat_i, Gamma_i)`, evaluate the full IMH acceptance ratio including the
Gaussian proposal density both directions. Per decision 9 it is wired into BOTH
kernel-dispatch sites:

- plain SAEM: alongside `do_mcmc` (`src/saem.cpp:1715+`);
- mixture SAEM: alongside `do_mcmc_msaem` (`src/saem.cpp:1181+`).

Degrade rule: on iterations where the fast kernel is not active (the tail of
`firstN`, or `fast=FALSE`) the loop calls exactly the existing random-walk
kernels, unchanged -- a regression fixture guards bit-identical plain-SAEM
output. Precompute the per-subject `psi_hat_i` / Cholesky of `Gamma_i` once per
active iteration; reuse across the M chain steps.

### 4.4 Orchestration: function-table swap on a shared allocation (chosen)

The validated proposal core (`fsaemInnerMap_`, done) runs while the FOCEi inner
model is the active global `rx` solve. Inside `saem_fit` the active global solve
is the SAEM model's own (`setupRx` + `rxInner`), used by `user_function`. Both
read `getRxSolve_()`, so they cannot be active simultaneously.

Chosen mechanism (user directive, 2026-07-10) -- do NOT swap the global
`rx_solve` pointer or add an rxode2 setter API. Instead, allocate ONE solve
sized for the LARGER problem (the inner/sensitivity-augmented model, a superset
of the SAEM prediction model) and swap the compiled-model FUNCTION TABLES
(`calc_lhs`/`update_inis`, the `rxUpdateFuns` / existing `saem_lhs`/`saem_inis`
pattern) between the SAEM model and the inner model on that shared allocation.
Per active SAEM iteration, in C++ segments:

1. Swap function table -> inner model; `updateTheta`/omega to the current SAEM
   estimate; run the per-subject MAP + `Gamma_i`; COPY the small results out
   (eta_hat_i: neta each; chol(Gamma_i): neta x neta each) into fsaem-owned
   storage.  The big shared solve memory is about to be overwritten.
2. Swap function table -> SAEM model; run the IMH kernel via `user_function`,
   using the copied-out proposals.

The first segment's solve buffers are gone after the swap (shared memory) --
hence the copy-out. Memory-efficient and API-change-free.

CORRECT mechanism (the focei precedent -- an earlier "cannot coexist" analysis
was WRONG; it mis-set-up TWO separate allocations). focei runs two models on ONE
shared allocation without any re-pointing:
- `innerOde(id)` / `predOde(id)` (inner.cpp:59-60) both call `ind_solve(rx, ...)`
  on the SAME `rx`, passing DIFFERENT model function pointers (`rxInner.dydt` vs
  `rxPred.dydt`). The ODE integrator uses the pointers handed to it, not the
  global.
- The compiled functions' `_solveData = _getRxSolve_()` (codegen2.h:614) returns
  the ONE shared `rx` for both models -- consistent, no clobber, because there is
  a single allocation.
- `IndNeqOverrideGuard` + `op_focei.predNeq` reconcile the different equation
  counts (the pred model has fewer states than the sensitivity-augmented inner).

So fsaem does NOT create a second solve or re-point. It sets up ONE allocation
sized for the LARGER problem -- the inner (sensitivity-augmented) model over the
SAEM design -- and runs both the fast kernel and the SAEM prediction on it via
function-table selection + neq override, exactly like focei's inner/pred. The
inner allocation is per physical subject (N); the nmc chains are handled by
evaluating `likInner0(eta, id)` at each chain's eta on those N slots (the kernel
already does this), so no allocation growth for chains.

Concrete P1-3 integration (validated cores: `fsaemInnerMap_`, `fsaemImhKernel_`):
set up the FOCEi inner as the shared solve for the fsaem fit (the `.fsaemInner*`
path, already built). Each `firstN` iteration, in C++: (1) refresh the inner at
the current SAEM estimate (theta from Plambda, omega from Gamma2); (2)
`fsaemInnerMap_` -> per-subject MAP + `Gamma_i`; `fsaemImhKernel_` -> accepted
etas for the chains; (3) write etas into `phiM` (phi = mprior + eta, mu-ref/cov
aware) and read the SAEM prediction for the M-step off the SAME shared solve via
the SAEM function table (predNeq override). iter >= N (and `fast=FALSE`): the
existing `do_mcmc` path -- degrade, guarded bit-identical. Mirror into
`do_mcmc_msaem` for mixtures. Open detail to settle in code: whether SAEM's
`user_function` (`par_solve` over all N*nmc) reads predictions off the shared
inner solve or keeps the SAEM model's own function table on that allocation.

### 4.4 R plumbing

- `R/saemControl.R`: new args + validation + `mcmc`-list packaging so the C++
  layer receives the schedule/cov/iter/lik settings. Roxygen documents each.
- `R/fsaem.R`: `nlmixr2Est.fsaem` (sugar) + `fsaemControl()` (inherits
  saemControl, forces `fast=TRUE`), `nmObjHandleControlObject.fsaemControl`,
  attrs mirroring `R/mufocei.R`.
- `R/nlmixr2Est.R`: register `fsaem` in the dispatch table.
- `src/init.c`: MANUAL update if any `[[Rcpp::export]]` arity changes
  (`compileAttributes` + hand-edit `callMethods[]`), per CLAUDE.md.

## 5. Control API (saemControl)

Flat args, defaulting to today's behavior when `fast=FALSE`:

- `fast = FALSE` -- master switch. `TRUE` enables the IMH kernel.
- `fastKernel = c("firstN", "throughout", "additive")` -- schedule (Section 2.3).
- `fastCov = c("auto", "jacobian", "hessian")` -- proposal covariance; `auto` =
  Jacobian for continuous endpoints, Hessian otherwise.
- `fastIter = 20L` -- N for the `firstN` schedule (ignored otherwise).
- `fastLik = c("focei", "foce", "focep")` -- which inner likelihood the Hessian
  path uses (also selects the SAEM likelihood when reported via that path).
- `fastNmc` (optional, phase 2) -- M chain steps for the IMH kernel if it should
  differ from `nmc`.

`fsaemControl()` = `saemControl(fast=TRUE, ...)` with `fast` un-overridable.

All decisions above are locked; see Section 2. New knobs must document sensible
defaults and be listed in the `saemControl()` roxygen `@param` block.

## 6. Phased implementation (small, compartmentalized, commit per checkpoint)

Commit after each milestone (project convention: commit per completed
checkpoint, not one batch). Merge `origin/main` before any push.

### Phase 0 -- scaffolding (no behavior change)
- P0-1: `saemControl()` args + validation + `mcmc`-list packaging; roxygen.
  Gate: `saemControl(fast=TRUE)` round-trips; `fast=FALSE` byte-identical to now.
- P0-2: `R/fsaem.R` sugar + dispatch registration; `est="fsaem"` runs today's
  RWM SAEM (fast path still a no-op). Gate: a trivial model fits via `fsaem`.

### Phase 1 -- reproduce the paper (acceptance gate)
- P1-1/P1-2 (DONE, commit 21dde70f): FOCEi-inner reuse producing the per-subject
  MAP + proposal precision `H = Gamma_i^-1` (`fsaemInnerMap_` in src/inner.cpp;
  R/fsaemInner.R).  Validated on theo_sd: MAP stationary to grad ~1e-6 and `H`
  matches the independent Eq-17 information to ~4e-5 relative
  (tests/testthat/test-fsaem-inner.R).  This is the standalone proposal core; it
  is NOT yet wired into the SAEM loop.
- P1-3: integrate the IMH kernel into the SAEM loop with the solve-swap
  orchestration (Section 4.4), the new IMH kernel (Section 4.3), the `firstN`
  schedule + Jacobian cov, degrading to the plain RWM kernels for the tail.
  Wire BOTH `do_mcmc` and `do_mcmc_msaem` (mixture) per decision 9.  Gate:
  warfarin PK converges to the saemix/FOCEI reference in far fewer iterations;
  parHist/counters show the IMH kernel actually fired (repo policy: fast-path
  tests must assert the mechanism was used); `fast=FALSE` bit-identical.
- P1-4: Hessian cov path via inner Hessian + `throughout`/`additive` schedules +
  `auto` selection; mixture (MSAEM) fast path exercised.  Gate: Weibull
  repeated-events TTE reproduces the paper's faster convergence.
- P1-5: tests + `NEWS.md` + roxygen `devtools::document()`; `R CMD check`.
  Gate: acceptance = both paper examples reproduced (Section 7).

### Phase 2 -- production hardening
- P2-1: robustness -- PD-safeguard `Gamma_i` (reuse `nearPD.cpp`), non-finite
  MAP fallback to RWM for that subject/iteration (loud, not silent).
- P2-2: broaden coverage -- multiple endpoints, mu-referencing edge cases, IOV,
  mixtures (MSAEM path), covariates.
- P2-3: full likelihood coverage -- surface foce/focep/focei likelihoods
  computed through the inner path as a reported SAEM objective.
- P2-4: performance -- parallelize the per-subject MAP/cov (OpenMP, thread count
  from rxode2 as elsewhere); cache Cholesky; avoid recompute in `firstN` tail.
- P2-5: docs/vignette + `fsaem` benchmark against `saem`.

## 7. Evaluation criteria (definition of high quality)

Correctness
- Jacobian `Gamma_i` matches paper Eq 17/19 recomputed independently (~1e-6) on
  warfarin.
- Hessian `Gamma_i` matches the FOCEI inner Hessian used for FOCEI covariance.
- IMH acceptance ratio includes the proposal density both directions (a linear
  model must accept with probability 1, paper Remark 3 -- an exact unit test).
- Final MLE from `fsaem` equals `saem` (same seed/data) within SAEM stochastic
  tolerance on warfarin PK and Weibull TTE; parameter values not systematically
  biased across a small Monte Carlo replicate set.

Reproduction (phase-1 acceptance gate)
- Warfarin PK (32 subjects, 1-cmt oral, lognormal ka/V/k): `fast=TRUE` reaches
  the reference estimate in substantially fewer iterations than `fast=FALSE`
  (paper: <10 vs ~50), monotone-ish convergence of the estimate sequence.
- Weibull repeated-events TTE: analogous acceleration via the Hessian path.

Mechanism-is-used (repo policy)
- Tests assert the IMH kernel fired (parHist type / counter / kernel-usage flag),
  not merely that results match -- silent fallback to RWM must be observable.
- Schedule selection is exercised: `firstN` reverts to RWM at N; `throughout`
  and `additive` run the kernel every iteration (assert per-iteration usage).

Robustness / no regression
- `fast=FALSE` is byte-for-byte the current SAEM (no numeric drift) -- a
  regression guard on an existing SAEM fixture.
- Non-finite MAP / non-PD `Gamma_i` handled with a loud, documented fallback,
  never a crash or a silently wrong proposal.

Engineering hygiene
- Main numerical kernel in C++ (`do_mcmc` method 4 + inner reuse); R layer thin.
- `src/init.c` hand-maintained if any export arity changes (CLAUDE.md).
- ASCII-only source/docs; camelCase R; terse single-bullet `NEWS.md`; short
  roxygen; `devtools::document()` + `R CMD check` clean.
- Commit per checkpoint; `origin/main` merged before every push.

## 8. Risks and open questions

- phi<->eta / mu-referencing mapping (P1-2) is the top correctness risk; validate
  numerically before building the kernel on top.
- Running a full FOCEI inner solve every early SAEM iteration is costly; the
  `firstN` default bounds this. Measure overhead vs iteration-count savings.
- Coupling SAEM to the FOCEI inner means both model representations must be
  compiled/consistent for one fit -- watch the 64-slot rxode2 model cache and
  concurrent-install pitfalls (see project memory).
- TTE MAP + Hessian through the inner path for non-continuous endpoints needs the
  inner likelihood to already support the endpoint type; confirm before P1-4.
- Interaction with the SA covariance phase (`nSaCov`) and MSAEM mixtures is out
  of phase-1 scope; keep the kernel dispatch guarded so those paths are
  unaffected when `fast=TRUE` is combined with them (or error clearly).

## 9. General-likelihood engine hook (Weibull TTE and beyond)

Phase-1 gate example 2 (Weibull TTE) is *not* reachable through SAEM's own
kernel: `nmGetDistributionSaemLines` implements only `.norm`, and the C++
E-step (`do_mcmc` + the `distribution==1/2/3` DYF blocks in `saem.cpp`) hard-code
normal/poisson/binomial observation losses. But fsaem's IMH kernel evaluates the
observation likelihood through the FOCEi inner (`likInner0`), which is fully
general (TTE, ordinal, dpois/dbinom, custom `llik`). So a general likelihood can
be estimated by fsaem (never plain saem) by driving the whole E-step off the
inner and never calling `do_mcmc`.

Key structural fact (verified in `saem.cpp` non-mixture M-step): the M-step
splits into (a) phi sufficient statistics `Statphi11/12/01/02` -> mu-ref means +
Omega, which are **likelihood-independent** (Gaussian phi hierarchy), and (b) the
residual `statr`->`sigma`/`ares`/`bres` update, which is the **only**
normal-specific M-step piece. For a general likelihood every population
parameter lives in phi (a mu-ref fixed effect, with or without a random effect,
phi1/phi0), so the M-step reduces to (a) alone -- no residual optimizer.

### DONE (2026-07-11) -- implemented the saemix way, for BOTH saem and fsaem

The first cut bypassed `do_mcmc` and drove the E-step off the inner only; that
was too invasive and broke the saem core.  The correct design (from `~/src/saemix`
`compute.LLy`) is far smaller: a `modeltype="likelihood"` model returns the
per-observation log-likelihood as its prediction, the observation loss is just
`DYF = -ll`, and the *standard RWM kernels run unchanged*.  So plain `saem` -- not
only `fsaem` -- fits `ll()` endpoints; the f-SAEM inner kernel is an optional
accelerator on top.  Implemented: `distribution==4` arm (`DYF/DYFm = -fk`) in the
main-loop DYF block and both `do_mcmc` switches; residual M-step pieces
(`statr`/`sigma`/residual FIM row/residual-param update) skipped for
`distribution==4`; the phi sufficient-statistics M-step is untouched.  R: `ll()`
endpoints route through the FOCEi line generator (`rx_pred_ ~ <ll>`), config
`distribution=4`, zeroed residual bookkeeping, and a fix for a spurious
`NA=THETA[]` line the saem model emitted when an endpoint had no residual theta.
Validated: exponential and Weibull TTE recover the population parameters for saem
and fsaem; `test-saem-loglik.R`.

### Design (determined by the math; no separate residual/error step)

- New distribution code `distribution == 4` ("inner/general"): C++ skips the DYF
  build and `do_mcmc`/`do_mcmc_msaem`, requires `fsaemStepFn` (error if fast is
  off), and forces the IMH kernel every iteration.
- fsaem forces `fastKernel="throughout"` for non-normal endpoints (a general
  likelihood must take the inner path on every iteration -- user directive).
- M-step: keep `Statphi11/12/01/02` (means + Omega) exactly as-is; skip `statr`,
  `sigma`, and the `d1_logsigma2`/`nb_param-1` FIM residual contributions when
  `distribution==4` (no residual parameter exists).
- `fsave`/predictions are not needed for the phi-stat M-step; for `distribution==4`
  skip the SAEM `user_fn` residual solve (the inner already produced the E-step).

### Phases

- **G1 (R model translation)**: let `nmGetDistributionSaemLines` accept the
  general endpoint for fsaem -- emit a benign pred line so `saemModel` compiles;
  mark the config `distribution=4`; set a "no residual" `res.mod`/`ares`/`bres`
  structure (`nendpnt` residual bookkeeping must not assume add/prop). Detect
  non-normal endpoint in `.fsaemSupported` -> allow only when fast; else error
  with a clear "saem cannot fit this endpoint; use est=fsaem" message. Force
  `fastKernel="throughout"`.
- **G2 (C++ E-step bypass)**: in the non-mixture branch, gate the DYF blocks and
  `do_mcmc` behind `distribution != 4`; for `distribution==4` run only
  `fsaemImhStep(phiM)` (already the general E-step) + phi-stat accumulation.
- **G3 (C++ M-step)**: skip residual accumulation / `sigma` update / residual FIM
  rows for `distribution==4`; verify Plambda + Omega updates still converge.
- **G4 (inner for TTE)**: ensure the fsaem inner is built on the TTE model so
  `likInner0` returns the survival likelihood and `fsaemInnerMap_` uses the
  Hessian proposal (`fastCov` auto -> hessian for non-normal). Reuse the existing
  inner-setup path.
- **G5 (validate)**: reproduce the paper Weibull TTE example; assert fsaem
  recovers the scale/shape and that the kernel fires every iteration; plain saem
  errors on the same model. Add a regression test.

### Open question (pick default, revisit if wrong)

Scope = the *general* FOCEi likelihood surface (any inner-supported endpoint) not
TTE-only, since the inner handles them uniformly; Weibull TTE is the validation
gate. Endpoints whose inner likelihood needs the Hessian proposal use
`fastCov` auto -> hessian (already validated on continuous).

## Phase R -- restoration (issue #845, 2026-08-05)

`est="fsaem"` was removed in 76a539b8f because it showed no early-convergence
advantage over plain SAEM.  Issue #845 asks for it back now that the shared ODE
pool (`odeSwap`) is in place.  Worktree `~/src/nlmixr2est-fsaem`, branch
`feat/fsaem-restore`.

- **R1 (revert)**: revert 76a539b8f onto main.  `.saemGeneralLik` keeps its
  post-removal name and home -- plain saem fits `ll()` endpoints now, so the
  detector is shared rather than duplicated.  NEWS goes under 7.0.3 (7.0.2 is the
  in-flight CRAN submission).
- **R2 (odeSwap audit)**: the pool itself is INERT for fsaem -- `.odeSwapInfo()`
  after a fit shows one registered slot (`inner`), `overrideNeeded=FALSE`,
  `scratchNlhs=0`, `lhsWidthMismatchN=0`, `activeOverride` all -1.  What IS
  shared and was wrong is the event-sensitivity ("jump") shape: setting up the
  FOCEi inner points rxode2's globals at the inner, and the SAEM model -- which
  has none -- was solved under it.  Fixed with `odeSwapEsOff()` at every SAEM
  solve setup/restore plus an explicit `OdeSwapEsBatch(odeSlotInner)` around the
  MAP+IMH step.  Only reachable for a model with modeled dosing (lag/f/rate);
  a plain bolus model has `eventSensInfo` NULL, which is why every existing
  fsaem test was blind to it.
- **R3 (stale chain likelihood)**: only the general-likelihood branch rescored
  the chain after the IMH move.  For a normal endpoint the kernel is additive and
  the random walks still run, so `do_mcmc` was scoring proposals against the
  PRE-IMH `U_y`/`DYF`/`fsave`.  The DYF rebuild is now a lambda re-run after the
  move for every distribution.
- **R4 (diagnostics)**: `fsaemDiag_()`/`fsaemDiagReset_()` expose steps,
  proposals, acceptances, MAP failures and non-PD proposal covariances -- the
  only evidence a test has that the kernel fired rather than silently degrading.
- **R5 (benchmark)**: the removal benchmark (theo_sd 1-cmt) cannot discriminate.
  It is well identified, 20 burn-in iterations already reach the MLE, and the
  f-SAEM paper itself reports "no regression" (i.e. no gain) there -- both
  methods sit at SAEM's Monte-Carlo noise floor (RMSE ~0.01).  On the collinear
  2-cmt N=60 model issue #845 cites, fsaem is clearly ahead at reduced
  iterations.  The "converges faster" unit test is rewritten accordingly.

## Phase S -- issue #845 options + verifying the over-acceptance fix (2026-08-05)

Continues Phase R on the same branch (`feat/fsaem-restore`, PR #852).

- **S1 (acceptance instrumentation)**: `do_mcmc` carried NO acceptance signal --
  `rmcmc` is a fixed 0.5 and nothing adapts from a measured rate -- so the stale
  `U_y` defect had nothing that could have caught it.  Adds per-kernel
  proposal/acceptance counters, plus `nRescore`/`uYStaleMax` on the post-IMH
  rescore, exposed as `fit$saemDiag`; `fit$fsaemDiag` moves onto the fit env too
  (the impmap `odeSwapInfo` pattern) instead of a `:::` global a later fit
  overwrites.

  **The predicted symptom was wrong.**  Omitting the rescore does NOT inflate the
  acceptance rate: `do_mcmc` writes `U_y(ind)=Uc_y(ind)` on every acceptance, so
  a stale entry self-heals the first time that subject accepts anything and the
  damage is diluted to noise (method-2 acceptance 0.2369 without the rescore
  against 0.2344 for plain saem, with the kernel firing on 40 of 45 iterations).
  The real symptom is a shifted answer -- tka 0.4362 against 0.4620.  So the
  regression test asserts `nRescore == nStep` and `uYStaleMax > 0` (~22 nats),
  NOT an acceptance bound.
- **S2 (options carrier)**: `fsaemMapImh` is reached through two registered entry
  points with hand-maintained arities in `src/init.c`, so one `fsaemSetOpts_()`
  setter carries all four new tunables as file statics -- one `init.c` row
  instead of four rounds of `compileAttributes` + `init.c` edits.  Reset per fit
  alongside the counters.
- **S3 (`nu[4]`)**: the sweep count, previously hardcoded 5 in four places.
  Sweep count ONLY -- the kernel is still switched on by `fast=TRUE`.  Relaxes
  the three `len=3` assertions, and fixes the `mcmc=` passthrough branch, which
  validated `mcmc$nu` but never assigned it.  The covariate path could not honor
  a sweep count at all before this (it returns before `cfg$fsaemNsweep` is
  attached and its closure is called with 8 arguments).
- **S4 (`fastFallback`)**: `"prior"` proposes from `N(centre, Omega)` for a
  subject with no usable Gamma.  Omega comes from the inner's own omega inverse,
  so no entry point widens.  Also counts a fourth failure mode that was
  invisible: a finite Cholesky that will not itself invert.

  **Not covered end to end**: `H = J' Sigma^-1 J + Omega^-1` is PD by
  construction, so no healthy fit produces a bad Gamma -- measured 0 on the
  exponential-TTE general likelihood, on 2-obs-per-subject theo_sd, and on the
  same with `fastCov="hessian"`.  Driving the prior term to zero fails earlier in
  SAEM's own nearPD.  The test asserts the option is INERT.
- **S5/S6 (`fastMode`, `fastHRefresh`)**: chain-mean centring, and reusing the
  MAP + Gamma between refreshes.  Both are safe because an independent MH
  proposal need not equal the target; they cost acceptance, not correctness.
  `fastHRefresh` must NOT skip the per-iteration re-parameterization -- likInner0
  has to score at the current theta/omega -- and does not.  Measured (theo_sd, 20
  fast iterations, all the same MLE): default accRate 0.882 / reuse 0;
  chainMean 0.594; hRefresh=5 0.671 / reuse 16; both 0.458 / reuse 16.
