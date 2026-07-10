# VAE-NLME estimation method for nlmixr2est -- implementation plan

Reference: Rohleff et al., "Redefining Parameter Estimation and Covariate
Selection via Variational Autoencoders: One Run Is All You Need", CPT:PSP 2025,
DOI 10.1002/psp4.70129. Paper PDF: ~/src/private/vae.pdf. Python original:
~/src/vae_nlme. R prototype (torch, not integrated): ~/src/vae-nlme-nlmixr2-main.

Work happens in git worktree ~/src/nlmixr2est-vae on branch `vae`.

## 1. Goal and definition of done

Add a first-class `est="vae"` estimation method to nlmixr2est. A VAE (LSTM
encoder over each subject's series -> Gaussian posterior q(eta|y); structural
ODE decoder) is trained on an ELBO / BICc-ELBO objective to estimate population
parameters AND select covariates in a single run.

Definition of done = full parity with both paper case studies:
- Case Study 1: theophylline single-dose (and multiple-dose), Table 1.
- Case Study 2: neonatal weight progression, 2425 subjects, Table 2 / Fig 4.

Built in small compartmentalized increments with a commit at every checkpoint;
full parity is the acceptance bar, not the first commit.

## 2. Locked decisions (from interview 2026-07-09)

| Decision | Choice |
|---|---|
| Deliverable | First-class `est="vae"` in nlmixr2est (nlmixr2Est dispatch, nlmixr2FitData) |
| Encoder + gradients | Native C++ LSTM with analytic backprop (RcppArmadillo/Eigen); no torch runtime dep |
| Decoder / ODE | rxode2 solver + its analytic sensitivity machinery for d(pred)/d(eta) |
| Covariate selection | R-level MIQP/selection solver (Suggests: highs/Rglpk/ompr), C++ hot loop feeds it sufficient stats |
| Standard errors | Linearization (FOCE-style) marginal-likelihood Hessian at convergence |
| torch | Test-only oracle (Suggests) to gradient-check the C++ encoder; never a runtime dep |
| Encoder architecture | Reverse-engineer compiled encoder.c as ground truth |

## 3. Reference facts already established

- Encoder is a **single-layer UNIDIRECTIONAL** LSTM (`encoder.c` shows
  `nn.LSTM(input_size, hidden_size, num_layers, batch_first)`, no
  `bidirectional`). The prototype README's "bidirectional" is incorrect.
- FC head: `Linear(h_last | covariates) -> [mu(z_dim) | log_sigma(z_dim) |
  L_offdiag(z_dim*(z_dim-1)/2)]`. L = tril, diag = exp(log_sigma).
  Reparam z = mu + L @ eps.
- Decoder: theophylline single-dose is a closed-form 1-cmt solution; neonates
  is an ODE solve (torchode). Under nlmixr2est BOTH go through rxode2.
- ELBO (per functions.py): loss = p(x|z) [gaussian nll] + DKL, with
  DKL = p(z) - q(z|x), z ~ N(z_pop = C beta, Omega). M-step closed-form
  updates z_pop, omega, a from EMA sufficient statistics s1..s4.
- Schedule: burn-in (encoder only, tiny KL) -> KL annealing (0.01->1) ->
  main EM with covariate selection -> EMA smoothing of z_pop.
- Covariate selection minimizes BICc-ELBO (L0 penalty on beta); reference uses
  Gurobi MIQP, prototype uses per-parameter subset enumeration.

## 4. Architecture / process split

Mirror the saem.cpp pattern: numerically heavy kernel in C++, thin orchestration
in R.

- C++ kernel (`src/vae*.cpp`): LSTM fwd/bwd (analytic), FC head, L build,
  reparam, ELBO + gradient, Adam step, closed-form M-step sufficient stats.
  Exposed as a stepper `.Call`ed per inner iteration.
- rxode2: compiled decoder ODE solve + analytic d(pred)/d(eta) at obs times,
  reusing the focei / foceiCovAnalytic sensitivity infrastructure.
- R orchestration (`R/vae.R`): EM loop, schedule, periodic R-level covariate
  selection solve, assembly of nlmixr2FitData.
- init.c updated by hand for every new `.Call` (per CLAUDE.md manual init.c rule
  + compileAttributes step).

## 5. Open technical questions to resolve in Phase 0/1 (not user decisions)

- Exact LSTM gate ordering and weight layout PyTorch uses (i,f,g,o) so C++
  backward matches the oracle bit-for-bit.
- Whether d(pred)/d(eta) from focei machinery is exposed per-observation in the
  form the ELBO gradient needs, or needs a thin adapter.
- Standardization of inputs (time, conc, covariates) the encoder expects -- must
  replicate reference scaling exactly for parity.
- EM boundary crossing cost: how often covariate selection runs and whether the
  R solve is cheap enough per outer iteration.

## 6. Phased plan (each phase = commit checkpoint)

**Phase 0 -- Reference capture & spec.**
- Fully reverse-engineer encoder.c/pop_parameter.c; write `plans/vae-encoder-spec.md`
  (exact layer config, init, gate order, FC layout, standardization).
- Generate golden fixtures: for a fixed seed + tiny input, capture encoder
  forward outputs, ELBO, and gradients (from prototype/torch) as `.rds`/CSV
  under tests/testthat/fixtures/vae/.
- Deliverable: spec doc + golden fixtures. No package code yet.

**Phase 1 -- Dispatch skeleton (R only).**
- `vaeControl()` + `nlmixr2Est.vae` + `vaeRxUiGet*` handlers; wire into
  rxPipeline, data prep, theta/eta ordering. Add to R/nlmixr2Est.R dispatch.
- Returns an explicit "not implemented" error but validates UI->rxode2 for a
  VAE model. Small; no math.

**Phase 2 -- C++ LSTM encoder (fwd + analytic bwd).**
- `src/vaeEncoder.cpp`: LSTM cell fwd/bwd, FC head, Cholesky L, reparam,
  analytic gradients. compileAttributes + init.c edit.
- Gradient-check vs torch oracle (Suggests) and finite differences against
  Phase 0 golden fixtures (tol ~1e-6).

**Phase 3 -- Decoder + analytic d(pred)/d(eta) via rxode2.**
- Solve the structural model in rxode2 and extract pred + per-obs eta-Jacobian
  from the focei sensitivity machinery.
- Validate on theophylline against the closed-form solution + FD.

**Phase 4 -- ELBO + Adam + closed-form M-step (Milestone A).**
- Assemble p(x|z), p(z), q(z|x); ELBO gradient chaining encoder bwd with decoder
  sensitivities. Adam in C++. Closed-form z_pop/omega/a from EMA stats.
- Schedule (burn-in, KL anneal, smoothing) in the thin R loop.
- MILESTONE A: theophylline single-dose pop params (no covsel) reproduce Table 1
  within tolerance. Tag/commit.

**Phase 5 -- Covariate selection (Milestone B).**
- Full covariate effect matrix beta; BICc-ELBO objective; R-level solver
  (highs/Rglpk/ompr, Suggests) with per-parameter enumeration fallback. Wire
  periodic selection into the EM loop.
- MILESTONE B: theophylline reproduces WT on ka and V, not ke.

**Phase 6 -- Fit object: SEs, OFV, BICc, diagnostics.**
- Linearization-Hessian covariance -> $parFixed SEs. OFV via linearization + IS;
  AIC/BIC/BICc. Populate nlmixr2FitData; EBEs from encoder posterior/MAP;
  addCwres/addNpde compatibility.

**Phase 7 -- Neonates parity (Milestone C = done).**
- Zero-order production + first-order elimination model; 5 params, 5 covariates;
  multiple obs/subject; 2425 subjects with rxode2 threading.
- Also validate theophylline multiple-dose.
- MILESTONE C: reproduce Table 2 (-2LL ~146370 lin; BICc ~146566) and Fig 4 VAE
  covariate set.

**Phase 8 -- Docs, tests, CRAN hygiene, merge/push.**
- Roxygen, vignette, terse NEWS.md bullet, ASCII-only check, camelCase R,
  R CMD check clean. Merge origin/main before push.

## 7. Evaluation criteria (quality gates)

Correctness:
- Encoder fwd/bwd: match torch oracle and FD to ~1e-6 on golden fixtures.
- Decoder sensitivities: match closed-form + FD to ~1e-6.
- Every C++ `.Call` registered in init.c with correct arity (no runtime arity
  errors).

Parity to the paper (tolerances; the method is stochastic so exact match is not
expected -- paper's own methods agree to <0.2-0.5% on -2LL):
- Fixed effects within ~2-3% of the paper VAE column.
- Omega SDs within ~0.02-0.05; additive error a within ~0.02.
- -2LL within ~0.5%.
- Covariate selection: exact set match on theophylline (WT->ka, WT->V);
  on neonates match all dark-green (all-method) cells in Fig 4, VAE column.
- Reproducible: deterministic given a seed.

Software quality:
- R CMD check clean (0 errors/warnings; notes justified).
- No Unicode in source/docs (CRAN ASCII policy).
- No new hard runtime dependency: torch and the MIQP solver are Suggests only;
  tests skip gracefully when absent.
- camelCase R, C++17, test thread policy respected (tests/testthat.R), Config
  parallel honored.
- Commit at every phase/milestone; merge origin/main before each push.

Performance (informational, not gating unless regressive):
- Theophylline single run in the same order of magnitude as the paper (~6s).
- Neonates completes in reasonable wall-time with threading.

## 8. Risks and mitigations

- Analytic LSTM backprop bugs -> torch oracle + FD gradient checks (Phase 2).
- focei sensitivity reuse may not expose the exact per-obs Jacobian -> thin
  adapter, validated against FD (Phase 3).
- Matching stochastic paper numbers exactly is impossible -> tolerance-based
  gates, seeded runs.
- R-level covsel per-iteration cost -> run selection periodically, not every
  inner step; profile in Phase 5.
- torch unavailable on some platforms -> Suggests + skip_if_not_installed.
