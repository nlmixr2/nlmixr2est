# VAE-NLME estimation method for nlmixr2est -- implementation plan

Reference: Rohleff et al., "Redefining Parameter Estimation and Covariate
Selection via Variational Autoencoders: One Run Is All You Need", CPT:PSP 2025,
DOI 10.1002/psp4.70129. Paper PDF: ~/src/private/vae.pdf. Python original:
~/src/vae_nlme. R prototype (torch, not integrated): ~/src/vae-nlme-nlmixr2-main.

Work happens in git worktree ~/src/nlmixr2est-vae on branch `vae`.

## 1. Goal and definition of done

Add a first-class `est="vae"` estimation method to nlmixr2est. A VAE (LSTM
encoder over each subject's series -> Gaussian posterior q(eta|y); structural
ODE decoder via rxode2) is trained on an ELBO / BICc-ELBO objective to estimate
population parameters AND select covariates in a single run.

Definition of done = full parity with both paper case studies:
- Case Study 1: theophylline single-dose and multiple-dose, Table 1.
- Case Study 2: neonatal weight progression, 2425 subjects, Table 2 / Fig 4.

Built in small compartmentalized increments with a commit at every checkpoint;
full parity is the acceptance bar, not the first commit.

## 2. Locked decisions (interview 2026-07-09)

Architecture / integration
- Deliverable: first-class `est="vae"` (nlmixr2Est dispatch, nlmixr2FitData).
- Encoder + gradients: native C++ LSTM with analytic backprop
  (RcppArmadillo/Eigen); no torch runtime dependency.
- Decoder / ODE: rxode2 solver + its analytic sensitivity machinery for
  d(pred)/d(eta) (reuse focei / foceiCovAnalytic infrastructure).
- Process split: heavy C++ kernel (LSTM fwd/bwd, ELBO+grad, Adam, M-step
  sufficient stats) called per inner iteration from a thin R EM loop
  (saem.cpp pattern).
- Encoder architecture: reverse-engineer compiled encoder.c as ground truth.
  Established: single-layer UNIDIRECTIONAL LSTM (prototype README's
  "bidirectional" is WRONG). FC head:
  Linear(h_last | covariates) -> [mu | log_sigma | L_offdiag]; z = mu + L eps.

Encoder inputs / priors
- Encoder input: FAITHFUL -- LSTM consumes (time, DV) sequence; subject
  covariates concatenated at the FC head; dose/regimen handled entirely by the
  rxode2 decoder (encoder is dose-blind). General dosing supported via rxode2.
- Encoder prior mean: from ini() thetas. Encoder prior SD (sigma0): a FIXED
  SMALL default (paper-style), exposed as a vaeControl() knob -- NOT derived
  from ini() omega.

Model scope
- Parameter transforms: nlmixr2 standard set (log default; logit/expit, probit
  for bounded; identity) as the invertible h.
- Error models: FULL nlmixr2 error-model parity -- additive first, then
  proportional/combined, transform-both-sides, and M2/M3/M4 censoring
  (reuse existing analytic-censoring work).
- No-IIV parameters supported: a structural parameter may have no random effect
  (not an encoder latent), and covariates may still be selected on its
  population value. Encoder emits etas only for IIV parameters.
- IOV / time-varying etas: OUT OF SCOPE (one posterior per subject); documented
  as future work.

Covariate selection
- Candidates: AUTO-DISCOVER every subject-level non-dose/non-time column as a
  candidate on every parameter (paper's automatic philosophy), with
  model/vaeControl() override to restrict or pin.
- Encoding: match the paper's parameterization via nlmixr2 conventions --
  auto center/scale continuous covariates; categoricals via model-order factor
  levels (reuse focei string/factor-ordering fix).
- Selection solver: R-level MIQP/selection (Suggests: highs/Rglpk/ompr) fed
  sufficient statistics by the C++ kernel; per-parameter subset enumeration
  fallback. BICc-ELBO (L0) objective.
- Toggle: `vaeControl(covariateSelection = TRUE)` default; FALSE fits the given
  fixed covariate structure only (fast pop-only mode; also isolates estimation
  bugs from selection bugs).
- Output: return the estimated covariate effect matrix (beta) AND an updated
  nlmixr2 model/ui with selected relationships baked in (inspectable/refittable).

Objective / uncertainty / diagnostics
- Objective: compute BOTH linearization (Lin) and importance-sampling (IS)
  -2LL; store both in the fit. DEFAULT active OFV (used with AIC/BIC) = IS.
  (nlmixr2 fits may hold multiple OFVs, one active at a time.)
- Standard errors: linearization (FOCE-style) marginal-LL Hessian at
  convergence -> $parFixed SEs.
- Post-fit: FULL parity -- addCwres/addNpde, VPC, GOF, $parFixed all work;
  benchmark against nlmixr2 SAEM/FOCEI on the same data as an evaluation gate.

Controls / performance / reproducibility
- vaeControl(): expose paper schedule knobs with paper defaults --
  iters_burn_in=100, kl_warmup=50, gamma_iter=250, iters=300, L_iter=5,
  h_dim=25, learning rate, seed, sigma0, covariateSelection.
- Performance: parallelize the C++ encoder batch and rxode2 decoder solves over
  subjects via rxode2's existing thread pool (no new knob; parallel S-matrix
  pattern). Timings informational, NOT a pass/fail gate.
- Reproducibility: deterministic given a seed (reference RNG stream is not
  matched exactly across C++ vs torch).
- torch: test-only oracle (Suggests) to gradient-check the C++ encoder; never a
  runtime dependency.

## 3. Reference facts established

- Encoder: single-layer unidirectional LSTM; FC head outputs mu, log_sigma,
  strictly-lower-triangular L entries; L diag = exp(log_sigma); reparam
  z = mu + L eps.
- Decoder: theophylline single-dose is a closed-form 1-cmt solution; neonates is
  an ODE solve. Under nlmixr2est BOTH go through rxode2.
- ELBO = p(x|z) [gaussian nll] + DKL, DKL = p(z) - q(z|x), z ~ N(z_pop=C beta,
  Omega). M-step: closed-form z_pop, omega, a from EMA sufficient stats s1..s4.
- Schedule: burn-in (encoder only, tiny KL) -> KL annealing (0.01->1) -> main EM
  with covariate selection -> EMA smoothing of z_pop.

## 4. Open technical questions for Phase 0/1 (not user decisions)

- Exact LSTM gate ordering / weight layout PyTorch uses (i,f,g,o) so the C++
  backward matches the oracle.
- Whether focei d(pred)/d(eta) is exposed per-observation in the form the ELBO
  gradient needs, or needs a thin adapter.
- Exact input standardization (time, DV, covariates) the reference encoder uses.
- EM boundary-crossing cost of the R-level covariate-selection solve; run it
  periodically, not every inner step.
- Mapping of no-IIV parameters and bounded transforms into the z_pop / C beta
  structure.

## 5. Process / module map

- C++ kernel `src/vae*.cpp`: LSTM fwd/bwd, FC head, L build, reparam, ELBO +
  gradient, Adam, M-step sufficient stats. init.c hand-updated per new .Call
  (+ compileAttributes).
- rxode2: compiled decoder solve + analytic per-obs eta-Jacobian.
- R `R/vae.R` (+ vaeControl.R, vaeRxUiGet*.R): dispatch, control, UI->rxode2
  wiring, EM orchestration, schedule, R-level covariate selection, fit assembly.

## 6. Phased plan (each phase = commit checkpoint)

Phase 0 -- Reference capture & spec. [DONE]
- Reverse-engineer encoder.c / pop_parameter.c; write plans/vae-encoder-spec.md
  (exact layer config, init, gate order, FC layout, standardization). DONE.
- Golden fixtures: fixed seed + tiny input -> encoder forward outputs, ELBO,
  gradients from the torch oracle, saved under tests/testthat/fixtures/vae/.
  DONE via tools/vaeGenEncoderFixtures.R -> encoder_golden.rds (torch 0.17.0
  oracle; forward + autograd grads for all encoder params, all finite).
- pop_parameter.c (M-step/covsel) detail deferred to a Phase 5 spec.

Phase 1 -- Dispatch skeleton (R). [DONE]
- R/vae.R: vaeControl() (paper defaults + shared plumbing), getValidNlmixrCtl.vae,
  nmObjHandleControlObject.vaeControl, nmObjGetControl.vae, rxUiDeparse.vaeControl,
  nlmixr2Est.vae. Validates UI (random-on-id, mu-ref) + control, then stops
  "not yet implemented". Verified: control defaults, unused-arg rejection, est
  in nlmixr2AllEst(), dispatch reaches stub on theo_sd. NAMESPACE/Rd via document().
- Remaining for later phases: vaeRxUiGet* handlers, rxPipeline data prep,
  theta/eta ordering, no-IIV handling (wired as the C++ kernel lands).

Phase 2 -- C++ LSTM encoder (fwd + analytic bwd). [DONE]
- src/vaeEncoder.cpp: vaeEncoderFwdBwd -- LSTM fwd (PyTorch i,f,g,o), FC head,
  L build, reparam; analytic backward = FC + reparam + LSTM BPTT. Backward takes
  the two ELBO upstream signals (dLoss/dz, direct dLoss/dlogSigma).
- Registration: compileAttributes; MANUAL init.c (prototype + CallEntries arity
  13); MANUAL src/Makevars.in SOURCES_CPP (explicit object list, NOT wildcard --
  src/Makevars is generated/gitignored).
- Verified (tests/testthat/test-vae-encoder.R, no torch at run time): forward vs
  torch oracle ~1e-7; analytic backward vs torch autograd ~1e-6; vs finite
  differences ~5e-7. 13/13 expectations pass.

Phase 3 -- Decoder + analytic d(pred)/d(eta) via rxode2. [DONE]
- R/vaeDecoder.R REUSES focei's analytic sensitivity machinery (per user):
  .vaeDecoderModel wraps .foceiAnalyticAugModelDirs(ui, dirs=etas) -> augmented
  rxode2 model emitting rx_predf_ (f), rx_f1_<eta> (df/deta), rx_rvarf_ (R),
  rx_rvar1_<eta> (dR/deta) as solve columns; .vaeDecoderSolveSubject reads them
  via .foceiAnalyticSolveFA. No ODE-sensitivity code re-implemented.
- .vaeDecoderPxz: ELBO data term p(x|z) + d(p_x_z)/d(eta) = sum(rf*a + rR*aR)
  reusing focei rho(f,R) partials; this IS the encoder's dLoss/dz upstream
  (z = z_pop + eta).
- Verified (test-vae-decoder.R, 12/12): f vs theophylline closed form ~5e-12;
  df/deta vs FD ~1e-9; R=add.err^2, dR/deta=0 additive; p(x|z) grad vs FD ~1e-6.

Phase 4 -- ELBO + Adam + closed-form M-step (MILESTONE A). [DONE]
- R/vaeFit.R: .vaeElboStep (encoder C++ fwd -> decoder p(x|z) -> KL -> encoder
  C++ bwd), .vaeAdam*, .vaeMStep (zPop=mean mu; omega_k=mean[(mu-zPop)^2+
  (LL')_kk]; a=sqrt(SSR/Nobs); EMA gamma), .vaeTrain (burn-in tiny-KL ->
  main KL-anneal+M-step -> smoothing), .vaeFitModel. nlmixr2Est.vae runs training
  (returnVae=TRUE returns raw fit; full FitData deferred to Phase 6).
- Verified: test-vae-elbo.R full-pipeline gradient vs FD ~1e-6; test-vae-train.R
  end-to-end. Short 30/80 run (58s): ka=1.59(1.63), ke=0.0863(0.0867),
  V=31.84(31.97), a=0.725(0.71) -- fixed effects + a inside the ~2-3% gate;
  omegas close (wKa 0.62 vs 0.53). Full 100/300 schedule tightens further.
- Orchestration is thin R; heavy work = C++ encoder BPTT + rxode2 decoder solve.

Phase 5 -- Covariate selection (MILESTONE B).
- Full beta matrix; auto-discovered candidates + override; paper-style covariate
  encoding; BICc-ELBO via R-level solver (enumeration fallback); periodic
  selection in the EM loop; updated-model output.
- MILESTONE B: theophylline reproduces WT on ka and V, not ke.

Phase 6 -- Fit object: SEs, dual OFV, diagnostics, benchmark.
- Linearization-Hessian covariance -> $parFixed SEs; both Lin and IS -2LL stored
  with IS active for AIC/BIC/BICc; EBEs from encoder posterior/MAP;
  addCwres/addNpde/VPC/GOF; benchmark vs SAEM/FOCEI on the same data.

Phase 7 -- Error-model parity.
- Proportional/combined, transform-both-sides, and M2/M3/M4 censoring in the
  decoder/ELBO gradient (reuse analytic-censoring work). Validate each.

Phase 8 -- Neonates parity (MILESTONE C = done).
- Zero-order production + first-order elimination model; 5 params, 5 covariates;
  multiple obs/subject; 2425 subjects with rxode2 threading. Also theophylline
  multiple-dose.
- MILESTONE C: reproduce Table 2 (-2LL ~146370 Lin; BICc ~146566) and the Fig 4
  VAE covariate set.

Phase 9 -- Docs, tests, CRAN hygiene, merge/push.
- Roxygen, vignette, terse NEWS.md bullet, ASCII-only, camelCase R, C++17,
  R CMD check clean. Merge origin/main before push.

## 7. Evaluation criteria (quality gates)

Correctness
- Encoder fwd/bwd match torch oracle + FD to ~1e-6 on golden fixtures.
- Decoder sensitivities match closed-form + FD to ~1e-6.
- Every C++ .Call registered in init.c with correct arity.

Parity to the paper (proposed band; method is stochastic)
- Fixed effects within ~2-3% of the paper VAE column.
- Omega SDs within ~0.02-0.05; additive error a within ~0.02.
- -2LL within ~0.5%.
- Covariate selection: exact set match on theophylline (WT->ka, WT->V); on
  neonates match all dark-green (all-method) cells in Fig 4, VAE column.
- Deterministic given a seed.

Cross-method
- VAE estimates and OFV benchmarked against nlmixr2 SAEM/FOCEI on the same data
  (sanity, not identical-number, agreement).

Software quality
- R CMD check clean; ASCII-only; camelCase R; C++17; test-thread policy
  respected.
- No new hard runtime dependency: torch and the MIQP solver are Suggests only;
  tests skip gracefully when absent.
- Commit at every phase/milestone; merge origin/main before each push.

Performance (informational, not gating)
- Theophylline single run in the paper's order of magnitude (~6s); neonates
  completes in reasonable wall-time with threading.

## 8. Risks and mitigations

- Analytic LSTM backprop bugs -> torch oracle + FD checks (Phase 2).
- focei sensitivity reuse may not expose the per-obs Jacobian directly -> thin
  adapter validated vs FD (Phase 3).
- Full error-model parity (esp. censoring/TBS) is large -> isolated Phase 7 with
  per-feature validation; additive-only through Milestones A/B.
- Matching stochastic paper numbers exactly is impossible -> tolerance gates,
  seeded runs.
- R-level covsel per-iteration cost -> periodic selection; profile in Phase 5.
- torch/solver unavailable on some platforms -> Suggests + skip_if_not_installed.
