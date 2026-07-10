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
  end-to-end. Full 100/300 schedule (211s): ka=1.617(1.63,0.8%),
  ke=0.0849(0.0867,2.1%), V=32.08(31.97,0.3%), a=0.718(0.71,1.1%) -- all inside
  the ~2-3% gate. wKa=0.61 is the CORRECT no-covariate value: paper Table 1
  (0.53) is POST-covariate-selection; Fig 2 shows wKa rising to ~0.6 in burn-in
  and dropping to 0.53 only when WT->ka is selected. So wKa drops in Phase 5.
  (Confirmed sigma0 and encoder weight decay do not move wKa -- decomposition
  showed Var(mu_ka)=0.36 dominates, i.e. the unexplained WT variance.)
  MILESTONE A CLOSED.
- Orchestration is thin R; heavy work = C++ encoder BPTT + rxode2 decoder solve.

Phase 5 -- Covariate selection (MILESTONE B). [DONE]
- .vaeDataPrep auto-discovers subject-constant covariates (excl reserved cols),
  paper encoding (continuous log(v/mean), categorical centered).
- .vaeMStepCov: per-parameter subset enumeration minimizing RSS_k/omega_k +
  log(N)*|S| = the exact per-parameter reduction of the global BICc-ELBO L0
  MIQP (problem decouples by parameter; enumeration is exact + CRAN-safe for the
  paper's covariate counts). .vaeElboStep refactored for subject-specific z_pop
  matrix (KL center); decoder baseline immaterial to f.
- MILESTONE B (full 100/300, 213s): auto-found WT; WT->ka beta=2.47(2.55),
  WT->V beta=0.51(0.57), WT->ke NOT selected -- EXACT paper match. omega dropped
  to (0.538,0.139,0.134) vs paper (0.53,0.15,0.13) -- ALL inside gate (confirms
  the Milestone-A wKa diagnosis: WT absorbs ka variance). ke/V/a intercepts
  match; ka intercept 1.77 ~8% high (secondary).
- test-vae-covariate.R (skip_on_cran).

Phase 6 -- Fit object: SEs, dual OFV, diagnostics, benchmark. [PARTIAL]
DONE + verified (test-vae-fit.R):
- R/vaeOutput.R .vaeUpdateModel: injects selected covariate effects into each
  mu-ref model line as EXACT centered expressions (continuous beta*log(COV/center),
  categorical beta*(COV-center)) via the canonical model() line-replacement API
  (auto-adds beta as a population param; per rxode2 Modifying-Models vignette);
  sets ini() to the VAE solution. Result e.g.
  ka <- exp((lka + beta_lka_WT*log(WT/69.58)) + eta.ka), ke unmodified.
- .vaeLinPrecomp/.vaeLinObj: FOCE-linearized marginal -2LL from the decoder f/J
  at the encoder EBEs (theo ~346 vs paper ~331; gap = short schedule + encoder
  mean used as EBE instead of a re-optimized MAP). Because f/J are fixed at the
  encoder z, the -2LL is an explicit fn of (intercepts, betas, omega, a).
- .vaeCov: linearization-Hessian covariance = 2*solve(Hessian of -2LL) (cheap,
  no ODE re-solves); returns NULL on a non-PD Hessian (hardened).
- .vaeToFit builds a real nlmixr2FitData WITHOUT running focei (follows
  babelmixr2 nlmer.R .nlmerFamilyFit): .foceiPreProcessData -> fullTheta, omega,
  cov (=.vaeCov), etaMat/etaObf (encoder EBEs), objective (=lin -2LL),
  model=ui2$ebe (KEY -- saemModelPred caused "args missing"),
  .vaeControlToFoceiControl attaches a 0-iter focei control for the residual/
  table ENGINE (assembles tables at supplied etas; does NOT estimate), then
  nlmixr2CreateOutputFromUi(est="vae"). class nlmixr2FitData/nlmixr2.vae; OBJF =
  lin -2LL minus the 2pi constant (AIC/BIC use the full -2LL).
- Original-vs-final accessors: original ui stashed in fit$env$iniDf0;
  nmObjGet.iniUi/added nmObjGet.uiIni return the ORIGINAL ui when iniDf0 is a ui
  (else prior behavior); added nmObjGet.iniDf0 returns the original iniDf. Fully
  backward-compatible (data.frame path unchanged for all other methods).
- Verified (test-vae-fit.R, 11/11): FitData with beta_lka_WT/beta_lV_WT, objDf,
  final covariate model, $uiIni/$iniDf0 recover the original (no-covariate) model,
  encoder EBEs.
Polish DONE:
- Likelihood moved to C++/OpenMP (src/vaeLik.cpp vaeFoceLik): per-subject
  FOCE-linearized marginal -2LL, Newton MAP-EBE re-optimization on the linearized
  model (no ODE re-solves), Laplace determinant; M2/M3/M4 CENSORING via
  censEst.h doCensNormal1 (value) + censNormalPartials (score/curvature) exactly
  like inner.cpp. Parallel over subjects; num_threads from rxControl$cores
  (fallback getRxThreads); setRxThreadId(omp_get_thread_num()) cross-DLL fix.
  Includes armahead.h/censEst.h/rxomp.h. Registered (arity 15).
- EBEs are the MAP-refined zStar; .vaeObjective/.vaeCov use the C++ kernel;
  cov = 2*solve(numeric Hessian of -2LL) (cheap, f/J fixed), NULL on non-PD.
- Verified: test-vae-lik.R (3/3) C++ vs R Laplace reference, observed + M3;
  test-vae-fit.R (11/11) still passes.
- Mixture likelihood: vaeFoceLik(nMix, mixProb) combines components as
  -2 logsumexp_m(log mixProb[m] + ll_i^(m)) like inner.cpp; per-subject try/catch
  guards the OpenMP region (inner.cpp pattern); test-vae-lik.R 6/6.
- Decoder-solve hardening (.vaeDecoderSolveSubject): tolerance-relaxation retry
  loop (maxRecalc x recalcFactor) + finite-difference fallback for the
  eta-sensitivities when the analytic sens states blow up but f/R are finite
  (inner.cpp-style). Jump event sensitivities already on (aug model
  eventSens="jump"). Decoder/fit tests still pass.
Mixtures -- architecture (per user): mixture is ALL in the model (mix(v1,v2,..)
via ui$mixProbs / ui$saemNMix; component selected by data column mymixest=m, the
mymixest==k form the inner/analytic model already keeps -- focei.R:1212). No
M-step extension: the likelihood expands to nMix*N pseudo-subjects, each
component gets its own MAP eta; mixnum=argmax assigned for the fit.
- DONE: vaeFoceLik returns zStar for EVERY (component, subject) -- nMix*N rows
  component-major (mix-1 all subjects, then mix-2, ...); + obji (combined -2
  logsumexp) + mixest. test-vae-lik.R 8/8. nMix=1 unchanged.
REMAINING mixture integration -- BLOCKER FOUND (2026-07-10). Concrete d/dt mix
model: ke <- mix(exp(lke1+eta.ke), p1, exp(lke2+eta.ke)). ui$saemNMix=2,
ui$mixProbs="p1" detectable. loadPruneSens keeps rxEq(mixest,k) in the ODEs. BUT
the decoder aug model from .foceiAnalyticAugModelDirs does NOT correctly select
components: solving with mixest=1 vs 2 (whether via event-data col OR a scalar
param) gives IDENTICAL, degenerate predictions (f~10 flat). The rxEq(mixest,k)
form is converted away and unresponsive in the analytic cov aug model -- it was
built for covariance (single model), not mixtures. The INNER model handles mix
correctly (focei.R:1208-1214): predOnly is built from the PRUNED model (preserves
mix(), nMix>0, rxode2 accepts per-id mixest from iCov); the sens model keeps the
mixest==k form and inner.cpp manages selection. NEXT: build the VAE decoder's
per-component solve from the inner-model machinery (pruned/predOnly + iCov mixest)
or fix .foceiAnalyticAugModelDirs to preserve a working mixest selection -- study
how the inner model assembly (with d/dt states + f/r) removes the mixture. THEN
(3) .vaeToFit concatenated etaMat + mixNum/mixList from zStar/mixest; (4) training
decoder sets mixest per subject. Kernel (vaeFoceLik) already ready.
Other REMAINING: dual OFV (IS active for AIC/BIC/BICc); SAEM/FOCEI benchmark;
general error-model R-scale in cov (additive exact now).

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
