# nlmixr2est 7.0.0

## New features

- FOCEi guards each `theta`'s scaling constant per transform, keeping the
  derivative-based `scaleC` where it is well-behaved and falling back only in that
  transform's singular / out-of-range region.  Each parameter keeps `1/|init|`
  (linear/additive), `1` (log-normal), or its transform-specific formula while the
  value stays inside a band tailored to that transform (the linear band is
  `foceiControl(scaleCband=)`, default `c(0.1, 10)`).  Outside the band it falls
  back to the parameter's native magnitude `|init|` (NONMEM7 Appendix K, eq 15.2);
  for a bounded transform (`logit`/`expit`/`probit`/`probitInv`), if `|init|` is
  also out of range it uses the geometric middle of the band.  This fixes the
  singular cases that froze or destabilized the fit -- `1/|init|` blowing up for a
  small covariate initial estimate (and the issue-641 large-additive case, whose
  special handling this subsumes), `log()` at init `1`, `logit` at the interval
  midpoint, `factorial`/`gamma` at a digamma zero -- while leaving the well-scaled
  common case, and its results, unchanged.

- Fixed FOCEi `scaleC` for a `gamma()`-transformed population parameter: rxode2
  reports it as `curEval="lgammafn"`, which the scaling setup did not recognize, so
  it silently received the linear `1/|init|` default instead of its `1/digamma`
  scaling.

- The FOCEi family nudges a structural population parameter (`theta`) initialized
  at exactly `0` off zero before estimation, controlled by
  `foceiControl(zeroTheta=)` (default `0.001`), since a zero initial estimate has
  no native scale to scale by.  `+zeroTheta` is used when within the parameter's
  bounds, otherwise `-zeroTheta`; if neither is within the bounds it errors.
  Fixed parameters (including those fixed at `0`) are left untouched.  Residual
  error parameters are also left untouched: they carry their own scaleC, so an
  error `sd` set to exactly `0` still disables that component and a combined
  error model reduces to the smaller model as before.

- `foceiControl()` gains `shi21hMax` and `shi21hMin` (defaults `2.0` and `1e-4`),
  the upper and lower bounds on the adaptive shi21 finite-difference step used for
  FOCEi gradients (both the inner eta and, when `shi21maxOuter != 0`, the outer
  theta/covariate finite differences).  A larger upper bound lets the gradient of
  a flat, small-magnitude parameter clear the ODE-solver noise floor.  The NLM
  family keeps its own fixed bounds.

- The `imp` / `impmap` / `qrpem` importance-sampling family is faster: the
  theta-score M-step, the Monte-Carlo covariance (`covMethod="imp"`, the default)
  and the per-subject proposal build in the E-step are now parallelized over
  subjects, using the `cores` set in the control's `rxControl` (defaulting to
  `rxode2::getRxThreads()`), joining the already-threaded E-step weight loop.  All
  are bit-identical to the single-threaded run at any thread count.  This also
  fixes a latent bug in the theta-sensitivity M-step where `d(V)/d(theta)` was read
  from an under-sized per-thread lhs buffer, so a residual-error variance that
  depends on a structural parameter now contributes the correct M-step gradient.

- `est="vae"` now runs multi-threaded.  The per-subject encoder forward pass and
  the exact branch-and-bound covariate M-step (previously serial, dominating the
  EM and covariate-selection phases) are parallelized over the `cores` set in
  `vaeControl(rxControl=rxode2::rxControl(cores=))` (defaulting to
  `rxode2::getRxThreads()`), joining the already-threaded decoder solve.  The
  encoder forward pass and the covariate branch-and-bound are bit-identical to the
  single-threaded run.  The encoder backward (gradient) pass is also parallelized
  by default (`vaeControl(parEncoderBackward=TRUE)`); its cross-subject sum cannot
  be reduced in parallel bit-identically, so it is deterministic for a fixed
  `cores` but differs slightly from the serial path.  A note is added to the fit's
  `$runInfo` when it is active.  For bit-identical, fully reproducible results set
  `options(nlmixr2.identical=TRUE)` (flips the default to serial) or
  `vaeControl(parEncoderBackward=FALSE)`.

- The SAEM Louis stochastic-approximation FIM (`covMethod="sa"`) and the
  importance-sampling Monte-Carlo observed information (`covMethod="imp"`) are no
  longer tied to `est="saem"`/`est="imp"`.  They can now be requested as the
  `covMethod` of any mixed-effects estimation method (computed post-fit at the
  converged estimates) and switched onto any completed fit with
  `setCov(fit, "sa")` / `setCov(fit, "imp")`.  (The population-only NLM family
  has no random effects, so `sa`/`imp` do not apply there.)

- Several estimation families changed their default `covMethod` now that any
  covariance can be applied to any mixed-effects method:
    - the FOCEI family (`focei`/`foce`/`laplace`/`agq`) now defaults to the
      `"r,s"` sandwich (was `"analytic"`);
    - `est="vae"` now defaults to `"r,s"` (was `"analytic"`);
    - `est="nlme"` now keeps nlme's own covariance (`"nlme"`) by default;
    - the nonparametric family (`npag`/`npb`) now defaults to the
      importance-sampling covariance (`"imp"`).
  `est="saem"` (`"sa"`), the importance-sampling family (`"imp"`), the NLM family
  (`"r"`/optimizer Hessian), `est="advi"` (`"advi"`) and `fo`/`foi` (no
  covariance) keep their previous defaults.

- `vaeControl(bnbStrategy=)` selects the frontier discipline for the exact
  branch-and-bound covariate selection in `est="vae"`: `"lifo"` (default, the
  existing last-in-first-out depth-first search), `"fifo"` (first-in-first-out)
  or `"lc"` (least cost / best-first).  The solver is exact, so the selected
  covariates are identical for every strategy; only the search order differs.

- `est="vae"` can now estimate structural population parameters that have no
  random effect (are not mu-referenced).  Previously such a `theta` was frozen at
  its `ini()` value because the VAE only estimates parameters in the latent space.
  `vaeControl(nonMuTheta=)` selects the treatment: `"regress"` (default, matching
  `saemControl(nonMuTheta=)`) injects no eta and estimates each such theta directly
  by a bounded `bobyqa` regression against the FOCEi inner likelihood every M-step
  (bounds from the `ini()` lower/upper, blended with the M-step gain), recovering a
  no-random-effect population parameter without adding a spurious random effect.
  The eta-injection alternatives estimate it as `theta + mean(eta)` (the temporary
  eta is dropped from the output model): `"eta"` estimates the injected omega and
  the typical value; `"fix"` holds both the injected omega AND the typical-value
  theta fixed at their `ini()` values (nothing about the parameter is estimated, so
  it does not appear in the iteration table); `"none"` keeps the old freeze
  behavior.  A `$runInfo` note lists which parameters were converted.

- The analytic observed-information covariance is now the preferred `covMethod`
  across the mixed-model estimation methods, falling back to each method's
  previous default when a model is out of analytic scope:
    - `est="saem"` keeps the stochastic-approximation FIM (`"sa"`) as
      the default `covMethod`, now followed by `"analytic"` and `"linFim"`.
      `covMethod="analytic"` computes the FOCEI analytic covariance at the
      converged SAEM estimates and falls back to the linearized FIM (`"linFim"`)
      with a message when out of scope or not positive definite; the `"linFim"`
      covariance stays selectable via `setCov(fit, "linFim")`.
    - `est="nlme"` gains a `covMethod` argument
      (`c("analytic", "r,s", "r", "s", "nlme", "")`, default `"analytic"`) that
      recomputes the covariance at the converged nlme estimates; `"nlme"` keeps
      nlme's own standard errors (also available via `setCov(fit, "nlme")`).
    - `est="npag"`/`"npb"` (and their `m`/`i` variants), which previously
      reported no covariance, now compute one post-fit at the converged
      estimates (default `"analytic"` with the finite-difference fallback
      chain).
    - `est="imp"`/`"impmap"`/`"qrpem"` gain a `covMethod` argument
      (`c("imp", "analytic", "r,s", "r", "s", "")`, default `"imp"`).  `"imp"`
      is the Monte-Carlo importance-sampling covariance that the old
      `impCov=TRUE` selected (the `impCov` argument is removed); the other
      tokens compute the post-fit FOCEI covariance.
    - `est="advi"` keeps its variational covariance (`"advi"`) as the default but
      now honors an explicit `covMethod` (e.g. `"analytic"`) without overwriting
      it with the variational covariance.
    - `setCov()`/`getVarCov()` accept `covMethod="analytic"` post-fit.

- A general FOCE-family per-subject log-likelihood can now be built from an
  `rxode2` UI model and used outside of a fit, for MCMC/SAMBA-style algorithms
  (issue #414).  `foceiLikLoad()` compiles the inner model and sets up the
  problem (including the data) in memory, `foceiLikRun()` evaluates the
  individual log-likelihood at a supplied population parameter vector and eta
  matrix -- in parallel per subject -- and `foceiLikUnload()` frees it.  The
  likelihood type may be `"focei"` (with interaction), `"focep"` (FOCE+) or
  `"foce"` (NONMEM-style), and `foceiLikRun(type=)` selects the individual joint
  density `log p(y_i, eta_i)` (`"joint"`, the default) or the conditional data
  log-likelihood `log p(y_i | eta_i)` alone (`"cond"`).  Only one likelihood
  system may be loaded at a time; loading a second errors until the first is
  unloaded.
- `fit$etaCI` returns per-subject confidence intervals for each individual's
  eta, complementing the existing `fit$etaSE` and `fit$etaRSE`.  The intervals
  are the empirical-Bayes estimate plus/minus a normal quantile times the eta
  standard error, using the fit's `ci` level (default 0.95).  Like `etaSE`, it
  requires `CWRES` in the fit (add with `addCwres()` for non-focei methods).

- `est="agq"` now supports the analytic outer gradient (`agqControl(fast=TRUE)`),
  which was previously available only to the FOCEi family.  The AGQ objective is
  the FOCEi objective with one term swapped -- `l(etahat)` becomes
  `log(sum_k a_k)` over the quadrature nodes, while the `log det`, Omega and tbs
  terms are unchanged -- so its gradient reuses the same sensitivity solve and
  adds the node terms plus `tr(Ht^-1 dHt/dp)` for the node placement.  As with
  FOCEi this replaces the finite-difference outer gradient, so it is exact rather
  than a difference approximation and costs one augmented solve instead of one
  extra solve per parameter.  The quadrature nodes solve a cheaper 1st-order
  model than the eta-hat point needs (they never read the 2nd-order block), which
  is where most of the node cost goes once the grid grows.  Requires
  `interaction=TRUE`; a fit that cannot use it falls back to finite differences
  rather than failing.

- `covType="analytic"` now covers `est="agq"` as well (it previously declined for
  `nAGQ > 1` and fell back to the finite-difference covariance).  The AGQ
  observed information is the FOCEi one with the same single term swapped, so the
  `log det` half is reused unchanged and only the data half becomes an expectation
  over the quadrature nodes plus a covariance between their score contributions.
  At `nAGQ=1` it reduces to the FOCEi observed information exactly, and the FOCEi
  and Laplace results are unchanged.  Validated against a finite-difference
  oracle (tight ODE tolerance, Richardson extrapolation): the AGQ standard errors
  agree to that oracle's own noise floor.  As with the gradient, a model outside
  its scope -- a general or multi-endpoint residual variance, censoring, IOV, a
  finite `agqLow`/`agqHi` clamp, `cholSECov=TRUE`, or `interaction=FALSE` --
  reports why and keeps the finite-difference covariance.

- Requesting an unsupported `est=` method (e.g. a typo) now prints the available
  estimation methods grouped by category (Linearized, Integral approximation,
  Stochastic EM, Nonparametric, Machine learning, Optimizer (NLM family)) with a short
  description of each, instead of a single flat list.  Calling `nlmixr2()` with no
  arguments prints the same grouped list (and invisibly returns it).  The new
  `nlmixr2AllEstType()`
  returns the same information as a data frame, and each built-in method carries `type`
  and `description` attributes (e.g. `attr(nlmixr2Est.focei, "type")`) that third-party
  methods can set to join the list.
- `est="npag"`/`est="npb"` now PIN the current ODE solve during the
  residual-error (`err`) parameter optimization instead of re-integrating.  Those
  parameters do not change the prediction `f`, so each subject's states are
  cached at its posterior etas and the ODE is frozen (`op_focei.freezeOde`) while
  only `r` is recomputed -- for a mixture the frozen recompute reuses each
  component's cached states rather than re-solving them.  A structural regressor
  (which does move `f`, including an estimated per-component clearance) still
  re-solves.  Results are unchanged.
- SAEM mixture models now fix per-subject membership by default
  (`saemControl(mixProbMethod="regress")`, the new default): each subject is
  hard-classified to its best component once, held fixed, and the soft-EM
  responsibility step is skipped (reusing the existing responsibility-weighted
  machinery via a 0/1 `mixWeights`).  This avoids the soft-EM collapse (a
  component running away to a degenerate value) and is lower-bias on both
  well-separated and overlapping component evaluations; on heavily overlapping
  components it can be higher-variance, so the previous soft-EM behavior remains
  available with `mixProbMethod="regularized"`.  Because membership is fixed, the
  S-step solves each subject once under its own component (a per-subject mixest
  regressor) instead of running one MCMC chain per component -- roughly an
  `nMix`-fold reduction in ODE solves per iteration.  Split-ETA mixtures (a
  separate eta per component, which start symmetric and must differentiate during
  the fit) automatically fall back to soft-EM (`regularized`).

- SAEM warm-starts its residual-error parameters from the observed per-endpoint
  moments at the initial predictions (additive SD from `sqrt(mean(err^2))`,
  proportional SD from `sqrt(mean((err/f)^2))`), the same moment estimate
  `est="npag"`/`est="npb"` use -- `saemControl(residWarmStart=TRUE)`, the default.
  Because SAEM forms this at the unconverged population prediction, the
  proportional moment excludes near-zero predictions (where between-subject
  variability dominates) and the warm-started value is clamped to a sane multiple
  of the `ini` value.  Set `residWarmStart=FALSE` to start from the `ini`
  residual values.  For mixture models the warm-start is disabled (the poor
  population initial fit would inflate the residual and stop the components from
  separating).

- The proportional residual moment used to warm-start `est="npag"`/`est="npb"`
  (and now SAEM) guards against a near-zero prediction: the ratio is
  `err / (abs(f) <= 1e-6 ? 1 : f)`, so an `f` at (or near) zero no longer blows
  up the proportional moment.

- SAEM now estimates population `theta` parameters that have no associated random
  effect (the SAEM `phi0` fixed effects) by a bounded direct optimization of the
  observation likelihood each iteration -- `saemControl(nonMuTheta="regress")`, now the
  DEFAULT -- keeping them as plain directly-estimated regressors instead of stochastic
  `phi0` draws with a shrinking variance.  The optimization uses robust coordinate
  descent within a local trust region, honoring each theta's `ini`-block bounds, and
  holds `phi0` fixed once the optimizer owns it.  On a simulated model with three no-eta
  thetas (`ka`, `V`, a Hill power) this recovered them far more accurately than the old
  handling (e.g. the absorption theta RMSE dropped ~16x), at some extra runtime (the
  objective re-solves the ODE).  The previous behavior is available with
  `saemControl(nonMuTheta="eta")`.  For mixture models this falls back to the
  stochastic `phi0` block (the direct optimizer cannot partition a per-component
  structural theta by subject membership).

- `est="npag"`/`est="npb"` now ESTIMATE a mixture (`mix()`) model's component structural
  parameters (e.g. a per-subpopulation clearance) instead of holding them at their
  initial values.  The residual/regressor step optimizes them against the exact mixture
  negative log-likelihood `-sum_i log(sum_m a_m exp(cll_m))` (NONMEM7 eq 1.182),
  marginalizing over the components with the current proportions `a_m` (which the
  proportion update step moves); each per-component conditional log-likelihood carries
  the `-0.5*log(r)` penalty, so the additive residual does not collapse.  Verified: a
  two-subpopulation clearance model recovers both component clearances and the mixing
  proportion, with a non-zero additive SD.

- The per-endpoint residual moment warm start now attributes each observation to its
  endpoint via a new rxode2 accessor (`getIndCmt`, reading the CMT time-varying
  covariate), so a multi-endpoint model warm-starts each endpoint's residual from its
  own moment.  Requires the matching rxode2 (function-pointer table index 82).

- `est="npag"`/`est="npb"` now estimate the residual-error parameters with EXTENDED
  LEAST SQUARES at the individual predictions instead of the marginal likelihood.  The
  marginal likelihood over a flexible nonparametric support rewards a vanishing residual
  (each support point can then fit its subjects arbitrarily well), so the residual could
  drift toward zero.  The residual step now minimizes the exact conditional normal
  negative log-likelihood `sum_obs(0.5*(f-dv)^2/r + 0.5*log(r) + 0.5*log(2*pi))` at the
  posterior-mean etas (equivalently extended least squares -- same minimizer) -- the
  `0.5*log(r)` term penalizes `r -> 0`, giving the saem/focei
  residual (e.g. theophylline add.sd ~ 0.73, prop.sd ~ 0.15) rather than a collapsed one.
  Each variance-scale parameter is warm-started (and, for a single scale per endpoint,
  set) from the saem-style per-endpoint moment: an additive SD from `sqrt(mean(err^2))`, a
  proportional SD from `sqrt(mean((err/f)^2))`, both on the transform-both-sides scale (so
  lognormal / box-cox are handled on the transformed residual).  A non-mu structural
  "regressor" is optimized in the same step, with the posterior-mean etas re-derived per
  candidate so the eta grid cannot stale-absorb the structural shift (this identifies it,
  e.g. recovering theophylline's clearance from a deliberately-wrong start).  After the
  residual + regressor thetas converge, a final adaptive-grid pass re-optimizes the
  support with those thetas held CONSTANT, so the support is the nonparametric MLE of the
  mixing distribution for the fitted residual and the D(F) global-optimality certificate
  is restored (~0).  npb runs the same residual/regressor step inside its sampler.  (A
  mix() model's structural component parameters are held at their initial values -- the
  ELS step is not mixture-aware; the components are handled by the mixture marginalization
  and proportion update.)

- `est="npag"` now picks the initial grid size automatically from the model's
  dimensionality when `npagControl(points=)` is not supplied: `max(2028, 512 * n_eta)`
  (2028 is the Pmetrics NPAG default, which covers a low-dimensional model but grows
  sparse and can collapse in high dimensions).  Theophylline (3 etas) resolves to
  2028 (matching Pmetrics); warfarin (8 etas) to 4096.  Supply `points` to override.

- `est="npag"` is more robust on high-dimensional models (many etas), validated by a
  golden comparison against Pmetrics NPAG on the Warfarin PK/PD model (transit
  absorption + Emax turnover, 8 parameters): the per-cycle Psi build is now per-row
  log-sum-exp normalized on the non-gamma path too, so a hard subject's conditional
  density cannot underflow a whole row to zero (which aborted condensation); the
  Burke interior-point solve ridges the Newton matrix and retries instead of aborting
  when it is ill-conditioned; and `npagControl()` exposes `gridWidth` and
  `gridBounds` (`"auto"`/`"ini"`/`"both"`) so a bounded, high-dimensional grid can be
  focused on the plausible region (an unbounded box collapses the support).  These
  are numerically transparent for well-conditioned fits (the normalization restores
  the exact objective; Burke weights are scale-invariant).

- `est="npag"`/`est="npb"` now estimate non-mu structural fixed-effect parameters
  (a theta with no eta, e.g. `ke <- exp(tke)`, which npag's grid otherwise does not
  cover -- it covers only mu-referenced and residual/likelihood parameters).  By
  default they are optimized as "regressors" in the residual step: the bounded
  `bobyqa` moves them alongside the residual parameters, re-solving the ODE per
  candidate (they feed the states, so the ODE freeze is turned off for that step).
  This identifies them sharply -- e.g. recovering theophylline's clearance from a
  deliberately-wrong start, and a bimodal mixture proportion (p1 = 0.70) that the
  grid alternative recovered only weakly.  The opt-in `npagControl(muExpand=TRUE)`
  instead uses the saem-style mu-expansion: inject a pseudo-eta
  (`ke <- exp(tke + eta.tke)`), grid-estimate, and recover it as a fixed effect at
  finalization (support-mean folded into the theta, injected random effect collapsed;
  the injected eta carries a FIXED omega, excluded from the free omega objective like
  IOV, so it also works in mixture models).  `residOptimize="none"` holds the
  structural regressors together with the residual parameters.  (A non-mu-referenced
  ETA -- an eta with no paired theta -- needs neither: the npag box already covers
  every eta, so it is a grid dimension estimated as a pure random effect.)

- `est="npag"` now supports generalized (non-normal) / user-`ll()` likelihoods.
  The nonparametric objective sums the inner per-observation llikObs, which for a
  non-normal endpoint is exactly the user's log-likelihood, so the objective is
  already correct; the residual/likelihood parameters (e.g. a Student-t's degrees
  of freedom, `iniDf$err` non-NA) are estimated with the same frozen-ODE bounded
  step as the residual parameters.  Freezing the ODE during that step is valid only
  when every optimized parameter feeds the post-solve f/r alone (err-tagged) -- if a
  non-err parameter ever enters the optimized set the step re-solves instead.  gamma
  is forced off (a non-normal endpoint has r == 1).  A non-mu-referenced structural
  fixed-effect parameter cannot be placed on the grid and is held at its initial
  value, reported in the fit's `$runInfo`.  `est="npb"` handles non-normal endpoints
  too (the Gibbs sweep sums the same llikObs).

- `est="npb"` now runs the residual/regressor optimization (previously it held the
  residual-error and non-mu structural "regressor" thetas at their initial values and
  only sampled the mixing distribution).  With the sampled mixing distribution held
  fixed, the same bounded `bobyqa` step npag uses fits the residual thetas (add/prop/
  lnorm/lambda/ar) and any structural regressor -- recovering, e.g., theophylline's
  clearance from a deliberately-wrong start.  `npbControl(residOptimize=)` selects it:
  `"alternate"` (default) re-fits during burn-in and then holds the thetas fixed for
  the sampling phase (so every collected draw shares the converged residual scale),
  `"final"` fits once at the converged draw, `"none"` holds them at their initial
  values.  Unlike npag, npb does not optimize the assay-error multiplier (gamma) -- the
  residual thetas are fit directly.

- `est="npag"` and `est="npb"` now support mixture (sub-population) `mix()` models.
  Each subject is split into per-component pseudo-subjects and the conditional
  likelihood is marginalized over the components using the mixture proportions
  (`p(y_i | phi) = sum_m mixProb_m * p(y_i | phi, component m)`).  npag updates an
  estimated proportion each cycle by an EM step (support points and weights held
  fixed); npb samples the proportions inside the blocked Gibbs sweep -- each subject
  draws a component from its posterior responsibility and the proportions are drawn
  from Dirichlet(1 + component counts), with the posterior-mean proportions reported
  in `$env$npbMixProb`.  A `fix()`ed proportion is held at its ini value in both.

- `est="npb"` now supports multiple independent chains (`npbControl(nchains=)`):
  the stick-breaking Gibbs sampler runs once per chain (seed offset per chain),
  the posterior-mean draws are pooled, and a Gelman-Rubin R-hat per eta is
  reported in `$env$npbRhat` (~1 at convergence; > ~1.1 flags non-convergence).

- `est="npb"` is faster: the two per-sweep loops that re-solve the ODE serially
  (the support-location Metropolis-Hastings step, and the mixture-proportion
  Gibbs step for `mix()` models) now solve their per-subject conditional
  likelihoods in parallel over subjects, matching the already-parallel Psi build.
  The proposal and accept/reject draws stay serial in their original order, so a
  fixed-seed fit is bit-for-bit identical regardless of thread count.

- `est="npag"` is faster: it no longer does a redundant full conditional-density
  build at the first cycle (the degeneracy check now reads the working build's
  per-subject maxima), and the one-time D(F) global-optimality scan is smaller by
  default and configurable via `npagControl(dfScan=)` (`-1` auto, `0` to skip the
  certificate, or an explicit scan size).  Neither change affects the fitted
  support, Omega, thetas, or objective.

- `npagControl(cores=)` and `npbControl(cores=)` set the number of threads used
  for the parallel per-subject conditional-likelihood solves.  The default
  (`NULL`) uses the current `rxode2` thread count (`rxode2::getRxThreads()`); an
  integer sets the thread count for the fit and restores it afterwards.

- `est="saem"` now fits general log-likelihood (`ll() ~ expr`) models the saemix
  way (the model returns the per-observation loglik; the standard MCMC kernels use
  `-ll` as the observation loss).  The solve event data keeps `DV` when the model
  references it (previously dropped, so the likelihood solve errored
  "parameter(s) required for solving: DV"); the fixed-effect-only (phi0)
  parameters are optimized with the bounded `bobyqa` honoring the ini-block bounds
  (so a likelihood SD stays non-negative).  Normal-endpoint saem is unchanged.

- Nonparametric engines (cont.): `est="npag"` optimizes the residual parameters
  with the bounded `minqa::bobyqa`, honoring the ini-block lower/upper bounds of
  each residual parameter (e.g. an additive SD stays >= 0, an AR correlation in
  (-1,1)).  An unbounded optimizer could wander into an invalid region, so newuoa
  / nelder-mead are no longer used for the residual step (the `residType` control
  is removed).

- SAEM general log-likelihood: the fixed-effect-only (phi0) refinement step
  (saemix "ind.fix10", `distribution=general`) is now optimized with the
  same derivative-free optimizers as the residual step (nelder-mead / newuoa,
  selected by `type`) instead of L-BFGS-B -- the model emits no analytic
  d(ll)/d(phi0), so the previous finite-difference-gradient L-BFGS was pure
  overhead.  phi0 does not enter the ODE, so the states are solved once and held
  fixed while phi0 is optimized (ODE-freeze), each evaluation recomputing only the
  log-likelihood.  The SAEM-side L-BFGS plumbing (phi0 gradient, trampolines,
  `lbfgs*` config) is removed; FOCEI's `outerOpt="lbfgsb"` is unaffected.

- Nonparametric engines (cont.): the `npag` residual-parameter optimization now
  freezes the ODE states -- the inner likelihood solves each (support point,
  subject) once and re-evaluates only the output `f`/`r` for each candidate
  residual theta, skipping the (costly) re-integration.  Results are identical to
  the full re-solve; on a combined-error theo fit it is ~35% faster, and much more
  for models with expensive ODEs.  Exposed as a general `freezeOde` option on the
  inner likelihood (off by default, so all other engines are bit-identical).

- Nonparametric engines (cont.): `est="npag"` now estimates the residual-error
  parameters generally.  A single variance-scale parameter (pure additive or
  proportional) is handled by the fast gamma up/down search folded into that
  theta; anything else -- combined additive+proportional (the add/prop ratio a
  single gamma cannot recover), multiple endpoints (each `add.sd`/`prop.sd`), and
  transform (`boxCox`/`yeoJohnson` lambda) or autocorrelation (`ar`) parameters --
  is optimized against the nonparametric objective with the support points and
  weights held fixed, using the same optimizers as SAEM (`residType`: `"newuoa"`
  default, or `"nelder-mead"`), with gamma as a warm start.  The `residOptimize`
  control selects `"alternate"` (default, every cycle), `"final"` (once at the
  converged support), or `"none"` (hold at ini).  On a simulated two-endpoint
  model npag recovers `add.sd1`=0.20 and `add.sd2`=1.47 (truth 0.20 / 1.50),
  matching FOCEI, where a single global gamma had forced them equal; on simulated
  AR(1) data (true `ar1.cor`=0.6) it recovers ~0.54 from a 0 start where gradient
  FOCEI stalls at the `ar1.cor`=0 saddle.  The reported residual reflects the
  estimate.  Note: because the support distribution is flexible it can absorb
  additive residual scatter, so the additive term of a combined error model may be
  smaller than a parametric fit (documented in `?npagControl`).

- Nonparametric engines (cont.): the `npag`/`npb` conditional likelihood now
  folds in the transform-both-sides (dTBS) per-observation Jacobian, so `lnorm`,
  `boxCox`, and `yeoJohnson` residual models are handled correctly and lambda-type
  transform parameters are estimable.  Proportional and combined additive +
  proportional error are supported, and the global-optimality certificate D(F) is
  now evaluated at the fitted gamma (so it reaches ~0 for proportional/combined
  models with gamma optimization on).  A model whose transform link sees a
  non-positive prediction (e.g. `lnorm` at an observation where the prediction is
  0) now raises a clear error instead of an Armadillo empty-matrix crash.

- Nonparametric engines (cont.): the `npag`/`npb` engines now support fixed
  parameters.  Fixed population `theta`s (including fixed residual parameters such
  as `add.sd <- fix(0.7)`) are held at their ini value.  Fixed-`Omega` etas -- for
  example a fixed inter-occasion variance `iov.ka ~ fix(0.05) | occ` -- remain
  support-point dimensions but keep their variance held at the fixed value instead
  of being estimated, so IOV models fit.

- Nonparametric engines (cont.): `est="npag"` now reports the global-optimality
  certificate D(F) (`$env$npagDF`; ~0 certifies the nonparametric maximum
  likelihood), records a per-cycle parameter-history trace through the shared
  scale.h printer (`$parHistData`), and installs the reported `Omega` masked by
  the model's sparsity so correlated-eta models keep their off-diagonal terms.
  AR(1) and other transform-both-sides / structured residual models are supported
  (the residual enters as `f + sqrt(r)*eps`, so any structure carried in `r` flows
  through the conditional likelihood).

- Validation: a bimodal-recovery test confirms `est="npag"` recovers a
  two-subpopulation (fast/slow absorption) parameter distribution -- both modes
  carry substantial weight and the recovered cluster means land near the
  simulated truth -- the defining nonparametric capability a single-mode
  parametric random-effect model cannot reproduce.

- `est="npb"` (nonparametric Bayes) is now a usable engine: a truncated
  stick-breaking Dirichlet-process mixture sampled by a blocked
  Metropolis-within-Gibbs sampler (cluster assignments, stick weights, MH support
  locations).  It reuses the same conditional-likelihood primitive as npag and
  returns a `nlmixr2FitData` with the posterior mixing distribution
  (`$env$npbSupport`/`npbWeights`), per-subject posterior-mean etas, and posterior
  draws of the population mean (`npbMeanDraws`) for Bayesian credible intervals.
  `npbControl()` exposes `points` (truncation K), `alpha`, `burnin`, `nsamp`,
  `propSd`, and `seed`.  (Gelman-Rubin multi-chain convergence is a follow-up.)

- `est="npag"` is now a usable engine: it returns a standard `nlmixr2FitData`
  object with the nonparametric population summary (mean + variance mapped to the
  reported `theta`/`Omega`), per-subject posterior-mean etas, and the discrete
  support-point distribution attached to the fit (`$env$npagSupport`,
  `npagWeights`, `npagPosteriorEta`, `npagGamma`, `npagNspp`).  `npagControl()`
  exposes `points`, `cycles`, and `gammaOptimize`.  (The reported `Omega` uses the
  support-point variances; correlated-Omega models and the global-optimality
  certificate are follow-ups.)

- Nonparametric engines (cont.): added the residual-error magnitude (gamma)
  optimization inside the NPAG cycle (per-cycle up/down search).  Gamma scales the
  residual variance inside the FOCEi inner likelihood, so censoring (BLQ/ALQ via
  the M3 censored likelihood -- the normal tail probability below/above the limit)
  and transform-both-sides are handled correctly at the scaled error.  The objective uses a log-sum-exp row normalization for
  numerical stability.  Generalized (non-normal) likelihoods are not supported and
  are rejected with an error.  Note: the npag/npb objective is the nonparametric
  marginal log-likelihood and is NOT comparable to NONMEM/FOCEI -2LL.

- Nonparametric engines (cont.): assembled the NPAG adaptive-grid cycle (Yamada
  Alg 1) -- Sobol grid, Psi, Burke IPM, weight/QR condensation, adaptive-grid
  expansion (`npExpandGrid`), and the eps/F convergence controller.  Runs
  end-to-end on Theophylline (exposed as `npagCycle_` ahead of the full
  fit-object wiring).

- Nonparametric engines (cont.): added the Sobol initial grid (`npSobolGrid`),
  weight-threshold and QR rank-revealing condensation (`npCondenseWeights` /
  `npCondenseQR`), and the eta-space support-point box (`.npEtaBox`,
  control-selectable via `gridBounds`/`gridWidth`).

- Nonparametric engines (cont.): added the conditional-likelihood primitive
  (`npEvalCondLik`) and the parallel Psi-matrix builder (`npBuildPsi`), reusing
  the FOCEi inner solve so residual-error models, transform-both-sides and
  censoring carry over unchanged.

- Scaffolding for two native nonparametric estimation engines, `est="npag"`
  (nonparametric adaptive grid) and `est="npb"` (nonparametric Bayes), plus their
  mu-referenced sugar variants `mnpag`/`inpag` and `mnpb`/`inpb` (OLS and IRLS
  covariate M-step).  Both reuse the FOCEI inner likelihood machinery; the
  estimation loop runs in C++.  The algorithm itself is added in subsequent
  releases (the drivers currently report that estimation is not yet implemented).

- Fix the covariance matrix (`$cov`) of a bounded-parameter fit run with an
  unbounded method (e.g. `saem`): the internal `rxBoundedTr.<name>` name leaked
  into `$cov` and the back-transform Jacobian was not applied to it, so the
  reported standard errors were on the internal (transformed) scale.  `$cov`
  (and the stashed full theta+Omega covariance) are now renamed to the original
  parameter names and Jacobian-corrected; Omega and residual terms are
  untransformed so they pass through unchanged.
- The nlm parameter-history machinery can now be driven by an external
  optimizer.  `nlmerSolveGrad()` gains a `record` argument that logs the
  evaluation's population parameter estimate (the per-subject mean of the
  `phi` columns) into the resident scale, and `nlmGetParHist()` is now
  exported so an externally-optimized engine (e.g. `babelmixr2`'s nlmer,
  driven by `lme4::nlmer`) can recover the accumulated parameter history
  before `.nlmFreeEnv()`.  A new optional `showOfv` field in the nlm solve
  control hides the objective column for these engines (they record
  parameters only).

### New estimation methods

- `est = "advi"` (`adviControl()`): automatic differentiation variational
  inference (Kucukelbir et al. 2017), mean-field or block full-rank family,
  point-estimate or full-Bayes mode.  The variational gradient comes from the
  FOCEi forward sensitivities and the whole optimization runs in one C++ call,
  reproducibly and independent of the thread count.

- `est = "impmap"` and `est = "imp"` (`impmapControl()` / `impControl()`):
  importance-sampling EM in the style of NONMEM `METHOD=IMP`, with the E-step
  proposal at each subject's MAP mode (`impmap`) or running conditional mean
  (`imp`).  Supports mu-referenced, mixture, bounded and `fix()`ed models; the
  reported objective is a FOCEi evaluation at the EM estimate.  Quasi-random
  (Sobol) importance sampling (`qr=`, Leary & Dunlavey 2012) and SIR M-step
  acceleration (`sir=`) are available and stay thread-count independent.

- `est = "qrpem"` (`qrpemControl()`): sugar for the impmap EM with `qr = TRUE`
  and `sir = TRUE`.

- Mu-referenced FOCEI family: `mfocei`/`ifocei`, `mfoce`/`ifoce`,
  `mfocep`/`ifocep`, `magq`/`iagq`, `mlaplace`/`ilaplace` (with matching
  `*Control()` functions).  Mu-referenced population and covariate-coefficient
  thetas are profiled out of the outer optimizer by an in-C++ OLS (`m*`) or
  IRLS (`i*`) regression; bounded mu parameters are regression-updated with a
  clamped step.  New `foceiControl()` options `muModel`, `muRefCovAlg`,
  `muModelTol`, `muModelMaxCycles`, `muModelClampRetries`.

- `focep`/`mfocep`/`ifocep`: the `foce`/`mfoce`/`ifoce` methods with
  `foce = "foce+"` forced.

- `*f` convenience methods (`focef`, `foceif`, `focepf` and the mu/irls
  variants): the base method with `foceiControl(fast = TRUE)` as the default.

### FOCEI / FOCE

- `foceiControl(fast = TRUE)`: analytic FOCEI/FOCE outer gradient from Almquist
  (2015) sensitivity equations, solved for all subjects in one threaded rxode2
  solve; out-of-scope models fall back to finite differences.  Covers censored
  M2/M3/M4, an estimated boxCox/yeoJohnson lambda, `matExp()`/`indLin()`, foce+,
  modeled dosing (`f()`/`lag()`/`rate()`/`dur()`), and mu-referenced covariate
  reuse.  Under `fast` the outer optimizer defaults to `lbfgsb3c` and `mceta`
  defaults to `-2` (Eq-48 warm-start of the next inner problem).

- `covMethod = "analytic"` (folding in the old `covType`): exact analytic
  observed-information covariance for FOCEI/FOCE matching NONMEM
  `$COV MATRIX=R`, covering additive/proportional/combined error, censored
  M2/M3/M4 (`censOption = "gauss"`), estimated lambda, foce+,
  `matExp()`/`indLin()`, and mu-referenced/covariate parameters; out-of-scope
  fits fall back to the finite-difference sandwich.

- `covFull = TRUE` (now the default) reports the full theta + residual + Omega
  covariance for both the analytic and finite-difference methods, with Omega
  rows named by the random effect (`om.eta.cl` / `cov.eta.cl.eta.v`).
  `covMethod = "r,s"` is a true sandwich `solve(Rfull) %*% Sfull %*% solve(Rfull)`,
  `"s"` is `solve(Sfull)`, `"r"` is `solve(Rfull)`; `covFull = FALSE` keeps the
  theta-only shape.

- `foceiControl(foce = c("nonmem", "foce+"))`: `"nonmem"` (default) freezes the
  FOCE residual variance at the `eta = 0` prediction to match NONMEM; `"foce+"`
  keeps the live conditional variance.

- `foceiControl(censOption = c("gauss", "laplace"))`: censored (M2/M3/M4/BLQ)
  inner-Hessian treatment; `"gauss"` (default) matches common tools, `"laplace"`
  uses the exact censored second derivative.

- `foceiControl(warm = c("calc", "save"))`: `"calc"` (default) warm-starts each
  `n1qn1` inner problem from the eta Hessian recalculated at the current theta.

- `sensMethod` (nlm-family controls and `foceiControl()`): forward or in-engine
  discrete adjoint (`"adjoint"`) ODE parameter sensitivities; `"default"` reads
  `getOption("nlmixr2est.adjoint")`.

- Residual (error-model) parameters are now included in the focei-family
  covariance (only fixed, IOV and mixture-probability thetas skip).

- Mixture (`mix()`) support for `focei`/`foce`/`fo`/`foi`.

### SAEM

- `saemControl(covFull = TRUE)` (default): full theta + residual + Omega
  covariance from the linearized FIM.  New `covMethod = "sa"`
  (stochastic-approximation Fisher information, Kuhn & Lavielle 2005).
  `parHistData` records off-diagonal Omega block covariances.

- `saem` fits general log-likelihood endpoints (`ll(name) ~ <expr>`,
  e.g. time-to-event); fixed-effect-only parameters are refined by bounded
  L-BFGS-B (`saemControl()` gains `lbfgsLmm`/`lbfgsFactr`/`lbfgsPgtol`/
  `lbfgsMaxIter`).

### matExp() / indLin()

- Matrix-exponential / inductive-linearization models estimate with the focei,
  nlm and SAEM families, matching the equivalent ODE model; compartments are
  ordered source-first from the `k_<from>_<to>` graph so default dosing is
  placed correctly.

### Output and utilities

- The nonparametric eta-space outputs now carry the eta names: for
  `est = "npag"` the support-point matrix (`fit$env$npagSupport`) and posterior
  eta matrix (`fit$env$npagPosteriorEta`) get eta column names; for
  `est = "npb"` the same two matrices plus the posterior mean draws
  (`fit$env$npbMeanDraws`) get eta column names, and the per-eta R-hat vector
  (`fit$env$npbRhat`) gets eta row names.

- `est = "npb"` now prints its per-sweep iteration history through the shared
  iteration printer (like every other method) and stores it on the fit as
  `parHistData`; the sampler's results are unchanged (bit-identical).

- The importance-sampling (`covMethod = "imp"`) covariance step now shows a
  progress bar over its finite-difference evaluations, like the focei
  covariance step (shown when iteration printing is on).

- New `vaeCovariates()` returns the covariates `est = "vae"` would explore.

- New `formatMinWidth()` for shorter `$parFixed` display; `$parFixed` is rebuilt
  with data.frame operations (#346, #516).

- All estimators share one iteration printer (`iterPrintControl()`,
  `src/scale.h`) with a common row layout; analytic gradients are tracked as
  their own `parHist` type and the fit header reports the gradient and mu-model
  used, e.g. `(outer: lbfgsb3c; grad: analytic; mu: irls)`.

- `est = "vae"` training runs entirely in C++ (`vaeTrainCpp_`) and
  reparameterizes the inner problem in place, so fits are substantially faster.

### Changed defaults

- `outerOpt = "nlminb"` for the finite-difference methods (`"lbfgsb3c"` when
  `fast = TRUE`), `sigdig = 4`, `mceta = -2`, `censOption = "gauss"`,
  `foce = "nonmem"`, `covMethod = "analytic"`.

## Bug fixes

### Estimation

- Fixed the FOCEi `scaleC` band guard corrupting `est="vae"` covariate selection.
  The guard only rescues a genuinely-computed derivative-based scaling constant
  (`> 0`) now; an uninitialized `scaleC` of exactly `0` is left for the usual
  min/max clamp instead of being overwritten with `|init|`.  The overwrite had
  broken VAE covariate discovery on theophylline (no covariates selected, betas
  collapsed to `0`).

- `est="vae"` with `covariateSelection=FALSE` now estimates the covariate
  coefficients written into the model -- both linear (`beta*WT`) and transformed
  (`beta*log(WT/70)`) effects -- rather than holding them at their `ini()` value.
  They are fit in place by the regress M-step regardless of `nonMuTheta`
  (previously fixed under `nonMuTheta="none"` and errored under `"fix"`/`"eta"`);
  a coefficient set with `ini(... ~ fix())` still stays fixed.

- `est="impmap"` now estimates the non-mu structural and residual-error thetas of a
  general (custom `ll()`) likelihood model.  For such an endpoint `rx_pred_` is the
  log-likelihood itself and `rx_r_` is `0`, so the Gauss-Newton M-step skipped every
  observation (`V<=0`) and left those thetas frozen at their initial values; the
  M-step now uses the analytic `d(ll)/d(theta)` directly (empirical-Fisher
  information), so a raw `ll()` fit recovers the same parameters as the equivalent
  `add()` model.

- `est="npag"`/`est="npb"` no longer error with `unused argument: 'dfScan'` when the
  post-fit importance-sampling covariance is recomputed (the `dfScan` field leaked
  into the down-converted `foceiControl`).

- `est="npag"`/`est="npb"` with a transform-both-sides (`lnorm`/log/box-Cox) endpoint
  whose model prediction is non-positive at some observation (e.g. a pre-dose
  observation where the structural prediction is `0`) now records a note in the fit's
  `$runInfo` instead of silently fitting the rxode2-floored value with no indication.

- `est="vae"` with `nonMuTheta="regress"` now shows the regressed non-mu-referenced
  thetas in the iteration table and parameter history.  The M-step `bobyqa`
  regression already estimated them, but they were omitted from the printed
  parameter walk (only the latent-space thetas, omega, and residual error were
  shown), so their progress was invisible; they are now appended to each row with
  the correct back-transform.

- `est="vae"` covariate selection no longer silently selects nothing at 32
  candidate covariates.  The best-subset step enumerated all `2^nCov` subsets,
  which is undefined behavior at `nCov = 32` (`1u << 32` wraps to `1`, so only
  the empty model was ever tried) and intractable well before that.  It now uses
  an exact branch-and-bound over the same L0/BIC objective, returning the
  identical optimum while scaling to a few dozen covariates.  The selection
  penalty now also follows the reference implementation's warmup ramp, tunable
  via `vaeControl(covSelectAlpha=)` (default `2`, ramped to `1` over `klWarmup`
  iterations); ramp iterations are labeled `CovSel ramp` in the iteration table.

- `est="vae"` no longer errors with `replacement has 0 rows` on data that has no
  `AMT` column (dose-free datasets such as the neonate weight data); such rows
  are now treated as observations (`EVID = 0`).

- `est="saem"` no longer dies with `argument is of length zero` when building the
  SAEM model list.  Some `rxode2` versions omit the `ar` column from a model's
  `predDf`, and the SAEM autocorrelation helpers indexed that column directly; they
  now fall back to the `iniDf` (`err == "ar"`) representation when the column is
  absent.

- A mu-referenced or method-variant FOCEi fit (`ifocei`, `mfocei`, `foce`,
  `focep`, `agq`, `laplace`, and the `*f` fast variants such as `ifoceif`) that
  needed to restart -- for example after a zero/bad-gradient theta reset -- died
  with `focei$control must be a focei control object`.  These controls are all
  built by `foceiControl()` and then reclassed to their own class, so they do
  not carry `"foceiControl"` in their class vector, and the restart-path
  environment check rejected them even though the fit had been set up from a
  valid control.  The check now recognises the whole FOCEi control family.
- Models that combine `linCmt()` with ODEs (for example a solved PK driving an
  effect-compartment ODE) now estimate correctly with the FOCEi and nlm
  families; the linear compartments are solved as ODEs for those methods.
  Previously the sensitivity compartments those methods add (one per eta for
  FOCEi, one per theta for nlm) shifted `depot`/`central` past the compartment
  numbers the data was translated against, so the dose silently landed in a
  sensitivity compartment, every prediction came back `0` and the objective
  function was meaningless.  Since the model is then no longer mixing a solved
  system with ODEs, these fits now warn (recorded in `fit$runInfo`) that the
  analytic `linCmt()` could not be used.  `est="saem"` was never affected, keeps
  the analytic `linCmt()` and does not warn, as do `linCmt()` models with no
  other ODE (#286).

- `est="saem"` no longer estimates a `fix()`ed theta that has no eta attached to
  it; such a parameter now stays at its initial estimate, as it already did for
  the FOCEi family.  The direct phi0 optimization (`nonMuTheta="regress"`, and
  general-likelihood models) takes over phi0 partway through the fit and skips
  the update that restores fixed values, so a fixed non-mu-referenced theta drifted
  off its initial estimate.  Estimates of non-fixed parameters are unchanged.
- `foceiControl(freezeResidGrad=TRUE)` (the default) no longer makes a fit die with
  "maximum number of theta resets (10) exceeded".  The base solve that caches the
  states/EBEs for the frozen gradient ran without the gradient flag set, so an
  ETA-drift theta reset raised inside a gradient restarted the whole fit -- on every
  gradient, until the reset limit tripped (#641).
- A model that combines an inter-occasion variability (IOV) term with a zero
  inter-individual variability eta on another parameter (for example
  `eta.ka ~ 0` alongside `iov.cl ~ 0.1 | occ`) no longer fails with "initial
  'omega' matrix inverse is non-positive definite".  With IOV present the omega
  is a per-condition list, so the zero-eta detector could not read the eta names
  and left the zero eta in the matrix, making it singular; the zero eta is now
  detected and removed as usual.  Restoring the original model after such a fit
  also no longer errors for `est="saem"` (including `table=list(cwres=TRUE)`),
  where the IOV eta is re-expressed as per-occasion id-level etas (#627).
- `est="saem"` no longer collapses subjects that combine two dosing episodes with
  overlapping clock times separated by an `evid=4` reset -- for example a crossover
  where an IV arm and a depot (`f(depot)`) arm share the same times.  SAEM solves
  each subject in the ODE solver's internal time-sorted order, which relocated the
  reset ahead of the first episode's observations and merged the two episodes into
  one trajectory; SAEM then reported a nearly constant `PRED` and a grossly inflated
  residual (`focei`/`posthoc` already handled this correctly).  The reset episodes
  are now offset internally so the solve times increase within a subject, matching
  `rxSolve()`/`focei`; predictions are unchanged because only time-since-reset
  matters (#455).
- The `est="fo"`/`est="foi"` linearization pass returned an intermediate fit
  object with an empty `control`, so `.updateParFixed()` silently fell back to
  default table settings (`ci`/`sigdigTable`) instead of the fit's control
  (#517).  The FO/FOI fit now carries its control, and an intermediate fit
  without a method-specific `nmObjGetControl` surfaces its stored control rather
  than returning `NULL`.
- `est="nlme"` now accepts the common `print` control alias, so
  `nlmixr2(..., "nlme", list(print=0))` no longer errors with
  `unused argument: 'print'`.  `nlme` prints through its own `verbose` option, so
  `print` maps to it (`print=0` runs quietly, any positive value is verbose);
  an explicit `verbose` is still honored when `print` is not supplied.
- FOCEi/FOCE models with a trigonometric term whose argument is a compound
  expression divided by something (for example a sinusoidal enterohepatic-cycle
  release `sin(2 * 3.14 * (time - mtime1) / period)`) no longer fail to build
  with "too few arguments to function 'sin'".  The fix is in `rxode2`'s
  `rxFromSE()` (which was dropping the whole argument, emitting `sin()`); a
  regression test is added here (nlmixr2/nlmixr2est#513).
- FOCEi now estimates a population parameter that is initialized at exactly `0`
  (e.g. a covariate effect or an additive term) instead of leaving it frozen at
  its starting value.  The default scaling constant is `1/|initPar|`, which is
  `Inf` when `initPar` is `0`; it clamped to `scaleCmax` and made the parameter
  effectively unoptimizable.  `getScaleC()` now falls back to unit scaling when
  the initial estimate is `0`.
- A single-subject / fixed-effect ("N of 1") model -- one whose only random
  effects are fixed to zero, which are dropped before estimation -- now gives an
  actionable error when a method that requires random effects (`fo`, `foi`,
  `saem`, `nlme`) is used, pointing to methods that can fit it (`focei`,
  `foce`, or a population method such as `nlminb`, `bobyqa` or `nls`).  The error
  also keeps the user's original model name instead of reporting the internal
  `.mod` (issue #493).
- A focei model whose predictions do not depend on any random effect (for
  example `y ~ dpois(rate)` where `rate` is a fixed population parameter rather
  than a model-predicted value) no longer reports the generic "Aborted
  calculation" message.  The underlying cause is raised directly with guidance on
  linking each endpoint's distribution parameter to an eta-varying model quantity
  (#515).
- `est="saem"`'s "mis-match in nbr endpoints in model & in data" error is now
  actionable: it reports the number of endpoints in the model versus the data,
  lists the observation compartments found in the data, and points the user to
  check that the `CMT`/`DVID` values match the number of model endpoints (error
  terms).  This is the common case of a dataset with extra `DVID` levels that
  the model has no matching endpoint for (issue #579).

- `est="advi"` now rejects a mixture (`mix()`) model up front with a clear
  message (`rxode2::assertRxUiNoMix`) instead of running a wrong fit that ignored
  the mixture structure and then failed late in the output tables with a cryptic
  "the probabilities in a mixture must sum to a number between 0 and 1, they sum
  to: 0".

### Covariance and standard errors

- `setCov(fit, "analytic")` no longer silently installs (and mislabels) the
  `"r,s"` finite-difference covariance when the analytic covariance cannot be
  computed for the model; the fit's covariance is left unchanged instead.

- `fit$etaSE` columns are now labeled `se(<eta>)` (matching `fit$etaRSE`'s
  `rse(<eta>)%`); the label was previously applied to a matrix's `names()`
  (a no-op) so the columns came back as bare eta names.

- `covMethod = "r"`/`"s"`/`"r,s"` standard errors were inflated by a constant
  factor (`sqrt(2)` for `"r"`, `2` for `"s"`) from using `2*R^-1`/`4*S^-1`; they
  now match NONMEM `$COV` (#666).

- A bounded-parameter fit under an unbounded method (e.g. `saem`) leaked the
  internal `rxBoundedTr.<name>` into `$cov` without the back-transform Jacobian;
  `$cov` is now renamed to the original parameters and Jacobian-corrected.

- The analytic FOCE/foce+ covariance no longer falls out of bounds (from dropped
  `eta = 0` solve slots) to the finite-difference Hessian; the general `(f,R)`
  covariance reports `covMethod = "analytic"` (was `"r"`), and
  `foceiCovAnalytic()`/`getVarCov()` reproduce it instead of falling back.

- Fixed a segfault in the analytic covariance for out-of-scope models (the
  augmented build freed the fit's solve before the finite-difference fallback
  ran), and the sign of the M2 upper-tail term in the censored inner gradient.

- The mu-referenced/irls FOCEI-family fits (`mfocei`/`ifocei`/...) now report
  `Condition#(Cov)`/`Condition#(Cor)` in `$objDf`; the post-fit covariance
  install skipped them because the fit tables were rendered before the
  full-model covariance was recomputed.

- Converting a fit to a different covariance (`setCov()`, `getVarCov()`)
  now refreshes `Condition#(Cov)`/`Condition#(Cor)` and the eigen
  diagnostics from the newly installed covariance instead of leaving the
  previous method's values in place.

- SAEM `covMethod = "fim"` adds the mu-block Hessian (was indefinite / NaN SEs),
  and `"fim"`/`"sa"` report off-diagonal Omega and combined residual SEs.
  Fixed `covMethod = "linFim"` and the SAEM covariance erroring for a single
  population/covariate parameter, and `cov2cor` for a one-nonzero-diagonal Omega.

### Estimation and convergence

- A FOCEI fit that hits a theta reset and then restarts no longer aborts with
  `Assertion on 'fitEnv$etaObj$ID' failed: Must be of type 'integer', not
  'factor'`.  The restart re-validated the previous attempt's `etaObf`, whose
  `ID` column is a factor of the original subject IDs; it is now coerced back to
  an integer so a genuinely non-converging fit reports its real reason instead of
  this spurious assertion (#470).

- Fixed the `est = "agq"` quadrature node scaling.  The adaptive Gauss-Hermite
  nodes were placed without the change-of-variable factor, so increasing `nAGQ`
  did not converge to the marginal likelihood -- it converged to a wrong value
  (still better than Laplace, so the objective looked reasonable).  The nodes are
  Gauss-Hermite for the `e^{-x^2}` kernel while the integral has an `e^{-z'z/2}`
  kernel, so they belong at `sqrt(2) * chol(Ht)^-1 * x` with an `exp(x'x)` untilt.
  With the fix the objective converges to the exact marginal likelihood as `nAGQ`
  grows.  Every `nAGQ > 1` objective value (and any standard errors derived from
  it) changes; `focei`/`foce`/`fo`/`laplace` are unaffected.

- The analytic covariance (`covType = "analytic"`) now falls back to finite
  differences under `cholSECov = TRUE`: the covariance step re-factors the eta
  Hessian with the generalized Cholesky, which for a non-positive-definite `Ht`
  differs from the `chol()` the analytic observed information assumes.

- Fixed the `fast = TRUE` analytic gradient for models whose residual variance
  depends on the prediction (`prop`, `add+prop`, `combined1`, `pow`, `add+pow`):
  a determinant chain-rule aliasing injected a spurious term.

- Fixed the `fast = TRUE` analytic gradient/covariance for a random effect
  shared across parameters, enabled sensitivity reuse for a covariate on an
  eta-less parameter, and fixed the gradient never being used live (it read
  finalize-only state and silently fell back to finite differences).

- Fixed the FOCE (`interaction = FALSE`) objective and empirical-Bayes
  estimates: the residual variance is now supplied at the `eta = 0` prediction,
  so ODE and `linCmt()` FOCE agree and match the NONMEM reference.

- Bounded the Shi (2021) finite-difference step so a curvature-free search can
  no longer corrupt the shared solver state.

- Fixed `muModel = "lin"`/`"irls"` erroring with two or more covariate
  expressions (#711) and the user-fixed covariate-coefficient regression bias.

- Fixed `impmapControl(impSeed = )` being ignored.

- FOCEI now updates additive mu-referenced population parameters with
  large-magnitude initial estimates (#641).

- FOCEI theta resets now keep every reset population parameter inside its
  bounds instead of restarting the optimization out of range, and stop with an
  informative error when a parameter's bounds are infeasible (#454).

### Crashes and stability

- Fixed a Windows heap-corruption segfault at more than one core (rxode2 saw
  every worker as thread 0); the inner loops now pass the real thread id.

- Fixed a segfault in `est = "vae"` (thread count capped at the solve's core
  count) and in `nlmSetup` on the first estimator call of a session.

- Fixed FOCEi aborting with `Cube::slice(): index out of bounds` when
  `mceta >= 1` and `maxInnerIterations == 0`, and a heap-buffer overflow / wrong
  back-transform in SAEM Box-Cox residual models.

- A non-positive-definite `Omega` is projected to the nearest PD matrix (SAEM
  mid-run, with a `fit$runInfo` warning; and the sym-inv-chol setup for a
  degenerate fit) so residual/table diagnostics still run; NPDE with a
  degenerate simulated covariance sets the subject's NPDE to `NA` instead of
  aborting.

- Fixed a segfault when a dataset has no observed subject at all (every subject
  is a placeholder with no `EVID==0` row, as in an aggregate-data output eval
  such as `babelmixr2`/`admixr2`).  The no-observation-subject drop now keeps
  the rows when there is no observed subject to fall back to, and `foceiSetup_`
  no longer reads an empty id vector out of bounds.  `.nlmSetupEnv()` also now
  supplies a default `iterPrintControl` when an external caller omits it,
  instead of erroring with `Index out of bounds: [index='iterPrintControl']`.

### Output, tables, and printing

- For models without etas, the `BSV(SD)` and `Shrink(SD)%` columns are no longer
  added to `$parFixed` and `$parFixedDf`; they were always blank for these models
  (#355).

- Model-defined variables (e.g. `ka`, `cl`, `v`, `tad`, `dosenum`, and any
  user-added line such as `WT.OUT <- WT`) are now included in the output table
  whether or not `cwres` is requested.  Previously `tableControl(cwres=FALSE)`
  dropped these columns while `cwres=TRUE` (the default) kept them, so the same
  model produced different output columns depending on the residual request
  (#497).

- A zero-fixed eta (e.g. `bsva ~ 0`) is again restored into the fitted model's
  `ini()`/`model()` blocks when the estimation makes a nested `nlmixr2()` call
  (e.g. adding the focei objective or CWRES), so `fit |> ini(bsva ~ 0.1)`
  works; the nested call used to wipe the restore info held in a global (#741).

- `augPred()` now works on a `focei` fit whose model has a zero-fixed eta that
  appears in the prediction (e.g. `eta.v ~ 0` used in both the ODE and the
  residual), instead of erroring with `parameter(s) are required for solving:
  eta.v`; the simulation model drops the zero eta consistently with `saem`
  (#514).

- `laplace`/`agq` family fits label their `$objDf` row `Laplace`/`AGQ<n>`
  (matching `$ofvType`) instead of `FOCEi`; previously the default
  `interaction=TRUE` made the interaction label win over the quadrature one.
  The quadrature objective stays the active one after CWRES; `setOfv(fit,
  "focei")` (and `addCwres()`) now evaluate the true focei objective on a
  quadrature fit instead of re-labeling its quadrature value.

- Restored the `Function Val.` objective column and the `$parFixed` shrinkage
  coloring; periodic headers now repeat only the column labels.

- `$parFixed` honors a user `sigdig`/`ci` for fits with literally-fixed
  parameters.

- Literally-fixed population parameters now report their back-transformed value
  (`exp`/`expit`/`probitInv`) in the `Back-transformed` column instead of the
  raw log/logit-scale estimate.

- `augPred()` now keeps the fit's original subject ids: the returned `id`
  factor carries the actual (character/factor) ids from the fit instead of the
  internal integer re-numbering (#450).

- `vpcSim(fit, pred=TRUE)` (and hence VPC plots with a `pred` line) now works
  for models with IOV.  With IOV the fit's `omega` is a list of matrices (`id`
  plus one per occasion level), which the `pred` path treated as a single
  matrix and errored with `invalid 'times' argument`; the population prediction
  now zeros every random effect across all omega levels (#629).

- `fit$time` again attributes model build/compile to `setup`/`configure` (and
  the nlm family times setup/optimize) instead of `other`.

- Aggregated ODE-solve warnings report the real subject id; `parHistData` shows
  mixture-probability parameters on the natural scale and `fit$mixList` returns
  all components; iteration printing labels the estimation phase (`Burn in`/
  `KL anneal`/`EM`/`Smooth` for vae, `SA`/`EM` for saem).

- `fast = TRUE` with a `linCmt()` model downgrades to `fast = FALSE` with a
  message instead of silently falling back per gradient call.

- `est = "vae"` with automatic covariate selection now reports the selected
  covariate coefficients (`beta_<par>_<cov>`) in `$parFixed`/`$parFixedDf`
  instead of dropping them when a population parameter is fixed, and the
  covariate-bearing mu-parameters back-transform (`exp`) instead of printing on
  the raw log scale.

- `est = "vae"` no longer errors with `cannot find parameter 'NA'` when a
  structural (mu-referenced) parameter is fixed with `fix()`; its random effect
  is kept (variance estimated) with the fixed value carried in the model.

### Data handling

- SAEM no longer errors with `No data with ID` for a dose-only subject;
  observation-less subjects are dropped before estimation and re-inserted into
  the output with a population `PRED` and `NA` individual columns, like FOCEi
  (#687).

- FOCEi no longer errors with `'names' attribute [n] must be the same length as
  the vector [m]` when a subject's records are all removed during data
  translation (e.g. every `TIME` is `NA`).  Such a subject vanishes from the
  processed data entirely rather than losing only its observations, so it is now
  detected and dropped from the subject index alongside observation-less
  subjects (#606).

- Fixed `nlmControl()` listing `eventSens`/`sensMethod` twice.  The "initial
  ETAs were nudged" warning fires only when a nudge actually happened, and a
  non-default `mceta` on a fully mu-referenced model falls back to the default
  with a warning.  `saemControl(covMethod = "")` (skip covariance) no longer
  errors.

### Internal

- Removed an unreachable duplicate `missingTable` default assignment in
  `nlmixr2Est0()` (issue #385); the earlier default already fixes the value, so
  the second block could never run.  No change to fit results.
- Removed the last bare `Rf_error` call from the C++ sources (issue #632):
  the `Rcpp::compileAttributes()` output now emits the parenthesized
  `(Rf_error)` form, and the internal `rxError` macro was switched to
  `(Rf_error)` as well, so the package no longer trips Rcpp's upcoming
  `Rf_error` deprecation warning (RcppCore/Rcpp#1247).  The C `.Call`
  entry-point validators keep their justified `Rf_errorcall` uses.

- Consolidated data preparation and the nlm-family control/fit functions, and
  the analytic-covariance augmented model now uses rxode2's chunked
  `rxOptExpr()`; no change to fit results.  The test suite runs a single
  testthat worker on CI/CRAN and parallel elsewhere, with within-solve threads
  capped to 2 on CRAN.

# nlmixr2est 6.1.0

- Added focei, foce, foi, fo mixture support in `nlmixr2est`

- Fix `focei` mixture models with llik residual distributions erroring
  when a model had exactly one mixture probability parameter

- Fix `fit$mixList` returning only the first mixture component

- `parHistData` Back-Transformed rows now show mixture probability
  parameters on the natural probability scale (0, 1) instead of the
  raw mlogit estimation scale.
- Fix issue 641: FOCEI now updates additive mu-referenced population
  parameters whose initial estimates are large in magnitude.
  Previously a missing branch in `.foceiOptEnvSetupScaleC()` let
  `scaleC` fall through to the C++ default of `1/|init|`, which mapped
  unit steps in scaled space to negligible steps in unscaled space and
  effectively pinned such parameters at their initial value (e.g.
  `tvemax <- -40` with no transform).
- When model estimation fails, all errors raised during the run are now
  collected and reported together, instead of only the last error. This
  is supported by a new `collectErr` argument to the internal
  `.collectWarn()` helper, which captures errors alongside warnings and
  returns them in the `error` element of its result list. As a result,
  errors hidden by `on.exit({rxode2::rxProgressAbort()})` handlers
  (such as the "Aborted calculation" message reported in issue 607)
  no longer mask the underlying cause; both the inner stop message and
  any follow-up error from `on.exit` are now reported to the user.
  parameters on the natural probability scale instead of the raw
  mlogit scale.
  parameters on the natural probability scale

- Hardened mixture-model (`mix()`) estimation: clearer errors for
  `est="nlme"` and invalid initial probabilities, warnings for
  underflowing/collapsing mixture probabilities, and a fix for the
  SAEM omega-diagonal floor being raised outside mixture fits

- Fix segfault in `nlmSetup` on the first estimator call of a fresh R
  session for pooled estimators

- Guard against null pointer arithmetic in inner.cpp

- Use OpenMP threading for S matrix calculation

- Use OpenMP threading while calculating NPDEs

# nlmixr2est 6.0.1

- Fix LTO violation as requested by CRAN by adding
  -DARMA_DONT_USE_OPENMP to PKG_CXXFLAGS in src/Makevars.in

- Require rxode2 5.1.2 which has the fixed M1-san issues observed
  here.


# nlmixr2est 6.0.0

- `focei`, `foce`, `fo`, `laplace`, and `agq` have all been
   successfully made thread safe and parallelized (for a single
   CPU). The default tolerance relaxation for difficult to solve ODEs
   has been changed to per individual instead of for the entire
   population (which is a breaking change, so major release).  This
   should allow more precision for a majority of the subjects in the
   optimization process.

- Add `predict(fit, level="ipred")`, `predict(fit,
  level="individual")` or `predict(fit, level=1)` to predict
  individual fits (with possibly a new dataset).

- Change test files to `.rds` files

- Drop magrittr `%>%` in favor of `|>`.

- **Breaking change:** Minimum R version increased from 4.0 to 4.1.0.
  This change is required to support the native pipe operator `|>`.
  Users on R < 4.1.0 will need to upgrade R to install this version
  of nlmixr2est.

- Bug fixes for deparsing nlmixr2 control objects

- `nlm` and related pooled methods now run in parallel (based on ID)

- Tests are optimized to reduce redundant fits and run in parallel.

- `nlm` (and related pooled optimizers: `bobyqa`, `newuoa`, `uobyqa`,
  `n1qn1`, `lbfgsb3c`, `optim`, `nlminb`) now support the same
  censoring behavior (M2/M3/M4) as FOCEI and SAEM.  The
  `$censInformation` field is populated for these fits in the same way
  as FOCEI/SAEM.

- `agqControl()` and `laplaceControl()` now have `rxUiDeparse()`
  methods so they can be saved better in packages like `nlmixr2save`
  and `shinyMixR`.

- Added new `outerOpt`; methods to `focei` and related methods (`agq`,
  `laplace`, `foce`, `fo`, `foi`): "uobyqa" and "newuoa".

- `saem` and other methods now respect bounds by default by internally
  adding the appropriate transform and then applying the
  back-transformation just before returning.

  For parameters that are mu-referenced, this breaks
  mu-referencing. When it breaks mu-referencing there is a warning
  issued.  The best practice is still to have unbounded parameters
  with mu-referencing.

  If you want to ignore this behavior you may
  use `control=list(boundedTransform=FALSE)` or for saem
  `control=saemControl(boundedTransform=FALSE)`

- The mu referencing covariate procedure was made less fragile to
  support mu referencing in conjunction with iov and bounded parameter
  transformations.

- Add some bench-marking capabilities and small speed fixes for focei/saem

# nlmixr2est 5.0.0

- Remove `qs` and change to `qs2`.  This breaks backward
  compatibility.

- Default to non-compressed nlmixr2 objects

# nlmixr2est 4.1.1

- Request nlmixr2est's pre-processing hooks for `augPred()`, `vpcSim()` and
  `$simInfo`, which fixes augPred in cases where `etas=0` are used in
  `nlmixr2` (#587)

- Fix scale.h so that `scaleType="none"` does not also require
  `scaleTo=0`

- Request Armadillo 15 with the special flag in the new `RcppArmadillo`

- Fix `focei` without etas (and without log-likelihood normal) to run
  `ELS` (See #590).

- Change the IOV implementation (#596):
   - Now shows estimates as `CV%` or `sd` without shrinkage calculation.
   - Allow different forms of `iov` estimation, controlled by
     `iovXform`.
   - Retains the `iov` parameter(s) in the output `data.frame`.
   - With `iov`, the `$omega` shows a list of variability by the
     conditioning variable(s).
   - `fit$iov` will show the IOV deviations by the conditioning
     variables(s) with the exception of `id`
   - IOV models can be used in other estimation methods and inherits
     the ETA values.

 - Added `$etaMat` method for `nlmixr2` fits to give the value that
   needs to be passed between each estimation method (related to iov #596)


# nlmixr2est 4.1.0

- Updated inferring the estimation method from the control
  object. Requires the control object to have a class of length one
  and match the estimation method.  For example `foceiControl()` would
  assume that the estimation method is related to `focei`.

- Changed Rstudio completion to not evaluate (in case it gets turned
  on for data.frames) (See #568)

- Turned on data completion for items like `$fitMergeInner`

- **Breaking change:** Changed the estimation method `posthoc` to add
  tables and calculate the covariance by default.  It is now a method
  with it's own control, `posthocControl()`.  As previously the
  default is not to include the interaction term (but you can turn it
  on with `posthocControl(interaction=TRUE)`).

- Added `foceControl()`, `foControl()` and `foiControl()` for the
  `foce`, `fo` and `foi` methods, respectively.  They try to convert
  the related control structures to the correct control structure for
  the estimation method.

- Added iov support for `focei`,  `foce`, and `saem` (#614)

- Added new estimation method `agq` which uses adaptive Gauss-Hermite
  Quadrature to fit a nonlinear-mixed effect model. In this method,
  you can choose the number of quadrature points to estimate the
  likelihood, with higher numbers giving more accurate likelihoods.
  The AGQ implementation in nlmixr2est allows you to specify the
  number of quadrature points via the `agqControl()` function, and
  supports both single and multiple subject models. This method is
  particularly useful for models where accurate likelihood estimation
  is critical.

- Also added a `laplace` method which is the same as
  `agq` with 1 node (and is numerically the same as `focei`, `foce` or
  log-likelihood `focei`/`laplace`, etc), but uses the `agq` routine.

- Fixed saem mu-reference display by not compressing the internal item
  `saem0`.

# nlmixr2est 4.0.2

- The loading and unloading of DLLs has been minimized in this version
  of nlmixr2est. This avoids loading/reloading the same DLLs and causing the
  CRAN mac m1 ASAN/USBAN false positive issue observed in CRAN.

- Additionally a new function `nlmixr2fix(fit)` has been added to
 `nlmixr2est`.  It attempts to make the fit loaded from a different
 version of nlmixr2 compatible with nlmixr2 4.0.  It also prints out
 the versions of `nlmixr2` that were used when creating this fit.
 With this information you are more likely to find a way to use the
 fit in your current session (or in an old session). (Issue #562)

# nlmixr2est 4.0.1

- Initialize lbfgsb3 error message to an empty string to address
  valgrind finding (as requested by CRAN).

# nlmixr2est 4.0.0

- When using a model to start a new focei model, the ETAs from the
  last fit are used as the starting point.  Now you can use
  `foceiControl(etaMat=NA)` to skip this and use `eta=0` for all
  items.

- When using `foceiControl(etaMat=fit)`, this will extract the ETAs
  from a fit for use in the next optimization.

- When using a `foceiControl(etaMat=)` option nlmixr2 no longer only
  evaluates the inner problem with the `etaMat` value.

- Add `mceta` option to `"focei"`.

  - `mceta=-1` is the default; the eta restarts at the best eta from
    the last step to start the inner optimization.
  - `mceta=0` the eta starts at `0` to start the inner optimization.
  - `mceta=1` the eta starts at either `0` or the best `eta`, which
    ever gives the lowest objective function to start the inner
    optimization.
  - `mceta=n` under the assumption of `omega` sample `n-1` `eta`
     values and use the lowest objective function of eta sampled, last
     best eta and eta=0 to start the inner optimization.

- Fix Rstudio print (issue #536)

- Support rxode2's new `+var()` definition in `saem`

- Support literal fixing of residuals (#524).  All methods that
  support a literal fix of residuals have an option `literalFixRes`
  which defaults to `TRUE`.  To get the behavior from older models you can use
  `literalFixRes=FALSE`
- More detailed error messages will be reported for models with errors

# nlmixr2est 3.0.4

- More robust covariance calculation in `focei`.

- Allow hook mechanism to handle piped arguments.

- Fix for when output message from optimizing doesn't print well
  (#325)

# nlmixr2est 3.0.3

- Moved data check for covariates and required data items to a
  pre-processing step. This fixes #499.  Each method that needs to
  have a covariate check needs to have a property `covPresent`. For
  example to apply the covariate data check to the `focei` method you
  need `attr(nlmixr2Est.focei, "covPresent") <- TRUE`.

- Bug fix for non-mu referenced etas when combined with mu referenced
  covariate values. (See #498)

- Changed option for `"saem"` to have `literalFix=FALSE`. This makes
  mu-referencing work better when fixing a population value.

# nlmixr2est 3.0.2

- Fix bug where models where omega boundary warnings caused problems
  in estimation (#490)

- Created a new api for pre-processing ui, allowing adding arbitrary
  hooks.  As written now, this includes literal fix and zero omega as
  well as added the new rxode2 ui processing.

- Fixed compilation to only use -I in most systems for maximum
  compatibility

# nlmixr2est 3.0.1

## New features

- Now when optimizing only a single parameter with `focei`-family,
  will change to use `stats::optimize()` for the outer problem (#481)

- When estimating with all fixed population parameters, do a posthoc
  estimation.

- Internally removed `assignInMyNamespace()` replacing with
  `nlmixr2global`, which fixes some edge case bugs where the nlmixr2
  environment was not reset properly.

- Treated edge case where all initial parameters are zero and change
  scaling from scaled to unscaled (#486)

- Added `mu`4 referencing that will change string expressions to
  `rxode2` numeric values.  This allows derived strings to also be
  treated as `mu` expressions (#484)

## Bug Fixes
- Fix `focei` covariance step when many `omega` values are fixed #482

# nlmixr2est 3.0.0

- No binary linking to `rxode2`, `lbfgsb3c` and `n1q1`, which means
  that updating these will not make `nlmixr2est` crash without
  recompiling.

- New `mu`3 referencing will take context from the model to see if the
  algebraic expression can be completed from defined model variables;
  These variable would have to be unique.

# nlmixr2est 2.2.2

## Breaking changes

- Saem non-mu reference input parameters/covariates were fixed so they
  work correctly with fixed parameters (Issue #445)

- Focei changed back to having a lower bound for standard deviations
  when not specified. This means that best model fits may change.  You
  can revert to the old settings by using
  `foceiControl(sdLowerFact=0.0)`.  You can also change the factors to
  other values than the default value, that is
  `foceiControl(sdLowerFact=0.000001)` for instance which would
  multiply the initial value by `0.000001` when either the lower bound
  isn't specified or the lower bound is specified as zero for the
  error estimates related to error-based standard deviations.

- In `nlmixr2`, expressions are optimized.  Because of that
  optimization, numerical rounding differences can cause different
  directions in optimization when fixing parameters in the model
  vs. fixing the parameters manually.

  This means that the fixed parameters in a model vs hard-coded fixed
  parameters could give different values in the final model.

  A new option `literalFix` was introduced which change the fixed
  population parameters to constants in the model while running the
  optimization.  This makes the output of fixing within the model and
  fixing manually the same (which is what is likely expected). The
  default is for this to be turned on (ie. `literalFix=TRUE`).  You
  can get back the old behavior by using the option
  `literalFix=FALSE`.

- In `saem`, the monte-carlo sampling occurs for all parameters
  including non-informative ETAs.  A fix ensure that non-informative
  etas in `saem` are fixed to zero while sampling the `phi` values.
  This may change results for models with uninformative etas. To
  ignore the uninformative etas with `saem` you ca use use the prior
  `saem` handling with `saemControl(handleUninformativeEtas=FALSE)`.

## New features

- Gracefully degrade when $cov is not in the right form (see #423)

- Add support for PopED in place solving (used in babelmixr2)

- If `est=foceiControl()` or other nlmixr2 control with the class
  `foceiControl` infer the estimation method is `focei`

- Add back the warnings when estimation methods ignore the boundaries

- When using `rxSolve`, now respects the values from `tableControl()`
  (#465 and #297)

## Bug fixes

- Will emit warnings when the return object is not a nlmixr2 fit
  (#453)

## Other things

- Moved actual code of some matrix libraries to `lotri` and import
  them via function pointers

# nlmixr2est 2.2.1

- Align with the possibility that linCmt sensitivities may not be
  present (like intel c++)


## Bug fix
- `focei` cache needs to be based on the parameter order as well as
  the model information (#415)

# nlmixr2est 2.2.0

## New Features

- Algebraic mu referencing has been implemented in `nlme` and `saem`.

- New estimation method "nlm" has been added to estimate population
  only likelihoods using `stats::nlm` and possibly return a
  standardized `nlmixr2` fit.

- New estimation method "nls" has been added to estimate population
  only problems.  This uses `minpack.lm::nlsNM` by default if
  present, or the `stats::nls`

- New estimation method "optim" has been added to estimate population
  only likelihoods.  This uses `stats::optim` and returns a
  standardized `nlmixr2` fit.

- New estimation method "nlminb" has been added to estimate population
  only likelihoods.  This uses `stats::nlminb` and returns a
  standardized `nlmixr2` fit.

- New estimation methods from the `minqa` package: "bobyqa", "uobyqa"
  and "newuoa" have been added to estimate population only
  likelihoods.  These methods returns a standardized `nlmixr2` fit.

- New estimation method "lbfgsb3c" to estimate population only
  likelihoods.  This returns a standardized `nlmixr2` fit.

- New estimation method "n1qn1" to estimate population only
  likelihoods.  This returns a standardized `nlmixr2` fit.

- Added new feature for `vpcSim()` where a minimum number of subjects
  are simulated from the model when trying to fill in ODEs that were
  not solved successfully.  By default this is `10`.  This also
  works-around a bug when there is only one subject simulated and the
  `data.frame` has a slightly different output.

## Breaking changes

- Removed `fit$saemTransformedData` since it isn't actually used in
  `saem` anymore (but will break anyone's code who is using it)

- Now the internal function `.foceiPreProcessData()` requires the
  rxode2 control `rxControl()` because some of the new steady state
  lag features need to translate the data differently based on
  `rxControl()` options.


## Bug fixes

- Printing models with correlated omega values and omega values fixed
  to zero no longer fails (#359)

- Add back values for $parHistData (#368)

- This requires a new `rxode2` which will fix multiple endpoint issues observed (#394)

- Manual back-transformed values in `$parFixed` are now displaying
  correctly and are calculated based on the confidence interval in the
  control instead of 95% confidence no matter what (#397)

## Other changes

- An `as.rxUi()` method was added for fit models (#377)

# nlmixr2est 2.1.8

- Version bump and a minor documentation update (same as nlmixr2est
  2.1.7).  This version bump is to simply allow correct binary linkage
  to rxode2 2.0.14. Otherwise `nlmixr2` models will crash R.

# nlmixr2est 2.1.7

- As requested by CRAN, remove `Rvmmin`

- Values in `$parFixed` for BSV without exponential transformation are now
  correctly shown (#366)


# nlmixr2est 2.1.6

## Breaking changes

- Since `rxode2` now allows simulation with `omega` having diagonal
  zero elements, `$omega` and `$omegaR` now reflects this information
  including the zero omega elements in the output. On the other hand,
  the other eta-information and standard error information for zero
  etas are still excluded in `$phiR`, `$phiSE`, `$eta` etc.

## Bug fixes

- `vpcSim()` works when an eta value is fixed to 0 (#341)

- `augPred()` now consistently uses the simulation model (instead of
  the inner model used for `CWRES` calculation).

## Other changes

- Dropped dependence on orphaned package `ucminf`

# nlmixr2est 2.1.5

- Add `$fitMergeFull`, `$fitMergInner`, `$fitMergeLeft`,
  `$fitMergeRight` as a complement to `$dataMergeFull`,
  `$dataMergInner`, `$dataMergeLeft`, `$dataMergeRight`.  The fit
  variants prefer columns in the fit dataset instead of the original
  dataset.  This is useful for goodness of fit plots with censoring
  since the `DV` in the fit simulates values under the ipred/residual
  assumption and will give more appropriate goodness of fits,
  otherwise these values are the limit of whatever censoring is
  applied

- Moved the mu reference fix for the split mu referenced model here
  (from babelmixr2)


# nlmixr2est 2.1.4

- Breaking change, now calculate condition number based on covariance
  and correlation, the names have changed to be more explicit.
  `conditionNumber` changed to `conditionNumberCov` and a new metric
  `conditionNumberCor` has been added.

- A bug in boundary value detection prevented automatic covariance calculation
  with FOCEi estimation (#318)

- Fix `vpcSim` so that it will be a bit more robust when it is
  difficult to simulate.

- A bug in model piping which did not allow models to be appended to was fixed
  (rxode2#364)

- An internal change was made in `nlmixr2.rxUi()` to better support the
  babelmixr2 PKNCA estimation method (babelmixr2#75)

- Fixed bug where `$iniUi` did not return the initial ui when running
  non `focei` related methods.  Also added alias of `$uiIni` to the
  same function.

- Dropped Stan headers for this package, also updated to C++17

# nlmixr2est 2.1.3

- Allows `$etaH` and related family to be integrated into a `saem` fit
  if `cwres` is calculated.

- Fixed a bug where `nlmixrLlikObs` in the merged dataset is sometimes
  named `llikObs`, now it is always named `nlmixrLlikObs`

- Fixed a bug where `nlmixrLlikObs` shows up in merged dataset when
  `cwres` is not calculated (it was always `0`), also allow `cwres`
  calculation to pick up `nlmixrLlikObs` in merged dataset.

- Dropped `dparser` dependency

# nlmixr2est 2.1.2

- Fixes `$etaH` memory corruption so the standard errors of etas are now correct

- Removed the memory requirements for focei by `neta*neta*nsub`

- Fixed character based covariates so the work correctly (again) with
  focei.  Added a test for this as well.

# nlmixr2est 2.1.1

- Fixes `$dataMergeInner` so that observation-based log-likelihoods
  work with infusions.  Should fix tests with `ggPMX`

- Fixes `$etaSE` and `$etaRSE` to work correctly when there is only 1
  eta.

- Fixes npde valgrind observed on CRAN machines

# nlmixr2est 2.1.0

## Breaking changes

### FOCEi

 - Gill forward differences will not repeat now (by default), You can
   change back to prior behavior with `foceiControl(repeatGillMax=3)`

 - Number of sticky recalculation is reduced to 4; to have the old
   behavior use `foceiControl(stickyRecalcN=5)`

 - `n2ll` has been changed to `ll` to specify individual
   log-likelihoods.  This was only used in simulation and was not well
   documented.

 - Generalized log-likelihood is only supported with `rxode2` `2.0.8` or later.

### FOCEi covariance calculation

 - The `S` matrix calculation was made a bit more robust to errors in
   individual gradients.  When there are errors in the individual
   gradient calculation, assume the gradient is the same as the
   overall gradient.  In the tests cases, were reasonable using this
   adjusted S matrix.  This means if some individuals do not have very
   much data to support a specific parameter, a `S` matrix calculation
   for the population will still be generated. When there is some
   patients/subject combinations that do not have sufficient data, we
   will add the following to the run information: `S matrix had
   problems solving for some subject and parameters`. The `S` matrix
   calculation will still fail if the percentage of parameters that
   are being reset is lower than `foceiControl(smatPer=0.6)` or
   whatever you specify.

 - The `r,s` covariance matrix will now also check for unreasonably
   small values (controlled by `foceiControl(covSmall=...)`) and
   select a different covariance estimate method even when the "r" and
   "s" matrices are calculated "correctly".

## New features

- What type(s) censoring (if any) is now stored in `fit$censInformation`

- Standard errors of `$etas` can now be obtained with `fit$phiSE`,
  also available are `fit$phiRSE` (relative standard error),
  `fit$phiH`, (individual hessian), `fit$phiC` (individual
  covariances), `fit$phiR` (individual correlation matrices)

- Can also use Shi 2021 differences in addition to Gill differences.
  In our tests (using the same datasets as CPT) these produced worse
  estimates than the Gill 1983, though it is unclear why since it
  should be a faster more accurate method.  A modified version is used
  in calculating the individual Hessians of numerically for the
  generalized likelihood approach.

- Generalized likelihood estimation is now present in `nlmixr2est` for
  `focei`, `foce` and `posthoc`

- `nmNearPD()` is a function you may use for nearest positive definite
  matrix.  This is derived from `Matrix::nearPD()` but is implemented
  in C/C++ to be used in (possibly threaded) optimization.

- Individual Hessians can be accessed by `$phiH`, covariance by
  `$phiC`, eta standard errors by `$phiSE` and eta RSEs can be
  accessed by `$phiRSE`.  There are `eta` aliases for these as well
  (`$etaH`, `$etaC`, `$etaSE`, and `$etaRSE`).

- Can now access the individual point's contribution to the overall
  likelihood when merging to the original dataset. These merges can be
  accessed with `$dataMergeFull`, `$dataMergeLeft`, `$dataMergeRight`,
  and `$dataMergeInner`.  The columns with the individual data column
  is `nlmixrLlikObs`.

  To calculate the total `focei`/`foce` objective function, the sum of the
  likelihoods still need to be adjusted by the omega/eta contribution,
  and the individual Hessians, and possibly the NONMEM objective
  function offset constant.

## Censoring fixes

 - Fixed bug where datasets with censoring that are not lower case `cens` and `limit` do not
   produce the correct table output (#180)

## FOCEi updates

- Resets now scale properly when a value is simulated outside the limit
- Models with zero gradients on the first step now switch to `bobyqa`
  by default.  With this, it is more important to examine the model
  parameters and fits for plausibility.

# nlmixr2est 2.0.8

## New features

- Add `pd`/`npd` as an output as well as `npd`/`npde`

## SAEM bug fix

- When loading a `nlmixr2` "saem" fit from another R session,
  `nlmixr2` will no longer crash with `fit$objf`

## NPDE/NPD fixes

- `NPDE` was identical to `NPD` even with correlated models, this was
  fixed (prior output was actually `NPDE`).

## Censoring fixes

- FOCEi censoring fixes:
  - M4 method equation bug fix
  - M4 method derivative change based on equation fix
  - M2 method added missing derivative
  - Censoring already dTBS

- SAEM Censoring fixes:
  - SAEM method M4 method equation bug fix
  - Censoring limit changed to dTBS

- Censoring handling was unified

## Internal changes

- Added `ui$getSplitMuModel` which is used in `babelmixr2` and will be
  used in the refined stepwise covariate selection of `nlmixr2extra`

- Added work-around to remove `_nlmixr2est_RcppExport_registerCCallable`
  since the registering of C callable are handled manually at the moment.

# nlmixr2est 2.0.7

- Use `.zeros()` for the matrices in armadillo in addition to relying
  on `calloc` to give zero matrices.

- Fixed one uninitialized object

- Fix for `augPred` so it works on population only models

- `nlme` no longer sets options to treat all covariates as non
  mu-referenced covariates, but directly calls a function that can
  turn on or off the mu-reference covariate selection.

- `vpcSim` now tries to simulate IDs that didn't simulate correctly (with a warning)

- Export nmObjHandleControlObject

# nlmixr2est 2.0.6 -- new package

`nlmixr2est` contains the estimation functions within `nlmixr2`.

## FOCEI family changes

- Remove lower level `foceiFit` function.  Focei, foce, fo, foi, and
  posthoc now directly takes rxode2 ui objects

- New error types are supported in focei including mixing theta and
  etas in residual errors and different types of proportional errors

- Different types of additive and proportional errors can be used for
  each endpoint using ` + combined1()` or `+ combined2()` otherwise it
  takes the supplied `addProp` option to figure out which type of
  combined model is run (by default `combined2()`)

- Focei model cache is now named `focei-md5Digest.qs` and uses `qs`
  compression/saving/loading.

- `foceiControl()` aligned between other methods.

- `foceiControl(adjLik=TRUE)` uses the NONMEM-style objective function
  throughout.  `foceiControl(adjLik=FALSE)` uses the adjusted
  objective function throughout, and adjusts it back to the NONMEM
  objective function.

- Lag time and other between subject variability differences no longer
  calculate an ideal relative step size, but an absolute step size
  when using Gill differences (default)

- Objective function checks for infinite/NaN/NA values for the entire
  solving space and ensures no overflow occurs when calculating the
  inner hessian

## SAEM changes

- mu referencing is no longer required for `saem`; Internally non
  mu-referenced values are converted to mu referenced values and the
  converted back when calculating the nlmixr2 object.

- `nlmixr2` forced the parameter ordering to (1) population effects,
  (2) non mu-referenced between subject effects (3) omega estimates
  and (4) residual effects. This changes the order that `nlmixr2` sees
  the parameters. Since this is based on a random number generator,
  the optimization trajectory will be different and have different
  results than `nlmixr`

- Components of `omega` can now be fixed.

- Residual error components can also be fixed.

- When optimizing only one residual value, nlmixr2's saem uses `nlm`
  from R, which is more efficient than the nealder-meade method.

- Lower level `saem` functions (like `configsaem()`) are not exported
  because they are increasingly difficult to use and convert to
  something standard; a few methods (like `print`, `summary` etc) are
  maintained to view the lower level object and for debugging it.

- Parameter history and print-out no longer includes fixed parameters.

- The model to calculate the residuals more closely matches the model
  used for estimation to remove small rounding differences that may
  occur in the models.

- Different types of additive and proportional errors can be used for
  each endpoint using ` + combined1()` or `+ combined2()` otherwise it
  takes the supplied `addProp` option to figure out which type of
  combined model is run (by default `combined2()`)

- Parameter history and printout now uses standard deviation for
  additive only components, matching the estimation of the components.

- `rxode2` solving options are now saved in the `rxControl` part of
  the `saemControl()`.  That is
  `saemControl(rxControl=rxControl(...))`; This fixes any conflicting
  option names as well as allowing alignment between the control
  structures in `focei`, `nlme` and `saem`

- `saemControl()` aligned between other methods.


## nlme changes

- `nlme` has been completely rewritten to directly run from the
  `rxode2` UI

- `nlme` always tries to use mu-referencing (when available)

- Internally `nlme` now uses parallel processing for solving so it
  should be faster.

- `nlmixr2NlmeControl()` (which will overwrite `nlmeControl()`)
  documents and adds more options to `nlme`. Also aligned with other
  methods.

- `weights`, `fixed`, `random` can be specified in
  `nlmixr2NlmeControl()`.  If so, then the `nlme` object will be
  returned.

- `returnNlme` is a new option that will return the `nlme` object
  instead of the traditional `nlme` object.

- `nlme_ode` and `lme_lin_cmpt` are both removed.

- `rxode2` solving options are now saved in the `rxControl` part of
  the `saemControl()`.  That is
  `nlmeControl(rxControl=rxControl(...))`; This fixes any conflicting
  option names as well as allowing alignment between the control
  structures in `focei`, `nlme` and `saem`

## nlmixr2 object change

- With `saem`, the nlmixr2 function now saves/compresses the `phiM`
  information.  This means the gaussian and Laplacians likelihoods can
  be calculated when you save the nlmixr object and then restore it
  later.

- The nlmixr2 object compresses infrequently used and removes many
  unneeded objects. Even with compression, the `saem` objects are
  often a bit bigger since they include the large `phiM` object.

- `nlmixr2` now supports non-mu referenced ETAs in the `fit$parFixed`
  and `fit$parFixedDf`

## nlmixr2 interface change

- `nlmixr2` interface changed to use `rxode2` UI

- `keep` and `drop` are added to `tableControl` to influence the end data-frame

- `$simInfo` uses a quoted expression for `$rx` instead of a string

- `$simInfo$sigma` is a diagonal matrix since now the normal
  simulation is controlled by the variability modeled as a population
  value.

- `nlmixr2` now allows etas that have initial omega estimates of zero
  to be dropped from the model (instead of issuing an error about a
  non-positive definite `$omega` matrix)

## NPDE changes

- Fixed a bug where the number of simulations for a NPDE calculation
  are correctly passed by `addNpde(fit, table=tableControl(nsim=500))`

## VPC changes

- `vpc` function rewritten and split out to `vpcSim()` and
  `vpcPlot()` (which is a replacement for `vpc()`).

- There were too many mismatches between `vpc::vpc` and `nlmixr::vpc`
  which caused inconsistencies in code based on load order of `vpc`
  and `nlmixr`.  This way both coexist, and you can use the `vpc`
  simulation for other packages more easily (like `ggPMX`) without
  creating or summarizing data since `ggPMX` has its own methods for
  summarizing and creating plots.

- VPC now directly uses `rxode2::rxSolve`

## augPred() changes

- `augPred()` has been written to use the new fit object.

- `nlmixr2AugPred` was changed to `nlmixr2AugPredSolve()`

- `augPred` uses the new interface and supports multiple endpoints.
  The endpoint name is now always on the `plot(augPred(fit))`.

## getFitMethod() change

- Internally, fit estimation method is saved in `fit$est`, and now
  `getFitMethod(fit)` simply returns `fit$est`

## Delete methods

- Many methods lower level utility functions have been deleted.

- `nmDocx`, `nmLst` and `nmSave` have been removed.

## Bug fixes

- Now will reset the cache when items cannot be loaded. In the past
  error messages like `function
  'rx_0ba247452048de33b1ffb8af516714fc__calc_lhs' not provided by
  package 'rx_0ba247452048de33b1ffb8af516714fc_'` would cause the
  estimation to stop.  Now `rxode2::rxClean()` is run when this occurs.
