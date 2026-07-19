# Control Options for SAEM

Control Options for SAEM

## Usage

``` r
saemControl(
  seed = 99,
  nBurn = 200,
  nEm = 300,
  nmc = 3,
  nu = c(2, 2, 2),
  print = 1L,
  trace = 0,
  covMethod = c("sa", "analytic", "linFim", "fim", "r,s", "r", "s", "imp", ""),
  covMethodDeferred = NA_character_,
  covFull = TRUE,
  nSaCov = 500L,
  calcTables = TRUE,
  logLik = FALSE,
  nnodesGq = 3,
  nsdGq = 1.6,
  optExpression = TRUE,
  literalFix = FALSE,
  adjObf = TRUE,
  sumProd = FALSE,
  addProp = c("combined2", "combined1"),
  tol = 1e-06,
  itmax = 30,
  type = c("newuoa", "nelder-mead"),
  powRange = 10,
  lambdaRange = 3,
  odeRecalcFactor = 10^(0.5),
  maxOdeRecalc = 5L,
  indTolRelax = TRUE,
  perSa = 0.75,
  perNoCor = 0.75,
  perFixOmega = 0.1,
  perFixResid = 0.1,
  compress = TRUE,
  rxControl = NULL,
  sigdig = NULL,
  sigdigTable = NULL,
  ci = 0.95,
  muRefCov = TRUE,
  muRefCovAlg = TRUE,
  handleUninformativeEtas = TRUE,
  iovXform = c("sd", "var", "logsd", "logvar"),
  boundedTransform = TRUE,
  eventSens = c("jump", "fd"),
  mixProbMethod = c("regress", "regularized", "annealed"),
  mixProbStepExp = 1,
  mixProbPriorN = 20,
  mixSampleMethod = c("parallel", "msaem"),
  nonMuTheta = c("regress", "eta"),
  residWarmStart = TRUE,
  censOption = c("gauss", "laplace"),
  fast = FALSE,
  fastKernel = c("firstN", "throughout", "additive"),
  fastCov = c("auto", "jacobian", "hessian"),
  fastIter = 20L,
  fastLik = c("focei", "foce", "focep"),
  lbfgsLmm = 5L,
  lbfgsFactr = NULL,
  lbfgsPgtol = NULL,
  lbfgsMaxIter = 20L,
  nRetry = 10L,
  ...
)
```

## Arguments

- seed:

  Random Seed for SAEM step. (Needs to be set for reproducibility.) By
  default this is 99.

- nBurn:

  Number of iterations in the first phase, ie the MCMC/Stochastic
  Approximation steps. This is equivalent to Monolix's `K_0` or `K_b`.

- nEm:

  Number of iterations in the Expectation-Maximization (EM) Step. This
  is equivalent to Monolix's `K_1`.

- nmc:

  Number of Markov Chains. By default this is 3. When you increase the
  number of chains the numerical integration by MC method will be more
  accurate at the cost of more computation. In Monolix this is
  equivalent to `L`.

- nu:

  This is a vector of 3 integers. They represent the numbers of
  transitions of the three different kernels used in the
  Hasting-Metropolis algorithm. The default value is `c(2,2,2)`,
  representing 40 for each transition initially (each value is
  multiplied by 20).

  The first value represents the initial number of multi-variate Gibbs
  samples are taken from a normal distribution.

  The second value represents the number of uni-variate, or multi-
  dimensional random walk Gibbs samples are taken.

  The third value represents the number of bootstrap/reshuffling or
  uni-dimensional random samples are taken.

- print:

  Either a scalar print-frequency (\`0\` = suppress, \`1\` (default) =
  every evaluation, \`N\` = every Nth), OR a pre-built
  \[iterPrintControl()\] object. Equivalent to \`iterPrintControl(every
  = print, ncol = printNcol, useColor = useColor)\`.

- trace:

  An integer indicating if you want to trace(1) the SAEM algorithm
  process. Useful for debugging, but not for typical fitting.

- covMethod:

  Method for calculating covariance. In this discussion, R is the
  Hessian matrix of the objective function. The S matrix is the sum of
  each individual's gradient cross-product (evaluated at the individual
  empirical Bayes estimates).

  "`sa`" (default) Use the stochastic-approximation Fisher Information
  Matrix. After estimation, a dedicated covariance phase (`nSaCov`
  iterations) holds the parameters at the converged estimate and keeps
  resimulating the individual parameters, Monte-Carlo averaging the
  Louis observed-information integrand into a converged FIM decoupled
  from the cooling schedule (the approach used by Monolix; Kuhn &
  Lavielle 2005). Always includes every estimated population parameter
  (theta, the `Omega` diagonal variances, and residual).

  "`analytic`" Compute the FOCEI analytic observed-information
  covariance at the converged SAEM estimates. When the model is out of
  analytic-covariance scope (e.g. `linCmt()`, a non-normal likelihood,
  or a non-SD IOV parameterization) or the result is not positive
  definite, it falls back to the linearized Fisher information
  (`linFim`) with a message.

  "`linFim`" Use the Linearized Fisher Information Matrix to calculate
  the covariance.

  "`fim`" Use the Fisher Information Matrix accumulated during SAEM
  estimation to calculate the covariance. Like `sa` it inverts the
  observed information to a full theta + `Omega` diagonal + residual
  covariance, but uses the (noisier) estimation-phase matrix rather than
  a dedicated cov phase.

  For both `fim` and `sa` the simulation-based Fisher information covers
  the structural theta, the `Omega` diagonal variances, and additive
  residual error. Off-diagonal `Omega` covariances and
  proportional/combined residual error are not estimated reliably by the
  simulation FIM (the complete-data correction is unstable when
  between-subject variability dominates the residual), so those
  variance-block standard errors are spliced in from the linearized FIM
  (`linFim`).

  "`r,s`" Uses the sandwich matrix to calculate the covariance, that is:
  \\R^-1 \times S \times R^-1\\

  "`r`" Uses the Hessian matrix to calculate the covariance as \\2\times
  R^-1\\

  "`s`" Uses the crossproduct matrix to calculate the covariance as
  \\4\times S^-1\\

  "" Does not calculate the covariance step.

- covMethodDeferred:

  Internal. When a foreign covariance ("sa"/"imp") is requested, it is
  stashed here and computed post-fit at the converged estimates by the
  decoupled recompute engine
  ([`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.md)
  uses the same path); `NA` otherwise.

- covFull:

  Boolean (default `TRUE`) indicating the covariance should include
  every estimated population parameter – the structural and residual
  thetas plus the `Omega` variance/covariance elements – named
  `om.<eta>` / `cov.<eta>.<eta>`. When `FALSE` the legacy
  structural-theta-only covariance is reported. Ignored by
  `covMethod="sa"`, which is always full.

- nSaCov:

  Number of iterations in the dedicated stochastic-approximation
  covariance phase used by `covMethod="sa"` (default `500`). These
  iterations run at the converged estimate (parameters frozen) and only
  resimulate the individual parameters to build the observed Fisher
  information; a larger value gives a less noisy covariance. Ignored by
  other covariance methods.

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- logLik:

  boolean indicating that log-likelihood should be calculate by Gaussian
  quadrature.

- nnodesGq:

  number of nodes to use for the Gaussian quadrature when computing the
  likelihood with this method (defaults to 1, equivalent to the
  Laplacian likelihood)

- nsdGq:

  span (in SD) over which to integrate when computing the likelihood by
  Gaussian quadrature. Defaults to 3 (eg 3 times the SD)

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- addProp:

  Type of additive-plus-proportional error: \`"combined1"\`, where
  standard deviations add: \$\$y = f + (a + b\times f^c) \times
  \varepsilon\$\$; or \`"combined2"\`, where variances add: \$\$y = f +
  \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon\$\$. Here y =
  observed, f = predicted, a = additive sd, b = proportional/power sd, c
  = power exponent (1 in the proportional case).

- tol:

  This is the tolerance for the regression models used for complex
  residual errors (ie add+prop etc)

- itmax:

  This is the maximum number of iterations for the regression models
  used for complex residual errors. The number of iterations is
  itmax\*number of parameters

- type:

  indicates the type of optimization for the residuals; Can be one of
  c("nelder-mead", "newuoa")

- powRange:

  This indicates the range that powers can take for residual errors; By
  default this is 10 indicating the range is c(-10, 10)

- lambdaRange:

  This indicates the range that Box-Cox and Yeo-Johnson parameters are
  constrained to be; The default is 3 indicating the range c(-3,3)

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- indTolRelax:

  When \`TRUE\` (default), only subjects whose ODE solve produced
  NaN/Inf have their tolerances relaxed, and the relaxed tolerance
  persists across optimizer calls (sticky). When \`FALSE\`, all subjects
  have their tolerances relaxed on each retry and tolerances are reset
  afterward.

- perSa:

  This is the percent of the time the \`nBurn\` iterations in phase runs
  runs a simulated annealing.

- perNoCor:

  This is the percentage of the MCMC phase of the SAEM algorithm where
  the variance/covariance matrix has no correlations. By default this is
  0.75 or 75 Monte-carlo iteration.

- perFixOmega:

  This is the percentage of the \`nBurn\` phase where the omega values
  are unfixed to allow better exploration of the likelihood surface.
  After this time, the omegas are fixed during optimization.

- perFixResid:

  This is the percentage of the \`nBurn\` phase where the residual
  components are unfixed to allow better exploration of the likelihood
  surface.

- compress:

  Should the object have compressed items

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- sigdig:

  Specifies the "significant digits" that the ode solving requests. When
  specified this controls the relative and absolute tolerances of the
  ODE solvers. By default the tolerance is `0.5*10^(-sigdig-2)` for
  regular ODEs. For the sensitivity equations the default is
  `0.5*10\^(-sigdig-1.5)` (sensitivity changes only applicable for
  liblsoda). This also controls the `atol`/`rtol` of the steady state
  solutions. The `ssAtol`/`ssRtol` is `0.5*10\^(-sigdig)` and for the
  sensitivities `0.5*10\^(-sigdig+0.625)`. By default this is
  unspecified (`NULL`) and uses the standard `atol`/`rtol`.

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- muRefCov:

  This controls if mu-referenced covariates in \`saem\` are handled
  differently than non mu-referenced covariates. When \`TRUE\`,
  mu-referenced covariates have special handling. When \`FALSE\`
  mu-referenced covariates are treated the same as any other input
  parameter.

- muRefCovAlg:

  This controls if algebraic expressions that can be mu-referenced are
  treated as mu-referenced covariates by:

  1\. Creating a internal data-variable \`nlmixrMuDerCov#\` for each
  algebraic mu-referenced expression

  2\. Change the algebraic expression to \`nlmixrMuDerCov# \*
  mu_cov_theta\`

  3\. Use the internal mu-referenced covariate for saem

  4\. After optimization is completed, replace \`model()\` with old
  \`model()\` expression

  5\. Remove \`nlmixrMuDerCov#\` from nlmix2 output

  In general, these covariates should be more accurate since it changes
  the system to a linear compartment model. Therefore, by default this
  is \`TRUE\`.

- handleUninformativeEtas:

  boolean that tells nlmixr2's saem to calculate uninformative etas and
  handle them specially (default is \`TRUE\`).

- iovXform:

  Transformation used on the diagonal of the IOV: one of `"sd"`,
  `"var"`, `"logsd"`, or `"logvar"`.

- boundedTransform:

  When \`TRUE\` (default), bounded parameters are transformed for
  unbounded optimization methods and back-transformed for final
  estimates. \`FALSE\` optimizes on the original scale with bounds
  passed to the optimizer. \`NA\` transforms for optimization but skips
  the final back-transform.

- eventSens:

  Controls how dosing/event-parameter (\`alag\`, \`F\`, \`rate\`,
  \`dur\`) sensitivities are computed for THETA/ETA gradients:
  \`"jump"\` (default) uses rxode2's analytic event sensitivities;
  \`"fd"\` uses the legacy finite-difference behavior.

- mixProbMethod:

  For mixture models (\`mix()\`, more than one component), stabilizes
  the mixing-probability estimate against collapsing onto a single
  component (the responsibility used to update it is itself weighted by
  the current mixing probability, which can create a runaway feedback
  loop). Three options:

  \* \`"regress"\` (default): treat per-subject mixture membership as a
  fixed regressor. Each subject is hard-classified to a component up
  front, held fixed, and fed into the solve (via the existing
  mixture-index regressor), skipping the per-iteration soft-EM
  responsibility step entirely. Avoids the responsibility feedback loop
  / collapse by construction and is lower-bias; on heavily overlapping
  components it is higher-variance (an early misclassification is not
  revisited), so prefer \`"regularized"\` when membership is genuinely
  uncertain.

  \* \`"regularized"\`: blend \`mixProbPriorN\` pseudo-subjects,
  distributed per the initial mixing probability, into the
  responsibility average each iteration (Dirichlet/MAP-EM-style).
  Prevents collapse even in difficult cases, at the cost of some bias
  toward the initial guess; may need larger \`nBurn\`/\`nEm\`.

  \* \`"annealed"\`: give the mixing-probability update its own decaying
  step-size schedule (\`mixProbStepExp\`) instead of the
  full-replacement step used during \`nBurn\`. Lower bias, but does not
  by itself fix a systematic (non-noise-driven) collapse.

- mixProbStepExp:

  Only used when \`mixProbMethod="annealed"\`. Decay exponent for the
  mixing-probability step size (\`1/iteration^mixProbStepExp\`), applied
  from iteration 1. Default 1; smaller values decay more slowly.

- mixProbPriorN:

  Only used when \`mixProbMethod="regularized"\`. Number of
  pseudo-subjects blended into the responsibility average each
  iteration. Larger values are more robust to collapse but bias the
  estimate more and need more \`nBurn\`/\`nEm\`. Default 20.

- mixSampleMethod:

  For mixture models with per-component etas (split-ETA, e.g. \`cl \<-
  mix(tcl1 + eta.cl1, p1, tcl2 + eta.cl2)\`), controls the
  MCMC/sufficient-statistic architecture for the individual random
  effects, independent of \`mixProbMethod\`. BSV (\`\$omega\`) for split
  components is unreliable under \`"parallel"\` regardless of
  \`mixProbMethod\`.

  \* \`"parallel"\` (default): one full MCMC chain per component per
  subject per iteration, blended post hoc by responsibility. Mirrors
  NONMEM's \`\$MIX\` and correctly estimates BSV shared across
  components, but cannot cleanly separate per-component BSV for
  split-ETA models (each "wrong-hypothesis" chain still explores its
  non-owned column(s) as unconstrained prior noise).

  \* \`"msaem"\` (experimental): the MSAEM algorithm (Lavielle &
  Mbogning 2014), as used by Monolix. Simulates one random-effects
  trajectory per subject per iteration (label marginalized out via a
  closed-form responsibility) instead of parallel per-component chains,
  so no post-hoc blending is needed. Not compute-matched to
  \`"parallel"\` at equal \`nmc\` – set \`nmc\` to roughly \`nMix\`
  times its default for a fair comparison. Uses a model-aware stratified
  initialization for split-ETA components that reliably achieves full
  theta/fixed-effect separation. Split-ETA BSV recovery is improved (two
  numerical bugs fixed: an \`IGamma2_phi1\` blowup that locked variance
  to exactly zero, and an inverted responsibility sign) but still not
  reliable – it often settles at a safety-floor value rather than the
  true variance. Prefer \`"parallel"\` unless specifically evaluating
  this method.

- nonMuTheta:

  Controls how a population \`theta\` that is not mu-referenced (does
  not appear linearly with an eta – the SAEM \`phi0\` fixed effects) is
  estimated.

  \* \`"regress"\` (default): keep the parameter as a plain
  directly-estimated \`theta\` regressor. Each iteration \`phi0\` is
  estimated by a bounded direct optimization of the observation
  likelihood (robust coordinate descent within a local trust region,
  honoring the \`ini\` bounds), held fixed rather than drawn
  stochastically with a shrinking variance. This recovers population
  parameters that have no associated random effect more accurately, at
  some extra runtime (the objective re-solves the ODE).

  \* \`"eta"\`: the historic SAEM treatment (the parameter is carried
  through the stochastic \`phi0\` block).

- residWarmStart:

  Boolean (default \`TRUE\`); warm-start the residual-error parameters
  from the observed per-endpoint moments at the initial predictions
  (additive SD from \`sqrt(mean(err^2))\`, proportional SD from
  \`sqrt(mean((err/f)^2))\`), the same moment estimate \`est="npag"\`
  uses. Gives the stochastic step a better starting residual scale. Set
  \`FALSE\` to start from the \`ini\`-block residual values instead.

- censOption:

  Treatment of the second derivative for censored (M2/M3/M4/BLQ)
  observations in the FOCEI family. `"gauss"` (the default) keeps the
  historic uncensored Gauss-Newton curvature, matching common PMx tools;
  `"laplace"` uses the exact censored second derivative of the objective
  (a proper Laplace inner Hessian and analytic covariance). Accepted by
  `saemControl`/`nlmControl` for a uniform interface but inert there –
  SAEM (stochastic EM) has no Laplace inner Hessian, and NLM uses a
  finite-difference Hessian that already reflects censoring exactly.

- fast:

  Boolean enabling the fast-SAEM (f-SAEM) simulation step (Karimi,
  Lavielle and Moulines 2020). When \`TRUE\`, the MCMC simulation of the
  individual random effects uses an independent Metropolis-Hastings
  proposal centered at each subject's conditional MAP estimate with a
  Laplace/linearization covariance, which converges in far fewer SAEM
  iterations than the default random-walk Metropolis. The
  \`est="fsaem"\` method is sugar for \`saemControl(fast=TRUE)\`. By
  default this is \`FALSE\` (standard SAEM). The \`fast\*\` options
  below are only consulted when \`fast=TRUE\`.

- fastKernel:

  Schedule for the f-SAEM independent Metropolis-Hastings (IMH) kernel:

  \* \`"firstN"\` (default): use the IMH kernel for the first
  \`fastIter\` iterations, then revert to the standard random-walk
  kernels. This is the recipe used in the f-SAEM paper – the early
  iterations only need an approximate posterior, so the fast kernel
  accelerates the initial convergence and the steady-state behavior is
  unchanged.

  \* \`"throughout"\`: use the IMH kernel on every iteration for the
  whole run. Simpler, but recomputing the MAP/covariance every iteration
  is costlier and unnecessary near convergence.

  \* \`"additive"\`: append the IMH kernel alongside the standard
  random-walk kernels on every iteration. Most mixing, most cost.

- fastCov:

  Covariance used for the IMH Gaussian proposal:

  \* \`"auto"\` (default): Jacobian linearization for continuous-data
  endpoints, Hessian (Laplace) for non-continuous endpoints.

  \* \`"jacobian"\`: \`Gamma_i = (J' Sigma^-1 J + Omega^-1)^-1\` from
  the structural-model Jacobian at the MAP (continuous data only).

  \* \`"hessian"\`: \`Gamma_i = (-H + Omega^-1)^-1\` from the Hessian of
  the individual log-likelihood at the MAP (any data type).

- fastIter:

  Integer number of initial iterations to run the IMH kernel when
  \`fastKernel="firstN"\` (default 20). Ignored by the other schedules.

- fastLik:

  Inner likelihood used for the Hessian proposal path, one of
  \`"focei"\` (default), \`"foce"\` or \`"focep"\`. Selects which
  FOCEI-family individual likelihood is reused to build the proposal
  (and, when the Hessian path is active, reported by SAEM).

- lbfgsLmm:

  Integer number of BFGS corrections (the L-BFGS-B \`lmm\` memory) used
  when refining the fixed-effect-only parameters of a general
  log-likelihood model (\`ll(name) ~ \<expr\>\`) by direct L-BFGS-B
  optimization of the observation likelihood. Default 5.

- lbfgsFactr:

  Convergence tolerance on the relative reduction in the objective for
  that L-BFGS-B refinement (the \`factr\` control, in units of machine
  epsilon). When \`NULL\` (default) it is derived from \`sigdig\` the
  same way as \`foceiControl()\` (\`10^(-sigdig - 1) /
  .Machine\$double.eps\`).

- lbfgsPgtol:

  Convergence tolerance on the projected gradient for that L-BFGS-B
  refinement (the \`pgtol\` control). When \`NULL\` (default) it is
  derived from \`sigdig\` (\`10^(-sigdig - 1)\`).

- lbfgsMaxIter:

  Integer maximum number of iterations for that L-BFGS-B refinement.
  Default 20.

- nRetry:

  Integer number of times a bounded log-likelihood parameter's f-SAEM
  IMH proposal is re-drawn when it lands outside the parameter's bounds
  before being clamped to the violated boundary. Default 10.

- ...:

  Other arguments to control SAEM.

## Value

List of options to be used in
[`nlmixr2`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md)
fit for SAEM.

## References

Kuhn E, Lavielle M (2005). "Maximum likelihood estimation in nonlinear
mixed effects models." Computational Statistics & Data Analysis, 49(4),
1020-1038.
[doi:10.1016/j.csda.2004.07.002](https://doi.org/10.1016/j.csda.2004.07.002)

Jiang L, Roy A, Balasubramanian K, Davis D, Drusvyatskiy D, Na S (2025).
"Online Covariance Estimation in Nonsmooth Stochastic Approximation."
arXiv:2502.05305.
[doi:10.48550/arXiv.2502.05305](https://doi.org/10.48550/arXiv.2502.05305)

## See also

Other Estimation control:
[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md),
[`nlmixr2NlmeControl()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2NlmeControl.md)

## Author

Wenping Wang & Matthew L. Fidler
