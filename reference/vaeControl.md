# Control for vae (variational autoencoder) estimation method in nlmixr2

Variational-autoencoder NLME estimation (Rohleff et al., CPT:PSP 2025):
an LSTM encoder learns the individual posterior q(eta\|y) and an rxode2
decoder reconstructs the observations, trained on an ELBO / BICc-ELBO
objective for simultaneous population-parameter estimation and covariate
selection.

## Usage

``` r
vaeControl(
  seed = 42L,
  itersBurnIn = 100L,
  klWarmup = 50L,
  gammaIter = 250L,
  iters = 300L,
  nGradStep = 5L,
  hiddenDim = 25L,
  learningRate = 0.005,
  burnInLearningRate = 0.008,
  sigma0 = NULL,
  covariateSelection = TRUE,
  covSelectAlpha = 2,
  bnbStrategy = c("lifo", "fifo", "lc"),
  nonMuTheta = c("regress", "eta", "fix", "none"),
  nonMuEtaOmega = 0.01,
  likelihood = c("focei", "foce", "focep", "laplace"),
  objf = c("importanceSampling", "linear"),
  nIsSample = 3000L,
  returnVae = FALSE,
  print = 1L,
  useColor = NULL,
  printNcol = NULL,
  covMethod = c("analytic", "r,s", "r", "s", ""),
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = FALSE,
  adjObf = TRUE,
  ci = 0.95,
  sigdig = NULL,
  sigdigTable = NULL,
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  indTolRelax = TRUE,
  eventSens = c("jump", "fd"),
  rxControl = NULL,
  ...
)
```

## Arguments

- seed:

  Random seed for the VAE training (encoder init, Adam,
  reparameterization sampling); default 42. Training is stochastic, so a
  fixed seed makes every fit reproducible.

- itersBurnIn:

  Number of burn-in iterations (encoder-only, tiny KL weight) before the
  main EM phase.

- klWarmup:

  Number of KL-annealing iterations over which the KL weight is ramped
  from a small value to 1 (prevents posterior collapse).

- gammaIter:

  Number of main iterations before the EMA-smoothing phase of the
  population-parameter update begins.

- iters:

  Total number of main-loop iterations (after burn-in).

- nGradStep:

  Number of Adam gradient steps per EM outer iteration (the reference
  \`L_iter\`).

- hiddenDim:

  LSTM hidden dimension (the reference \`h_dim\`).

- learningRate:

  Adam learning rate used in the main training phase.

- burnInLearningRate:

  Adam learning rate used during burn-in.

- sigma0:

  Encoder prior standard deviation(s) at initialization (a small value
  giving a sharp initial posterior). \`NULL\` uses a small default per
  individual parameter. This is distinct from the \`ini()\` omega.

- covariateSelection:

  When \`TRUE\` (default) perform automated BICc-ELBO covariate
  selection during training; when \`FALSE\` fit the given fixed
  covariate structure only (faster population-only mode).

- covSelectAlpha:

  Starting multiplier for the covariate-selection L0 penalty, ramped
  linearly from \`covSelectAlpha\` down to \`1\` over the \`klWarmup\`
  warmup iterations and held at \`1\` afterward (matching the reference
  implementation's \`linspace(alpha, 1, kl_iter)\`). Values \`\> 1\`
  penalize covariate entry more heavily early in training; \`1\`
  disables the ramp.

- bnbStrategy:

  Frontier discipline for the exact branch-and-bound covariate
  selection: \`"lifo"\` (default, last-in-first-out depth-first search),
  \`"fifo"\` (first-in-first-out) or \`"lc"\` (least cost / best-first).
  The solver is exact, so the selected covariates are identical for
  every strategy; only the search order (and thus efficiency) differs.

- nonMuTheta:

  How to treat a structural population \`theta\` that has no random
  effect (is not mu-referenced) so it can still be estimated by the VAE
  (which only estimates parameters that occupy the latent space). For
  the eta-injection modes a small eta is injected so the parameter
  enters the latent space, and the reported fixed effect is \`theta +
  mean(eta)\` with the temporary eta dropped from the output model.

  \* \`"regress"\` (default, matching \`saemControl(nonMuTheta=)\`): no
  eta is injected; instead each such theta is estimated directly,
  re-optimized every M-step by a bounded \`bobyqa\` regression against
  the FOCEi inner likelihood (bounds from the \`ini()\` lower/upper),
  blended with the M-step gain. This recovers a no-random-effect
  population parameter without adding a spurious random effect.
  \`nonMuEtaOmega\` is unused in this mode. \* \`"eta"\`: inject the eta
  with an ESTIMATED omega (starting at \`nonMuEtaOmega\`); the typical
  value is estimated and appears in the iteration table. \* \`"fix"\`:
  inject the eta with omega held FIXED at \`nonMuEtaOmega\` AND hold the
  typical-value theta fixed at its \`ini()\` value. Nothing about the
  parameter is estimated, so it is not shown in the iteration table (it
  is reported at its \`ini()\` value, marked fixed, with the injected
  eta dropped). \* \`"none"\`: leave non-mu-referenced thetas frozen at
  their \`ini()\` value (the historic behavior).

- nonMuEtaOmega:

  Variance of the eta injected for a non-mu-referenced theta (starting
  value for \`nonMuTheta="eta"\`, fixed value for \`nonMuTheta="fix"\`;
  unused for \`"regress"\`).

- likelihood:

  Inner likelihood used for the objective, EBEs, and gradients, all run
  through the same FOCEi inner interface: \`"focei"\` (default, with
  eta-epsilon interaction), \`"foce"\` (no interaction, NONMEM FOCE with
  R frozen at the population prediction), \`"focep"\` (FOCE+, no
  interaction but R evaluated at the live conditional eta), or
  \`"laplace"\`.

- objf:

  Which objective-function value is active for AIC/BIC/BICc. Both the
  linearization and importance-sampling -2LL are always computed and
  stored; this selects the default active one.

- nIsSample:

  Number of importance-sampling draws for the IS -2LL.

- returnVae:

  When \`TRUE\` return the raw VAE training object instead of the
  nlmixr2 fit.

- print:

  Either a scalar print-frequency (\`0\` = suppress, \`1\` (default) =
  every evaluation, \`N\` = every Nth), OR a pre-built
  \[iterPrintControl()\] object. Equivalent to \`iterPrintControl(every
  = print, ncol = printNcol, useColor = useColor)\`.

- useColor:

  Logical (or \`NULL\`) emit ANSI bold/color escapes in the iteration
  print. \`NULL\` (default) defers to \[crayon::has_color()\].

- printNcol:

  Integer (or \`NULL\`) parameter columns per row before wrapping.
  \`NULL\` (default) uses \`floor((getOption("width") - 23) / 12)\`.

- covMethod:

  Method for calculating the covariance at the VAE estimates, run
  through the FOCEi covariance step; the same choices as
  [`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md):
  `"analytic"` (default), `"r,s"`, `"r"`, `"s"`, or `""` to skip.

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- literalFixRes:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- addProp:

  Type of additive-plus-proportional error: \`"combined1"\`, where
  standard deviations add: \$\$y = f + (a + b\times f^c) \times
  \varepsilon\$\$; or \`"combined2"\`, where variances add: \$\$y = f +
  \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon\$\$. Here y =
  observed, f = predicted, a = additive sd, b = proportional/power sd, c
  = power exponent (1 in the proportional case).

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

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

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

- indTolRelax:

  When \`TRUE\` (default), only subjects whose ODE solve produced
  NaN/Inf have their tolerances relaxed, and the relaxed tolerance
  persists across optimizer calls (sticky). When \`FALSE\`, all subjects
  have their tolerances relaxed on each retry and tolerances are reset
  afterward.

- eventSens:

  Controls how dosing/event-parameter (\`alag\`, \`F\`, \`rate\`,
  \`dur\`) sensitivities are computed for THETA/ETA gradients:
  \`"jump"\` (default) uses rxode2's analytic event sensitivities;
  \`"fd"\` uses the legacy finite-difference behavior.

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- ...:

  Other arguments to control SAEM.

## Value

vae control structure (class \`vaeControl\`)

## Details

Covariate selection – MIQP vs. branch-and-bound. Per latent parameter
the selection step minimizes the same L0/BIC objective \`RSS_S/omega +
log(N)\*\|S\|\` over subsets \`S\` of the candidate covariates
(\`RSS_S\` is the residual sum of squares of the ordinary-least-squares
fit on the intercept plus \`S\`). The reference implementation (Rohleff
et al.) writes this as a Mixed-Integer Quadratic Program (MIQP) – binary
include/exclude indicators with big-M constraints – and solves it with
the commercial Gurobi solver through \`cvxpy\`. No MIQP-capable solver
is freely available in R: Gurobi is commercial/licensed, and the open QP
solvers on CRAN (e.g. \`osqp\`) are continuous-only and cannot represent
the binary selection. A continuous convex relaxation (L1 / lasso) would
be solvable but only approximates best subset.

This package instead solves the identical L0/BIC objective EXACTLY with
a self-contained branch-and-bound: each candidate support's coefficients
are the closed-form OLS fit and branches are pruned by a valid lower
bound (the RSS of the OLS fit using all still-free covariates). It
therefore returns the same optimum the MIQP would – no commercial
dependency and no relaxation/accuracy loss – and scales to a few dozen
covariates. The search is worst-case exponential in the number of
covariates, but the pruning makes the practical (sparse) case fast (e.g.
32 candidate covariates in a fraction of a second).

## Author

Matthew L. Fidler
