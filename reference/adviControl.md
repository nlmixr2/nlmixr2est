# Control for advi (automatic differentiation variational inference) in nlmixr2

Variational-inference NLME estimation following Kucukelbir et al.
(2017): the latent variables are transformed to an unconstrained real
coordinate space, a Gaussian variational family is posited there, and
the ELBO is maximized by stochastic gradient ascent using the
reparameterization trick. The gradient of the log-joint is obtained from
the FOCEi forward sensitivities (inner per-subject eta gradient and the
outer population sensitivity contraction), not from automatic
differentiation. The whole optimization loop runs in C++.

## Usage

``` r
adviControl(
  seed = 42L,
  iters = 300L,
  nMc = 1L,
  adviFamily = c("fullRank", "meanField"),
  pointEstimate = TRUE,
  optim = c("advi", "adam"),
  adaptEta = TRUE,
  etaCandidates = c(0.01, 0.025, 0.05, 0.1, 0.25),
  tau = 1,
  alpha = 0.1,
  tol = 1e-04,
  likelihood = c("focei", "foce", "focep", "laplace"),
  returnAdvi = FALSE,
  resume = NULL,
  print = 1L,
  useColor = NULL,
  printNcol = NULL,
  covMethod = c("advi", "analytic", "r,s", "r", "s", ""),
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

  Random seed for the ADVI optimization (reparameterization sampling);
  default 42. The Monte-Carlo gradient is stochastic, so a fixed seed
  makes every fit reproducible. Reparameterization noise is drawn from a
  counter-based stream keyed by the global iteration index, so a shorter
  run is a bit-for-bit prefix of a longer one and results are
  independent of the number of cores.

- iters:

  Total number of ADVI (stochastic gradient ascent) iterations.

- nMc:

  Number of Monte-Carlo samples used to approximate the ELBO gradient at
  each iteration (the paper's \`M\`; typically 1-10).

- adviFamily:

  Variational family in the unconstrained space. \`"fullRank"\`
  (default) uses a block full-rank Gaussian: a dense \`neta x neta\`
  Cholesky factor per subject plus a dense block over the population
  vector (mean-field across blocks). \`"meanField"\` uses a fully
  factorized (diagonal) Gaussian. Mean-field is faster but is known to
  underestimate marginal variances.

- pointEstimate:

  When \`TRUE\` (default) run a variational-EM hybrid: the variational
  posterior covers the per-subject etas only, and the population
  parameters (thetas / omega / residual error) are point estimates
  maximized by the ADVI gradient (stochastic maximum likelihood); output
  semantics match FOCEi/SAEM. When \`FALSE\` run full Bayes: the
  variational posterior also covers the unconstrained population vector,
  with flat priors.

- optim:

  Stochastic optimizer. \`"advi"\` (default) uses the paper's adaptive
  step-size sequence (Eqs 10-11); \`"adam"\` uses Adam.

- adaptEta:

  When \`TRUE\` (default) adaptively choose the step-size scale \`eta\`
  by a short search over \`etaCandidates\` before the main loop; when
  \`FALSE\` use a fixed \`eta\` (the first \`etaCandidates\` entry).

- etaCandidates:

  Candidate step-size scales searched when \`adaptEta\` is \`TRUE\` (the
  paper searches \`c(0.01, 0.1, 1, 10, 100)\`).

- tau:

  Stabilizing constant \`tau \> 0\` in the step-size denominator (paper
  Eq 10); the step-size is insensitive to it.

- alpha:

  Weighting \`alpha\` in (0, 1) of new vs old gradient information in
  the step-size memory recursion (paper Eq 11).

- tol:

  Convergence tolerance on the relative change in the ELBO; the loop may
  stop early once the change stays below this. \`0\` disables early
  stopping (run all \`iters\`).

- likelihood:

  Inner likelihood used for the per-subject objective and gradient, run
  through the FOCEi inner interface: \`"focei"\` (default), \`"foce"\`,
  \`"focep"\`, or \`"laplace"\`.

- returnAdvi:

  When \`TRUE\` return the raw ADVI optimization object instead of the
  nlmixr2 fit.

- resume:

  Optional warm-resume state: a previous \`est="advi"\` fit (or its
  \`\$env\$adviState\`). The optimization continues from that state for
  \`iters\` more iterations, bit-for-bit identical to a single fresh run
  of the combined length (the counter-based RNG is keyed by the global
  iteration index).

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

advi control structure (class \`adviControl\`)

## Author

Matthew L. Fidler
