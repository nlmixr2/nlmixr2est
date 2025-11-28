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
  print = 1,
  trace = 0,
  covMethod = c("linFim", "fim", "r,s", "r", "s", ""),
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
  type = c("nelder-mead", "newuoa"),
  powRange = 10,
  lambdaRange = 3,
  odeRecalcFactor = 10^(0.5),
  maxOdeRecalc = 5L,
  perSa = 0.75,
  perNoCor = 0.75,
  perFixOmega = 0.1,
  perFixResid = 0.1,
  compress = FALSE,
  rxControl = NULL,
  sigdig = NULL,
  sigdigTable = NULL,
  ci = 0.95,
  muRefCov = TRUE,
  muRefCovAlg = TRUE,
  handleUninformativeEtas = TRUE,
  iovXform = c("sd", "var", "logsd", "logvar"),
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

  The number it iterations that are completed before anything is printed
  to the console. By default, this is 1.

- trace:

  An integer indicating if you want to trace(1) the SAEM algorithm
  process. Useful for debugging, but not for typical fitting.

- covMethod:

  Method for calculating covariance. In this discussion, R is the
  Hessian matrix of the objective function. The S matrix is the sum of
  each individual's gradient cross-product (evaluated at the individual
  empirical Bayes estimates).

  "`linFim`" Use the Linearized Fisher Information Matrix to calculate
  the covariance.

  "`fim`" Use the SAEM-calculated Fisher Information Matrix to calculate
  the covariance.

  "`r,s`" Uses the sandwich matrix to calculate the covariance, that is:
  \\R^-1 \times S \times R^-1\\

  "`r`" Uses the Hessian matrix to calculate the covariance as \\2\times
  R^-1\\

  "`s`" Uses the crossproduct matrix to calculate the covariance as
  \\4\times S^-1\\

  "" Does not calculate the covariance step.

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

  specifies the type of additive plus proportional errors, the one where
  standard deviations add (combined1) or the type where the variances
  add (combined2).

  The combined1 error type can be described by the following equation:

  \$\$y = f + (a + b\times f^c) \times \varepsilon\$\$

  The combined2 error model can be described by the following equation:

  \$\$y = f + \sqrt{a^2 + b^2\times f^{2\times c}} \times
  \varepsilon\$\$

  Where:

  \- y represents the observed value

  \- f represents the predicted value

  \- a is the additive standard deviation

  \- b is the proportional/power standard deviation

  \- c is the power exponent (in the proportional case c=1)

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

  This is the transformation used on the diagonal of the \`iov\`. The
  possibilities are:

  - `sd` Estimate the IOV as the standard deviation for IOV

  - `var` Estimate the IOV as the variance for IOV.

  - `logsd` Estimate the IOV as the log(sd) instead of sd.

  - `logvar` Estimate the IOV as the log(var) instead of variance.

- ...:

  Other arguments to control SAEM.

## Value

List of options to be used in
[`nlmixr2`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.md)
fit for SAEM.

## See also

Other Estimation control:
[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md),
[`nlmixr2NlmeControl()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2NlmeControl.md)

## Author

Wenping Wang & Matthew L. Fidler
